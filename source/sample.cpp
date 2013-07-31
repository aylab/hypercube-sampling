/*
Deterministic simulator for zebrafish segmentation
Copyright (C) 2013 Ahmet Ay, Jack Holland, Adriana Sperlea, Sebastian Sangervasi

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <cstring>
#include <unistd.h>
#include <sys/wait.h>
#include <cmath>
#include <sys/resource.h>
#include <fcntl.h>

#include "sample.hpp"
#include "ihs.hpp"
#include "file_io.hpp"
#include "analyze.hpp"

using namespace std;

/*
	This program is used for sampling values from a parameter space, sending them to a simulation program,
and intelligently processing the results it sends back to find parameter sets that are successful. The
program by default uses recursive Latin Hypercube Sampling. The particular implementation it uses for the
LHS sampling is Improved Hypercube Sampling as designed by John Burkhardt -- see ihs.cpp and ihs.hpp for
his code. The program performs repeated LHS, choosing random values from within the cubes specified by LHS
and pipes these sets to child processes that run zebrafish somitogenesis simulations. Based on the resulting
scores of each parameter set, the simulation chooses which cubes within the space were most successful
(depending on a user-given theshold) and recurses, performing further LHS within the successful cubes. The
program stops when either it reaches a user-given maximum recursion depth, or it finds a cube that has mostly 
successful values (based on a user-given threshold). It writes the final successful parameter sets to a ".params"
file and the successful cube ranges to a ".ranges" file.
	The code also has an alternate mode, with input "-S" which enables a simple sampling from the parameter 
space, i.e. choosing purely random values without LHS. Based on simulation results it then picks sets with scores
within 10% (or user specified value) and cuts the sample range by up to two standard deviations from the max
and min ranges of those sets. It repeats this process until the top 3% have perfect scores or until it reaches
a user-specified number of reruns.  
	The executable file that does the simulation is "sogen-deterministic/deterministic" for which the source code 
is located in "sogen-deterministic/source/". All sampling source code is in "sogen-deterministic/sampling/",
which compiles to the executable "sogen-deterministic/sampling/sampler".	
*/
int main(int argc, char** argv){
	//Initializing the program.
	set_resources(); // (AAy: ??? Setting resources for what?)
	input_params ip; 
	accept_params(argc, argv, ip); // (AAy: ??? Getting commandline arguments)

	// Getting the parameter ranges, as doubles, from a file.
	// The double array has the format: [min1, max1, min2, max2, ...]
	char* indata = NULL; 
	indata = read_file(ip.ranges_file); // (AAy: ??? Read the ranges file)
	// Keeps a count of how many ranges should be zeroed out.
	count_dims(indata, ip.dims);
	if (ip.dims[0] < 0){ // (AAy: ??? What is this looping doing?)
		cout << "There is a bracket mismatch!" << "\n" << "Please correct entry " << ip.dims[0]+ip.dims[1] << "\n";
		free(indata);
		return 1;
	}
	
	double ranges[ip.dims[0]*2]; 
	int range_map[ip.dims[0]];	//Used for mapping the generated ranges to the correct output format with chosen values zeroed out.
	parse(ip.dims[0], indata, ranges, range_map);
	free(indata); //Move this down if indata is needed again for some reason. 
	
	//This actually does the work of sampling and running the simulations.
	recurse_params rp(ip); 
	timestamp();
	//If the user specified they would just like some random parameter sets, just give them that.
	if(ip.raw_data > 0){ 
		make_raw(ip, rp, range_map,ranges);
		
	//Defalt behavior is to actually run the LHS code. 
	}else if(ip.lhs){
		range_list* valid_cubes = new range_list(); 
		sample_recurse(ip, rp, ranges, range_map, valid_cubes);
		
		// Writing the successful cubes (ranges):
		range_node* to_write = valid_cubes->head; 
		int append = 0;
		for(; to_write != NULL; to_write = to_write -> next){
			range_write (append, 2*ip.dims[0], ip.segments, ip.duplication, to_write->grade, to_write -> ranges, ip.cubes_file );
			append++;
		}
		cout << "Successful cubes found: " << append << "\n"; 
		delete valid_cubes;
		
	//If lhs was turned off, this will just run simple sampling.	
	} else{ 
		simple_sample(ip, rp, ranges, range_map); 
	}
	timestamp(); 

	if(ip.quiet) cout_switch(false, ip);
	if(rp.failure != NULL) usage(rp.failure, rp.failcode); 
	return 0;
}

/*
	This recursive function:
	 -gets IHS cubes,
	 -scales them to the current parameter ranges
	 -randomly selects values from within each cube to simulate
	 -sends the parameter sets to the simulation program
	 -analyzes the results an repeats accordingly.
	Successful cubes are added to the linked list valid_cubes.  
*/
bool sample_recurse(input_params& ip,  recurse_params& rp, double* ranges, int* range_map, range_list* valid_cubes){
	rp.depth++;
	//Base case: if the mutant threshold was not met for this cube, give up on this range.
	if(rp.depth > ip.max_depth){
		cout << "Reached maximum recursion depth. " << "\n";
		rp.depth--;
		return false; 
	}
	
	//Looping to get many points from within the space.
	int success_count = 0;
	//This list is used to hold on to all cube that were successful enough at this depth and should be recursed upon.
	range_list* new_cubes = new range_list();
	for(int rerun = 0; rerun < ip.lhs_runs; rerun++){
		//Getting the IHS choice of segments.
		ih_sample(ip, rp.segs);

		//Applying the sample values to our actual ranges.
		//This array also gets zeroes where specified (by range_map).
		scale_sample(ip, rp.segs, range_map, ranges, rp.samples);
		//Running simulation(s).
		simulate_samples(ip, rp); 
		if(rp.failure != NULL){
			rp.depth--;
			return false;
		}
		cout << "Received simulation results from deterministic. (Depth " << rp.depth << ", Run: " << rerun <<")" << "\n";
	
		//Analyzing the resulting grades 
		success_count += analyze(ip, rp, ranges, new_cubes); 
		if(ip.verbose){
			sample_write (rerun+rp.depth-1, ip.dims[0]+ip.dims[1], ip.segments, ip.duplication, ip.lhs_runs, rp.passing_grade, rp.samples, rp.grades, ip.verbose_file );	
		}
	}
	
	//If this granularity was successful enough, don't recurse deeper.
	double percent_success = (double)success_count / (double)(ip.segments*ip.lhs_runs); 
	if(percent_success >= ip.cube_threshold){
		cout<< "At depth "<< rp.depth << " this space had " << percent_success*100 << "% successful cubes." << "\n";
		if(rp.depth ==1) valid_cubes->concat(new_cubes);
		//sample_write (rp.good_samples, ip.dims[0]+ip.dims[1], ip.segments, ip.duplication, ip.lhs_runs, rp.passing_grade, rp.samples, rp.grades, ip.params_file );
		rp.good_samples++;
		rp.depth--;
		return true;
	}
	
	//Otherwise, recurse using more precise ranges.
	range_node* ay = new_cubes->head;
	range_node* next; 
	while( ay != NULL){ 
		next = ay->next;
		cout << "Increasing recursion depth." << "\n";
		if(sample_recurse(ip, rp, ay->ranges, range_map, valid_cubes)){
			valid_cubes->take(new_cubes, ay);
			//range_write (1, 2*ip.dims[0], ip.segments, ip.duplication, ay->grade, ay -> ranges, ip.cubes_file );
		}
		if(rp.failure != NULL) break;
		ay = next; 
	}
	delete new_cubes; 
	rp.depth--;
	return false;
}

/*Getting the (improved) latin hypercube sample
 This funciton takes:
	dim_num = number of parameters
	segments = number of sets that will be generated
	fout = file to write to.
It returns:
	segs = the IHSampled segment numbers.
And it writes segs to a file.
*/
void ih_sample( input_params& ip, int* segs){
	ihs ( ip.dims[0], ip.segments, ip.duplication, &ip.random_seed, segs ); 

	//Write out the unscaled segment values (less important)
	if(ip.write_segments) ihs_write (ip.dims[0], ip.segments, ip.duplication, ip.random_seed, ip.random_seed, segs, ip.segments_file);	
	
}

/* This funciton takes:
	ip = input parameters
	segments = number of sets that will be generated
	segs = the IHS random sample of segments in the parameter space
	ranges = the ranges
	samples = the array in which the resulting, properly-scaled values will be stored
*/
void scale_sample(input_params& ip, int* segs, int* map, double* ranges, double* samples){
	int param_num = ip.dims[0] + ip.dims[1];
	int i = 0;
	int j = 0;
	int mloc = 0;
	int zshift = 0;
	double min, max, low, high;
	// Cycles through every parameter set.
	for(; i < ip.segments; i++){
		//Cycles through every dimension/parameter in each set.
		for(; j < ip.dims[0]; j++){
			//Gets the appropriate range values.
			low = ranges[j*2];
			high = ranges[j*2 + 1];

			//Calculates the correct min/max edges of the hypercube that was selected by the segment sampling
			min = low + (high - low)*(double)(segs[i*ip.dims[0] + j] - 1) / (double)ip.segments;
			max = low + (high - low)*(double)segs[i*ip.dims[0] + j] / (double)ip.segments;
			
			while(map[mloc] > j+zshift){
				samples[i*(param_num) + j+zshift] = 0;
				zshift++; 
			}
			mloc++;
			//Chooses a random value from within the hypercube.
			samples[i*(param_num) + j+zshift] = rand_range(min, max); 
		}
		while(param_num > j+zshift){
			samples[i*(param_num) + j+zshift] = 0;
			zshift++;	
		}
		mloc = 0;
		zshift = 0;
		j=0;	
	}
}

/* A simple function that generates a random number from the given range.*/
double rand_range(double min, double max){ 
	double rand = randouble();
	return (max-min)*(rand)+min;
}

/*	Establishing a communication pipe from the parent (sampler) to each simulation child for the passing
of parameter sets and results. */
bool make_pipes(int processes, int** pipes){ 
	for(int i = 0; i < processes; i++){
		pipes[i] = new int[2];
		if (pipe(pipes[i]) == -1) {
			del_pipes(i+1, pipes, true);
			return false;
		}
	}
	return true;
}

/* 	Closes the reading end of each pipe and deletes the data needed to store them.
	Child processes (../deterministic) are responsible for closing the writing end, but this will try to close them in case of child failure.
*/
void del_pipes(int processes, int** pipes, bool close_write){ 
	for(int i = 0; i < processes && pipes[i] != NULL; i++){
		close(pipes[i][0]);
		if(close_write && fcntl(pipes[i][1], F_GETFD ) != -1) close(pipes[i][1]);
		delete[] pipes[i];
	}
} 

/*	Mallocs an array of char** that are needed for passing arguments to the execv call for each child simulation. 
It fills in the appropriate argument space with the file descriptor of the pipe the children will read/write with.*/
char*** make_args(input_params& ip, int** pipes){ 
	int strlen_num;
	int pipe_loc;
	char*** child_args = (char***)malloc(sizeof(char**)*ip.processes);
	for(int i = 0; i < ip.processes; i++){
		pipe_loc = 0; 			
		child_args[i] = (char**)malloc(sizeof(char*)*ip.sim_args_num);
		for(int j = 0; j < ip.sim_args_num; j++){
			if(j == 2 || j == 4){
				strlen_num = log10(pipes[i][pipe_loc]+1)+1;
				child_args[i][j]= (char*)malloc(sizeof(char)*(strlen_num+1));
				sprintf(child_args[i][j], "%d", pipes[i][pipe_loc]);
				pipe_loc++;
			} else if (ip.simulation_args[j] == NULL){
				child_args[i][j] = NULL;
			}else{
				child_args[i][j] = strdup((const char*) ip.simulation_args[j]);
			}
		}
	}
	return child_args;
}

/*	Deleting the arguments that were used.*/
void del_args(input_params& ip, char*** child_args){ 
	for(int i = 0; i < ip.processes; i++){
		for(int j = 0; j < ip.sim_args_num; j++){
			if(child_args[i][j] != NULL) 
				free(child_args[i][j]);	
		}
		free(child_args[i]);	
	}
	free(child_args);
}

/*	Determines how many parameter sets shoudl be passed to each child by distributing them evenly. If the
number of children does not divide the number of sets, the remainder r is distributed among the first r children.*/
void segs_per_sim(int segments, int processes, int* distribution){ // (AAy: ??? What is this code doing?)
	int remainder = segments % processes;
	for(int i = 0; i < processes; i++){
		distribution[i] = (remainder > 0 ? segments/processes + 1 : segments/processes);
		remainder--; 
	}
}

void print_samples(int start, int stop, int dim_num, double* samples){ // (AAy: ??? Printing the parameter sets that are used?)
	int rstart = dim_num*start;
	int rstop = dim_num*stop+rstart;
	for(int i = 0; i+rstart < rstop; i++){
		cout << samples[i+rstart-1];
		if(i%dim_num == 0){
			cout << "\n";
		}else{
			cout << ",";
		}
	}
}
/*
	This function takes care of running the simulation by forking and executing (execv)
by calling ../deterministic. The parameter sets are passed to child processes and the results of
the simulations are passed back via a read/write pipe pir for each child.*/
void simulate_samples(input_params& ip, recurse_params& rp ){ 	
	int param_num = ip.dims[0] + ip.dims[1];
	
    int* pipes[ip.processes];
    if(!make_pipes(ip.processes, pipes)){ 
    	rp.failure = strdup("!!! Failure: could not pipe !!!\n"); 
    	return;
    }
    char*** child_args = make_args(ip, pipes); 
    
    pid_t simpids[ip.processes]; 
    for(int i = 0; i < ip.processes; i++){
    	simpids[i] = fork();
		if (simpids[i] == -1) {
			rp.failure = strdup("!!! Failure: could not fork !!!\n");
   			break;
        }
        //Child runs simulation.  
		if (simpids[i] == 0) {  
			if (-1 == execv(ip.sim_exec, child_args[i])){ 
				rp.failure = strdup("!!! Failure: could not exec deterministic !!!\n");
				break; 
			}
		}   	
    }
    if(rp.failure != NULL){
    	del_args(ip,child_args); 
		del_pipes(ip.processes, pipes, true); 
		return;  
    }

    // Parent gives sets and processes results. Writes params to simpipe[1], reads results from simpipe[0].
    void* param_str; //Used as a writing buffer.
    void* sets_str;
    
    //Send the number of parameters that are used for each simulation (should be 95 as of 6/21/13).
	param_str = malloc(sizeof(int));
	sets_str = malloc(sizeof(int)); 
	memcpy(param_str, &param_num, sizeof(int));
	
	//This array is used for distributing the number of simulations sent to each child process.
	int distribution[ip.processes]; 
	segs_per_sim(ip.segments, ip.processes, distribution);
	
	int samples_sent = 0;
    for(int i = 0; i < ip.processes; i++){
		write(pipes[i][1], param_str, sizeof(int));
	 
		//Send number of sets that will be sent.
		memcpy(sets_str,  distribution + i, sizeof(int)); 
		write(pipes[i][1], sets_str, sizeof(int)); 
	
		//Send the array of parameter values.
		write(pipes[i][1], rp.samples + samples_sent, sizeof(double)*param_num*distribution[i]);
		samples_sent += param_num*distribution[i];
	}
	free(param_str);
	free(sets_str);
	
	//Loop for waiting on children and checking their exit status.		
    ssize_t result; //For reading purposes.
	samples_sent = 0;
    for(int i = 0; i < ip.processes; i++){ 
		int status = 0; 
		waitpid(simpids[i], &status, WUNTRACED);

		if(WIFEXITED(status)){
			cout << "Child (" << simpids[i] << ") exited properly with status: " << WEXITSTATUS(status) << "\n";	
		} else if (rp.failure == NULL){
			if(WIFSIGNALED(status)){
				rp.failcode = WTERMSIG(status);
				#ifdef WCOREDUMP // (AAy: ??? What is this code doing?)
				if(WCOREDUMP(status)){
					rp.failure = strdup("!!! Failure: child experienced core dump !!!");
				} else{
					rp.failure = strdup("!!! Failure: child received signal stop !!!");
				}
				#else
					rp.failure = strdup("!!! Failure: child received some kind of signal. Exact status is uncertain because your OS does not support core dump reporting.\n");
				#endif
			} else if(WIFSTOPPED(status)){
				rp.failure = strdup("!!! Failure: child stopped by delivery of a signal !!!");
				rp.failcode = WSTOPSIG(status);
			} else{
				rp.failure = strdup("!!! Failure: child process did not exit properly !!!");
				rp.failcode = simpids[i];
			}
		}
		if(rp.failure != NULL) print_samples(samples_sent, distribution[i], param_num, rp.samples);
		samples_sent += distribution[i];
	}
	//Children are done, so we know we can delete the argument array.
	del_args(ip,child_args);
	if(rp.failure != NULL){
		del_pipes(ip.processes, pipes, true);   
		return ; 
	}
	
	int samples_got = 0;
	int passing_grade;
	for(int i = 0; i < ip.processes; i++){ 
		//Read the test results.
		result = read(pipes[i][0], &passing_grade, sizeof(int));
		if(result != sizeof(int)){
			rp.failure = strdup("!!! Failure: could not read passing_grade from deterministic !!!"); 
			rp.failcode = i;
			break;
		}
		
		result = read(pipes[i][0], rp.grades + samples_got, sizeof(int)*distribution[i]);
		if(result == -1 || (unsigned int) result !=  sizeof(int)*distribution[i]){
			rp.failure = strdup("!!! Failure: could not read simulation grades from deterministic !!!\n"); 
			rp.failcode = i;
			break;
		}
		
		samples_got += distribution[i];
	}
	rp.passing_grade = passing_grade;
	del_pipes(ip.processes, pipes, true); 
}

/*	This function takes the results of simulation and determines which cubes of the LHS should be sampled
from at the next level of recursion. By default the passing condition is ip.condition_threshold but the 
user can also turn on ip.take_best which will recurse upon the cube with the highest score even if it 
did not satisfy the threshold.*/
int analyze(input_params& ip, recurse_params& rp, double* ranges, range_list* new_cubes){
	int success = 0;
	int best_run = 0;
	bool desperate = true;
	double best_grade = 0;
	double grade;
	int* new_seg;
	double* new_range;
	for(int i = 0; i < ip.segments; i++){
		grade = (double)rp.grades[i] / (double)(rp.passing_grade);
		if(grade > best_grade){
			best_run = i;
			best_grade = grade;
		}
		cout << "Grade: " << rp.grades[i] << " ----- Passing grade: " << rp.passing_grade << "\n";
		//Check to see if this segment yeilded a successful parameter set.
		if( grade  >= ip.calc_threshold(rp.depth)){
			desperate = false; 
			new_range = new double[ip.dims[0]*2]; 
			new_seg = new int[ip.dims[0]];
			rescale(i*ip.dims[0], ip.dims[0], ip.segments, rp.segs, new_seg, ranges, new_range);
			//Check to see if the segment has already been discoverd and added to the list
			if(new_cubes->find(ip.segments, new_seg)){ 
				delete[] new_range;
				delete[] new_seg;
			} else{
				//If it's new, put it onto the list.
				new_cubes->add(new_seg, new_range, grade);
				//Write the successful parameter set to file.
				single_sample_write (rp.good_samples, ip.dims[0]+ip.dims[1],rp.passing_grade, rp.grades[i],  rp.samples+(i*(ip.dims[0]+ip.dims[1])),ip.params_file);
				rp.good_samples++;
			}
			success++;
		}
	}
	//This check is used if take_best is enabled -- in which case the program will recurse on the best segment found even if it did not pass the threshold.
	if(ip.take_best && desperate && best_grade>0){ 
		cout << "No successful cubes. Taking the best cube which satisfied " << best_grade*100 <<"% of conditions" << "\n";
		new_range = new double[ip.dims[0]*2];
		new_seg = new int[ip.dims[0]];
		rescale(best_run*ip.dims[0], ip.dims[0], ip.segments, rp.segs, new_seg, ranges, new_range);
		new_cubes->add(new_seg, new_range, best_grade );
	}	
	return success;	
} 

/*	Takes an LHS cube (segment coordinates), the range from which it was sampled, and an array for where to
store the new range. It then puts into the new range the bounds of the cube, which will be used for deeper
recursion.*/
void rescale(int start, int dims, int segments, int* segs, int* new_seg, double* ranges, double* new_range){
	double low;
	double high;
	for(int i = 0; i < dims; i++){
		new_seg[i] = segs[start+i];
		low = ranges[2*i];
		high = ranges[2*i+1];
		if(low == 0 && high == 0){
			new_range[2*i] = 0;
			new_range[2*i+1] = 0;
		}
		else{
			new_range[2*i] = low + ((high - low)/(double)segments)*((double)segs[start+i] - 1) ;
			new_range[2*i +1] = low + ((high-low) / (double)segments)*((double)segs[start+i]);
		}
	}
}

/*	Getting command line arguments from the user, which are stored in input_params& ip.*/
void accept_params (int num_args, char** args, input_params& ip) {
	int sim_args_index = 0;
	bool just_help = false;
	if (num_args > 1) { // if arguments were given and each argument option is followed by a value
		for (int i = 1; i < num_args; i += 2) { // iterate through each argument pair
			char* option = args[i];
			char* value;
			if (i < num_args - 1) {
				value = args[i + 1];
			} else {
				value = NULL;
			}
				
			/*
			Check for each possible argument option and overrides the default value for each specified option. If the option isn't recognized or the value given for an option doesn't appear valid then the usage information for the program is printed with an error message and no simulations are run. The code should be fairly self-explanatory with a few exceptions:
			1) atoi converts a string to an integer, atof converts a string to a floating point number (i.e. rational)
			2) strings should always be compared using strcmp, not ==, and strcmp returns 0 if the two strings match
			3) usage(true) prints the usage information with an error message while usage(false) prints it without one
			*/
			if(ip.sim_args){
				ip.simulation_args[sim_args_index] = option;
				sim_args_index++;
				i--;
				if (strcmp(option, "-h") == 0 || strcmp(option, "--help") == 0){
					just_help = true;
				}
			} else if (strcmp(option, "-e") == 0 || strcmp(option, "--exec") == 0) {
				ensure_nonempty(option, value);
				ip.sim_exec = value;
			} else if (strcmp(option, "-r") == 0 || strcmp(option, "--ranges-file") == 0) {
				ensure_nonempty(option, value);
				ip.ranges_file = value;
			} else if (strcmp(option, "-o") == 0 || strcmp(option, "--output-file") == 0) {
				ensure_nonempty(option, value);
				ip.params_file = value;
			} else if (strcmp(option, "-c") == 0 || strcmp(option, "--cubes-file") == 0) {
				ensure_nonempty(option, value);
				ip.cubes_file = value;
			} else if (strcmp(option, "-m") == 0 || strcmp(option, "--segments-file") == 0) {
				ensure_nonempty(option, value);
				ip.segments_file = value;
				ip.write_segments = 1;
			} else if (strcmp(option, "-v") == 0 || strcmp(option, "--verbose-file") == 0) {
				if(value == NULL || value[0] == '-'){
					i--;
				} else{
					ip.verbose_file = value;
				}
				ip.verbose = true;
			} else if (strcmp(option, "-n") == 0 || strcmp(option, "--segments-number") == 0) {
				ensure_nonempty(option, value);
				ip.segments = atoi(value);
				if (ip.segments < 1) {
					usage("You must use a postivie, non-zero, integer number of segments.",0);
				}
			} else if (strcmp(option, "-d") == 0 || strcmp(option, "--duplication-factor") == 0) {
				ensure_nonempty(option, value);
				ip.duplication = atoi(value);
				if (ip.duplication < 1) {
					usage("You must use a postivie, non-zero integer for the duplicaiton factor.", 0);
				}
			} else if (strcmp(option, "-D") == 0 || strcmp(option, "--recursion-depth") == 0) {
				ensure_nonempty(option, value);
				ip.max_depth = atoi(value);
				if (ip.max_depth < 1) {
					usage("You must use a postivie, non-zero integer for the recursion depth.", 0);
				} else if (ip.max_depth > 5){
					usage("I see you like to live dangerously. \n But seriously, that amount of recursion will break things.", 0);
				}
			} else if (strcmp(option, "-l") == 0 || strcmp(option, "--lhs-runs") == 0) {
				ensure_nonempty(option, value);
				ip.lhs_runs = atoi(value);
				if (ip.lhs_runs < 1) {
					usage("You must use a postivie, non-zero integer for the number of lhs runs you would like to perform.", 0);
				}
			} else if (strcmp(option, "-s") == 0 || strcmp(option, "--random-seed") == 0) {
				ensure_nonempty(option, value);
				ip.random_seed = atoi(value);
				if (ip.random_seed < 1) {
					usage("You must use a postivie, non-zero integer for the ranodm seed you would like to perform.", 0);
				}
			} else if (strcmp(option, "-R") == 0 || strcmp(option, "--raw-data") == 0) {
				ensure_nonempty(option, value);
				ip.raw_data = atoi(value);
				if (ip.raw_data < 1) {
					usage("You must use a postivie, non-zero integer for the raw data quantity.", 0);
				}
			} else if (strcmp(option, "-p") == 0 || strcmp(option, "--processes") == 0) {
				ensure_nonempty(option, value);
				ip.processes = atoi(value);
				if (ip.processes < 1) {
					usage("I doubt you want a negative amount of processes to run.", 0);
				}
			}else if (strcmp(option, "-t") == 0 || strcmp(option, "--condition-threshold") == 0) {
				ensure_nonempty(option, value);
				ip.condition_threshold = atof(value);
				if (ip.condition_threshold > 1 || ip.condition_threshold <0) {
					usage("You must use a positive fraction less than or equal to one for the mutant threshold", 0);
				}
			}else if (strcmp(option, "-T") == 0 || strcmp(option, "--cube-threshold") == 0) {
				ensure_nonempty(option, value);
				ip.cube_threshold = atof(value);
				if (ip.cube_threshold < 0 || ip.cube_threshold >1) {
					usage("You must use a positive fraction less than or equal to one for the cube threshold", 0);
				}
			} else if (strcmp(option, "-S") == 0 || strcmp(option, "--simple") == 0) {
				ip.lhs = false;
				i--;
			} else if (strcmp(option, "-b") == 0 || strcmp(option, "--take-best") == 0) {
				ip.take_best = true;
				i--;
			} else if (strcmp(option, "-i") == 0 || strcmp(option, "--increase-threshold") == 0) {
				ip.increase_threshold = true;
				i--;
			} else if (strcmp(option, "-w") == 0 || strcmp(option, "--write-segments") == 0) {
				ip.write_segments = true;
				i--;
			} else if (strcmp(option, "-q") == 0 || strcmp(option, "--quiet") == 0) {
				ip.quiet = true;
				i--;
			} else if (strcmp(option, "-a") == 0 || strcmp(option, "--sim-args") == 0) {
				ip.sim_args = true;
				ip.simulation_args = new char*[num_args - i  + 5];
				ip.sim_args_num = num_args - i  + 5;
				sim_args_index = 5;
				i--;
			} else if (strcmp(option, "-h") == 0 || strcmp(option, "--help") == 0) {
				const char* mess = "Welcome to the help options.\n Possible command line arguments are:\n"; 
				usage(mess,0);	
				i--;
			}
		}
	}
	//Ensures that there won't be processes made with nothing to do.
	if(ip.segments < ip.processes) ip.processes = ip.segments;	
	
	if(just_help){
		ip.segments = 1;
		ip.lhs_runs = 1;
		ip.max_depth = 1;
		ip.processes = 1;
		ip.quiet = true; 
	}
	//Setting up quiet mode.
	if(ip.quiet) cout_switch(true, ip); 
	
	//Use minimum system resources if the user just wants deterministic help menu.
	
	//Initializing some arguments that are always passed into the simulation program.
	if(!ip.sim_args){
		ip.simulation_args = new char*[6];
		ip.sim_args_num = 6;
		sim_args_index = 5;
	}	
	ip.simulation_args[0] = ip.sim_exec;
	ip.simulation_args[1] = (char*)"--pipe-in";
	ip.simulation_args[3] = (char*)"--pipe-out";
	ip.simulation_args[sim_args_index] = NULL;
	//Initializing the random seed.
	init_seed(ip);
}

void ensure_nonempty (const char* flag, const char* arg) {
	if (arg == NULL) {
		char* message = (char*)malloc(strlen("Missing argument for '' flag.") + strlen(flag) + 1);
		sprintf(message, "Missing the argument for the '%s' flag.", flag);
		usage(message, 0);
	}
}

//This function enables/disables quiet mode by redirecting cout to devnull.
void cout_switch(bool turn_off, input_params& ip){ 
	if(turn_off){
		ip.cout_orig = cout.rdbuf();
		ip.null_stream = new ofstream("/dev/null");
		cout.rdbuf(ip.null_stream->rdbuf());
	} else{
		cout.rdbuf(ip.cout_orig);			
	}
}

/*	This function is used to ensure that this program can use the maximum number of file descriptors allowed by the system.
	This may be necessary for large amounts of parallelization, but the max limit should only be a problem if files are not being properly closed.*/
void set_resources(){ 
	rlimit file_limit; 
	getrlimit(RLIMIT_NOFILE, &file_limit);
	file_limit.rlim_cur = file_limit.rlim_max;
	setrlimit(RLIMIT_NOFILE, &file_limit);
}

void usage (const char* message, int error) {
	if (strlen(message) != 0) { // if there is an error message to print then print it
		cout << message << "\n" << "Exit status: " << error << "\n";
	}	
	cout << "Usage: [-option [value]]. . . [--option [value]]. . ." << "\n";
	cout << "-a, --sim-args              [arguments]      : INCLUDING THIS WILL PASS ALL SUCCEEDING ARGUMENTS TO THE SIMULATION PROGRAM.\n"; 
	cout << "-e, --exec                  [path]           : The path of the executable that should run simulations of parameter sets. Default is \"../deterministic\". \n";
	cout << "-r, --ranges-file           [filename]       : The .ranges file from which the ranges to sample from will be taken. Default is \"sampleData/initial.ranges\".\n";
	cout << "-o, --output-file           [filename]       : The .params file to which the output of the sampling will be written. Default is \"sampleData/initial.params\". \n";
	cout << "-c, --cubes-file            [filename]       : The .cubes file to which the successful ranges will be written. Default is \"sampleData/initial.cubes\".\n";
	cout << "-m, --segments-file         [filename]       : The .segments file to which the indicies of the IHS-chosen segments (hypercubes) can be written. Default is \"sampleData/initial.segments\".\n";
	cout << "-v, --verbose-file          [filename]       : The .params file you would like to write all parameter set that were run. Default is disabled.\n"; 
	cout << "-d, --duplication-factor    [int]            : Duplication factor. A positive, non-zero integer. More duplication means more attempts at well-distributed LHS, but also means more time+memory used. \n";
	cout << "-n, --segments-number       [int]            : Number of segments (hypercubes) that should be chosen from the parameter space. This is equivalent to the number of points you would like to sample. \n" ;
	cout << "-l, --lhs-runs              [int]            : Number of times you would like run the hypercube sampling. While the number of segments determines how many cubes are in the space, the number of runs determines how many times to resample from the same space.\n";
	cout << "-s, --random-seed           [int]            : Seed for random number generator, a postivite non-zero integer. If unset, a seed will be generated based on system time. Set this value to get the same output from multiple runs of the program \n"; 
	cout << "-p, --processes             [int]            : For parallelizing, this specifies how many processes the simulation should be spread accross. Default is 1, which is just serial simulation.\n"; 
	cout << "-D, --recursion-depth       [int]            : How deep the program should recurse looking for parameter sets. The range of the resulting values will be (ranges)*n^D.\n"; 
	cout << "-R, --raw-data              [int]            : If you would like to just generate raw parameter sets and write them to the -o file, use this.\n"; 
	cout << "-t, --condition-threshold   [float]          : Percentage of mutants which must pass for a cube to be considered valid. \n";
	cout << "-T, --cube-threshold        [float]          : Percentage of cubes (ranges) within a space that must pass the mutant threshold to make that space valid. \n";	
	cout << "-S, --simple                                 : Including this will cause the sampling to act without LHS and instead use simply random values within the range. \n";   
	cout << "-b, --take-best                              : Including this will cause the simulation to analyze the best cube found when no cube in the space passed the mutant threshold. Disabled by default. \n"; 
	cout << "-i, --increase-threshold                     : Including this will increase the threshold for how many mutants have to pass (-M) at each depth of the recursive search. Disabled by default. \n"; 
	cout << "-w, --write-segments                         : Including this will enable writing the segment selections to a file. Automatically activated if a --segments-file is used. \n"; 
	cout << "-q, --quiet                                  : Turn off printing messages to standard output (for sampling program only). Disabled by default. \n"; 
	cout << "-h, --help                                   : This help menu. \n";
	cout << "\n" << "Example: ./ihsample -r data.ranges -o data.params --segments-number 10 -p 4 -t .66 -T .75 -i -random-seed 314 --sim-args -q"<< "\n" << "\n";
	exit(error);
}

/*
	The following is the non-LHS code.
*/

void simple_sample(input_params& ip,  recurse_params& rp , double* ranges , int* range_map){
	int param_num = ip.dims[0] + ip.dims[1];
	double* successes[param_num];
	int i = 0;
	for(; i < param_num; i++){
		successes[i] = new double[ip.num_needed];
		for(int j = 0; j < ip.num_needed; j++){
			successes[i][j] = 0;
		}
	}
	i = 0;
	for(; i < ip.lhs_runs && rp.good_samples < ip.num_needed; i++){
		best_sets(ip, rp, ranges, range_map, successes);
		simple_write(i, param_num, ip.segments*i, ip.num_needed, successes, ip.params_file);
		cut_ranges(ip, successes, ranges);		
	}
	simple_write(i, param_num, ip.segments*i, ip.num_needed, successes, ip.params_file);
	
}

void rand_sample(input_params& ip, recurse_params& rp, int* map, double* ranges){
	int param_num = ip.dims[0] + ip.dims[1];
	int i = 0;
	int j = 0;
	int mloc = 0;
	int zshift = 0;
	double min, max;
	//Cycles through every parameter set.
	for(; i < ip.segments; i++){
		//Cycles through every dimension/parameter in each set.
		for(; j < ip.dims[0]; j++){
			//Gets the appropriate range values.
			min = ranges[j*2];
			max = ranges[j*2 + 1];

			while(map[mloc] > j+zshift){
				rp.samples[i*(param_num) + j+zshift] = 0;
				zshift++;	
			}
			mloc++;
			//Chooses a random value from within the hypercube.
			rp.samples[i*(param_num) + j+zshift] = rand_range(min, max); 
		}
		while(param_num > j+zshift){
			rp.samples[i*(param_num) + j+zshift] = 0;
			zshift++;	
		}
		mloc = 0;
		zshift = 0;
		j=0;	
	}
}

void best_sets(input_params& ip, recurse_params& rp, double* ranges, int* range_map, double** successes){
	int param_num = ip.dims[0] + ip.dims[1];
	int i = 0;
	int j = 0;
	int len = 0;
	int sofar = 0;
	rp.good_samples = 0;
	//Creats random sets and simulates them until it gets enough that are close to the maximum per set.
	while(sofar < ip.num_needed){
		rand_sample(ip, rp, range_map, ranges);
		if(ip.verbose) sample_write(1, param_num, ip.segments, ip.duplication, ip.lhs_runs, rp.passing_grade, rp.samples, NULL, ip.verbose_file );
		simulate_samples(ip, rp );
		/*Print out the samples if there was a child error.
		if(rp.failure != NULL){
			for(int s = 0; s < ip.segments; s++);
		}*/
						
		len = count_success(ip, rp);
		i = 0;
		for(; i < len && i+sofar < ip.num_needed; i++){
			j = 0;
			if(rp.grades[i] >= ip.cube_threshold*(double)(rp.passing_grade)){
				for(; j < param_num; j++){
					successes[j][i+sofar] = rp.samples[i*param_num + j];
				}
			}	
		}
		sofar = i;
	}
}

void cut_ranges(input_params& ip, double** successes, double* ranges){
	double mean;
	double s_dev;
	double min_max[2];
	int i=0;
	for(; i < ip.segments; i++){
		mean = get_mean(successes[i], ip.num_needed);
		s_dev =  std_dev(successes[i], mean, ip.num_needed);
		get_min_max(successes[i], ip.num_needed, min_max);
		ranges[2*i] = (ranges[2*i] < min_max[0] - ip.s_devs*s_dev ? min_max[0] - ip.s_devs*s_dev : ranges[2*i]); 
		ranges[2*i] = (ranges[2*i+1] > min_max[1] + ip.s_devs*s_dev ? min_max[1] + ip.s_devs*s_dev : ranges[2*i+1]);
	}		
}

int count_success(input_params& ip, recurse_params& rp){
	int count = 0;
	int best_grade = 0;
	for(int i = 0; i < ip.segments; i++){
		best_grade = (rp.grades[i] >= best_grade ? rp.grades[i] : best_grade);
		if(rp.grades[i] == rp.passing_grade ) rp.good_samples++;
	}
	for(int i = 0; i < ip.segments; i++){
		if(rp.grades[i] >= ip.cube_threshold*(double)best_grade){
			count++;
		}	
	}
	rp.passing_grade = best_grade;
	return count;
}

void make_raw(input_params& ip, recurse_params& rp, int* range_map, double* ranges){
	int excess = (ip.raw_data % ip.segments);
	for(int i = 0; i < ip.raw_data/ip.segments; i++){
		rand_sample(ip, rp, range_map, ranges);
		sample_write(i, ip.dims[0] + ip.dims[1], ip.segments, ip.duplication, ip.lhs_runs, rp.passing_grade, rp.samples, NULL, ip.params_file );
	}
	ip.segments = excess;
	for(int i = 0; i < excess; i++){
		rand_sample(ip, rp, range_map, ranges);
		sample_write(i, ip.dims[0] + ip.dims[1], ip.segments, ip.duplication, ip.lhs_runs, rp.passing_grade, rp.samples, NULL, ip.params_file );
	}
}
