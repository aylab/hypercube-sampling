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

#ifndef SAMPLE_HPP
#define SAMPLE_HPP

#include <iostream>
#include <streambuf>

using namespace std;
/*
	Input parameters struct that is passed throughout most of the program.
	All of this is put on the stack except simulation_args, which is allocated in accept_params
	The memory that simulation_args contains pointers to is only on the stack, except items [2] and [4] 
	which are malloc'd and free'd appropriately in simulate_samples.
	Simulation_args is deleted in the destructor. 
*/
struct input_params{	
	bool lhs;
	bool write_segments;
	bool sim_args;
	bool increase_threshold;
	bool take_best;
	bool quiet;
	bool verbose;
	int segments; //A value of 1 makes LHS equivalent to getting a single random sample.
	int duplication; //More duplication means more attempts at well-distributed LHS, but also means more time+memory used.
	int lhs_runs;
	int random_seed;
	int max_depth;	
	int dims[2]; // (AAy: ??? dimensions of the parameter sets to take into account or ignore)
	int processes;
	int sim_args_num;
	int s_devs;
	int num_needed;
	int raw_data; 
	double condition_threshold;
	double cube_threshold;
	streambuf* cout_orig;
	ofstream* null_stream;
	char* ranges_file;
	char* params_file;
	char* segments_file;
	char* cubes_file;
	char* verbose_file;
	char* sim_exec;
	char** simulation_args;
	
	input_params(){
		lhs = true;		
		quiet = false;
		sim_args = false;
		write_segments = false;
		increase_threshold = false;
		take_best = false;
		verbose = false;
		segments = 10;
		duplication = 5;	
		lhs_runs = 1;
		random_seed = 0;
		max_depth = 3;	
		processes = 2;
		s_devs = 2;
		raw_data = 0;
		num_needed = (int)(.03*(double)segments) + 1;
		condition_threshold = .70;
		cube_threshold = .70;
		ranges_file = (char*)"sampleData/88.ranges";
		params_file = (char*)"sampleData/passed_88.params";
		segments_file = (char*)"sampleData/passed_88.segments";
		cubes_file = (char*)"sampleData/passed_88.cubes";	
		verbose_file = (char*)"sampleData/all_sets_used_88.params";			
		sim_exec = (char*) "../sogen-deterministic/deterministic";
		simulation_args = NULL;
	}
	
	~input_params(){
		if(simulation_args != NULL) delete[] simulation_args;
	}

/*	This function is used for increasing how many mutants need to pass depending on the depth of the recursion. 
	This is disabled by default, in which case the function just returns ip.conditions_threshold.
*/	
	double calc_threshold(int depth){
		if(increase_threshold){
			 return condition_threshold + (1 - condition_threshold )*((double)depth / ((double)depth + (double)max_depth));
		}
		return  condition_threshold;
	}
};

struct recurse_params{
	int depth; //Keeps track of which level of recursion the struct is currently in.
	int good_samples;
	int passing_grade; //The score that means a simulation passed all mutant tests.
	int* grades;  //The scores of each simulation.
	int* segs;
	double* samples;
	
	char* failure;
	int failcode;
	
	recurse_params(input_params& ip){
		depth = 0; 
		good_samples = 0;
		passing_grade = 0;
		grades = new int[ ip.segments ];  //This will store the scores of each simulation.
		for (int i = 0; i < ip.segments; i++) {
			grades[i] = 0;
		}
		segs = new int[ip.dims[0]*ip.segments]; //This will store the location of each box
		for (int i = 0; i < ip.dims[0]*ip.segments; i++) {
			segs[i] = 0;
		}
		samples = new double[(ip.dims[0]+ip.dims[1])*ip.segments]; //Initially zero out the array.
		for (int i = 0; i < (ip.dims[0]+ip.dims[1])*ip.segments; i++) {
			samples[i] = 0;
		}
		failure = NULL;
		failcode = 0;
	}
	~recurse_params(){
		delete[] grades;
		delete[] segs;
		delete[] samples;
		if(failure != NULL){
			free(failure);
		}
	}
	
};

struct range_node{
	double* ranges; //  Ranges of the parameter
	int* segment; //  Location of the box 
	double grade; // Score for this box
	int count; // Counter for the successful hits in a box
	
	range_node* next;
	range_node* prev;
	
	range_node(double* arr, double g){
		ranges = arr;
		next = NULL;
		prev = NULL;
		count = 1;
		grade = g;
	}
	range_node(int* seg, double* arr, double g){
		segment = seg;
		ranges = arr;
		next = NULL;
		prev = NULL;
		count = 1;
		grade = g;
	}
	~range_node(){
		delete[] ranges;
		delete[] segment;
	}
};

struct range_list{
	range_node* head;
	range_node* tail;
	range_list(){
		head = NULL;
		tail = head;
	}
	~range_list(){
		tail = head;
		while(head!=NULL){
			tail = head->next;
			delete head;
			head = tail;
		}
	}
	void add(double* arr, double g){
		if(head == NULL){
			head = new range_node(arr, g);
			tail = head;
		} else{
			tail->next = new range_node(arr, g);
			tail->next->prev = tail;
			tail = tail->next;
		}		
	}
	void add(int* seg, double* arr, double g){
		if(head == NULL){
			head = new range_node(seg, arr, g);
			tail = head;
		} else{
			tail->next = new range_node(seg, arr, g);
			tail->next->prev = tail;
			tail = tail->next;
		}		
	}
	//A simple linked-list concatonation function
	void concat(range_list* other){
		if(head == NULL){
			head = other->head;
			tail = other->tail;
		}else{
			tail->next = other->head;
			tail->next->prev = tail;
			tail = other->tail;
		}
	}
	//A function used to cut a successful box out of a temporary list and add it to the successfull list.
	void take(range_list* source, range_node* victim){ 
		if(victim == source->head){
			source->head = victim->next;
			source->head->prev = NULL;
		} else{
			victim->prev->next = victim->next;	
		}
		if(victim->next != NULL) victim->next->prev = victim->prev;
		if(this->head == NULL){
			this->head = victim;
			victim->prev = NULL;	
		} else{		
			this->tail->next = victim;
			victim->prev = this->tail;
		}
		victim->next = NULL;
		this->tail = victim;
	}
	//Used to check for redundant segments in the successful linked list. If a segment is successful multiple times it gets its count increased.
	bool find(int seg_size, int* seg){ 
		bool result = false;
		bool coord;
		range_node* subject = head;
		for(; subject != NULL && !result; subject = subject->next){
			coord = true;
			for(int i = 0; i < seg_size; i++){
				coord = coord && (seg[i] == subject->segment[i]);
			}
			result = result || coord; 
		} 
		if(result){
			subject->count +=1;
			this->rank_up(subject);
		}
		return result;
	}
	//Implemented this function so that the nodes will be ordered from most-found to least-found.
	void rank_up(range_node* bubble){
		if(bubble == this->head) return; //Check to see if the node is already in first place.
		
		range_node* superior = bubble->prev; //Search for the node to insert bubble next to.
		for(; superior->count < bubble->count && superior != NULL; superior = superior->prev);
		if(superior == bubble->prev ) return; //If this is true, bubble need not move.
		
		bubble->prev->next = bubble->next; //Cut bubble out of its current position
		if(bubble->next == NULL){
			 this->tail = bubble->prev;
		} else{
			bubble->next->prev = bubble->prev; //Only change the next node if it's not NULL.
		}
		//If you get to NULL, then bubble needs to become the new head.
		if(superior == NULL){
			bubble->next = this->head;
			this->head->prev = bubble;
			this->head = bubble;
		} else{
			bubble->next = superior->next;
			bubble->next->prev = bubble;
			bubble->prev  = superior;
			superior->next = bubble;
		}
	} 
};	

void init_seed (input_params& ip) {
	if (ip.random_seed == 0) {
		ip.random_seed = abs((((int)time(0) * 181) * ((getpid() - 83) * 359)) % 805306457); // the final seed value, calculated using the given seed and the process ID
	}
	srand(ip.random_seed);
}
double randouble(){
	return (double)rand() / (double)RAND_MAX;
}

void scale_sample(input_params&, int* , int*, double*, double* );
double rand_range(double, double);

void ih_sample(input_params& , int*);

void simulate_samples(input_params&, recurse_params&);
bool make_pipes(int, int**);
void del_pipes(int, int**, bool);
char*** make_args(input_params&, int**);
void del_args(input_params&, char***);
void segs_per_sim(int , int, int* );

void set_resources();
void accept_params (int, char** , input_params& );
void ensure_nonempty (const char* , const char* );
void cout_switch(bool , input_params& );
void usage (const char* , int ) ;

bool sample_recurse(input_params&,  recurse_params& , double* , int* , range_list*);
int analyze(input_params& , recurse_params& , double* ,range_list* );
void rescale(int , int , int , int* , int*, double* , double* );

void simple_sample(input_params&,  recurse_params& , double* , int* );
void best_sets(input_params& ip, recurse_params& rp, double* ranges, int* range_map, double** successes);
void rand_sample(input_params& , recurse_params& , int*, double* );
void cut_ranges(input_params& ip, double** successes, double* ranges);
int count_success(input_params& , recurse_params& );

void make_raw(input_params& ip, recurse_params& rp, int* , double* );	
#endif
