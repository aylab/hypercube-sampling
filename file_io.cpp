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

#include <cerrno>
#include <cstring>
#include <iostream>
#include <sys/stat.h>
#include <fstream>
#include <cstdlib>
#include <string.h>
using namespace std;

#include "file_io.hpp"

char* read_file(char* fname){ // (AAy: ??? Reading the file to the buffer. Which file are we reading here? Ranges file probably.)
	FILE* file_pt;
	long long_size;
	size_t result;
	
	file_pt = fopen(fname, "rb");
	if (file_pt == NULL) {
		cout << "Couldn't open file " << fname << "!" << endl;
		exit(1);
	}
	
	// Obtain file size and check that it is under 400MB 
	fseek(file_pt, 0, SEEK_END);
	long_size = ftell(file_pt);	
	rewind(file_pt);
	
	// Allocate memory to contain the whole file 
	char* buffer = (char*)malloc(sizeof(char) * long_size + 1);
	if (buffer == NULL) {
		cout << "no memory"<< endl;
		exit(2);
	}
	
	// Copy the file into the buffer
	result = fread(buffer, 1, long_size, file_pt);
	if ((int)result != long_size) {
		cout << "Reading error" << endl;
		exit(3);
	}
	buffer[long_size] = '\0';
	// The whole file is now loaded in the memory buffer
	fclose(file_pt);
	return buffer;
}

int count_dims(char* raw, int* dim_ig){ // (AAy: ??? How many ranges should be zeroed out. How does this function work?)
	int i = 0, dims = 0, ignore = 0, open = 0;
	while(raw[i] != '\0'){
		while(raw[i] == '#'){
			while(raw[i] != '\n' && raw[i] != '\0') i++;			
		}
		if(raw[i] == '\0') break;
	
		while(raw[i] != '[' && raw[i] != '\0') i++;
		if(raw[i] == '\0') break;
		if(raw[i] == '['){
			while(raw[i] != '-' && raw[i] != ']' && raw[i] !='\0') i++;
			if(raw[i] == ']'){
				dims++;				
			} else if(raw[i] == '-'){
				ignore++;
			} else{
				open++;
				break;
			}
		} 
		i++;
	}
	dim_ig[0] = dims;
	dim_ig[1] = ignore;
	//This just returns negative dims if there was a bracket mismatch.
	if(open!=0) return -1*(1+dims+ignore);
	return dims;
}

void parse(int dim_num, char* raw, double* parsed, int* map){ // (AAy: ??? Parsing the parameter set. What is this function doing, how does it work?)
	if(raw == NULL) return;
	/*
	Assumes the file format:
	#this line will be ignored

	parameter1 [min,max]
	parameter2 [min,max] 
	parameterToZeroOut [ -min, max]

	Only the brackets and negation for zeroing are essential. 	
	*/
	int i = 0;
	int ploc = 0;
	int mloc = 0;
	int zshift = 0;
	double temp = 0;
	for(; raw[i]!='\0' ; i++){
		while(raw[i] == '#'){
			while(raw[i] != '\n' && raw[i] != '\0') i++; // (AAy: ??? Why do we have i++ twice?)
			i++;			
		}
		
		while(raw[i] != '[' && raw[i] != '\0') i++;
		if(raw[i] == '\0') break;
		i++;
		
		temp = atof(raw+i); // (AAy: ??? Returns the value at raw[i])		
		if(temp >= 0){ 
			parsed[ploc] = temp;
			map[mloc] = mloc+zshift; // (AAy: ??? What is this part of the code doing?)
			mloc++;
			ploc++;	
		
			while(raw[i] != ',') i++; 
			i++;
			parsed[ploc] = atof(raw+i);
			ploc++;	
		} else{
			zshift ++;
		}	
		while(raw[i] != '\n' && raw[i] != '\0') i++;		 
	}
}


//Write the scaled samples to a file.
void sample_write ( int append, int dim_num, int segments, int duplication, int lhs_runs, int passing_grade, double* samples, int* grades, char *file_out_name ){
	ofstream file_out;
	int i;
	int j;
	double x;
	
	if(append > 0){ // (AAy: ??? Open the file in the append mode)
		file_out.open ( file_out_name, ios::out | ios::app );
	} else{
		file_out.open(file_out_name);
	}
	//cout << "Sample file: " << file_out_name;
	if ( !file_out ) 
	{
		cout << "\n";
		cout << "  Could not open the output file.\n";
		exit ( 1 );
	}

	if(append == 0 && grades!=NULL){
		file_out << "#  " << file_out_name << "\n";
		file_out << "#  Spatial dimension DIM_NUM = "  << dim_num      << "\n";
		file_out << "#  Number of points N =        "  << segments     << "\n";
		file_out << "#  Duplication factor D =      "  << duplication  << "\n";
		file_out << "#\n";
	}
	for ( j = 0; j < segments; j++ )
	{
		if(grades != NULL) file_out << "#Grade: " << grades[j] << "/" << passing_grade << " = " << grades[j]/(double)passing_grade << "% of conditions.\n";
		for ( i = 0; i < dim_num; i++ )
		{
			x = ( double ) samples[i+j*dim_num];
			file_out << x;
			if(i <  dim_num -1) file_out << ",";
		}
		file_out << "\n";
	}
	
	file_out.close( ); // (AAy: ??? Writing out the parameter sets. How does it work?)

return;
}

//Write the scaled samples to a file.
void single_sample_write (int append, int dim_num, int passing_grade, int grade, double* samples, char *file_out_name ){
	ofstream file_out;
	int i;
	double x;
	
	if(append > 0){
		file_out.open ( file_out_name, ios::out | ios::app );
	} else{
		file_out.open(file_out_name);
	}
	//cout << "Sample file: " << file_out_name;
	if ( !file_out )
	{
		cout << "\n";
		cout << "  Could not open the output file.\n";
		exit ( 1 );
	}

	if(append == 0){
		file_out << "#  " << file_out_name << "\n";
		file_out << "#  Spatial dimension DIM_NUM = "  << dim_num      << "\n";
		file_out << "#\n";
	}

	file_out << "#Grade: " << grade << "/" << passing_grade << " = " << 100*(double)grade/(double)passing_grade << "% of conditions.\n";
	for ( i = 0; i < dim_num; i++ )
	{
		x = ( double ) samples[i];
		file_out<< x;
		if(i <dim_num -1) file_out << ",";
	}
	file_out << "\n";

	file_out.close( );
	return;
}


//Write the valid ranges to a file.
void range_write ( int append, int dim_num, int segments, int duplication, double grade, double* ranges, char *file_out_name )
{
	ofstream file_out;
	int i;
	double x;
	
	if(append > 0){
		file_out.open ( file_out_name, ios::out | ios::app );
	} else{
		file_out.open(file_out_name);
	}
	if ( !file_out )
	{
		cout << "Range write error: \n";
		cout << "  Could not open the output file.\n";
		exit ( 1 );
	}

	if(append == 0){
		file_out << "#  " << file_out_name << "\n";
		file_out << "#  Spatial dimension DIM_NUM = "  << dim_num      << "\n";
		file_out << "#  Number of points N =        "  << segments     << "\n";
		file_out << "#  Duplication factor D =      "  << duplication  << "\n";
		file_out << "#\n";

	}
	file_out << "#\tGrade: " << grade << endl;
	for ( i = 0; i < dim_num; i++ )
	{
		x = ( double ) ranges[i];
		file_out << "[" << x << ",";
		i++;
		x = (double) ranges[i];
		file_out << " " << x << "]\n";	
	}
	file_out << "\n";

	file_out.close( );

return;
}

//For the simple random sampling method, write the scaled samples to a file.
void simple_write ( int append, int dim_num, int num_runs, int num_samples, double** samples, char *file_out_name ){
	ofstream file_out;
	int i;
	int j;
	double x;
	
	if(append > 0){
		file_out.open ( file_out_name, ios::out | ios::app );
	} else{
		file_out.open(file_out_name);
	}
	if ( !file_out )
	{
		cout << "\n";
		cout << "  Could not open the output file.\n";
		exit ( 1 );
	}

	if(append == 0){
		file_out << "#  " << file_out_name << "\n";
		file_out << "#  Number of runs performed =        "  << num_runs     << "\n";
		file_out << "#\n";
	}
	for ( j = 0; j < num_samples; j++ )
	{
		for ( i = 0; i < dim_num; i++ )
		{
			x = ( double ) samples[i][j];
			file_out<< x;
			if(i <dim_num -1) file_out << ",";
		}
		file_out << "\n";
	}
	
	file_out.close( );

return;
}

