/*
 Deterministic simulator for the zebrafish segmentation clock.
 Copyright (C) 2012 Ahmet Ay, Jack Holland, Adriana Sperlea
 
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

#include <math.h>

/*Generate the mean of a set, given its length*/
double get_mean(double *list, int listlen){
	double avgsum = 0;
	
	int i = 0;
	for(;i<listlen;i++){
		avgsum += list[i]/((double) listlen);
	}
	return avgsum;
}

void get_min_max(double* list, int listlen, double* min_max){
	int i = 1;
	double min = list[0];
	//printf("MINIMUM: %f\n",min);
	double max = list[0];
	for(;i<listlen;i++){
		//printf("List[%d]: %f\n",i,list[i]);
		if(list[i] < min){
			min = list[i];
		}
		if(list[i] > max){
			max = list[i];
		}
	}
	
	min_max[0] = min;
	min_max[1] = max;
}

/*Generate standard deviation of a given set when given the mean of that set and its length*/
double std_dev(double *list, double mean, int listlen){
	double sum = 0;
	//Assuming a normal distribution for these ranges
	
	int i = 0;
	for(;i<listlen;i++){
		int tmp = sum;
		sum += pow((list[i] - mean),2)/(listlen - 1); //Sample std dev including Bessel's correction
		
		if(sum<tmp){ //overflow
			return -1; //real std deviation is always an absolute value, this indicates an error
		}
	}
	
	return sqrt(sum);
}

