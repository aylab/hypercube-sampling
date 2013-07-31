Legal:

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


Project: 
	Zebrafish Somitogenesis Deterministic Model, 2013
Git repository (private as of June, 2013):
	https://github.com/aylab/sogen-deterministic
Authors: 
	Ahmet Ay, Jack Holland, Adriana Sperlea, Sebastian Sangervasi*
	*sampling code.
Description:
	This directory contains the code used for sampling parameter sets for the deterministic model. The code uses Latin Hypercube Sampling. 
	The sampling program borrows an implementation of Improved Hypercube Sampling from John Burkardt:
		ihs.cpp
		ihs.hpp
	The above code was made available by GNU LGPL licencse on J Burkardt's website. This code has not been changed from its original state. (7/1/2013) 
	
	The sampling code that uses this IHS implementation is specifically designed for use with the simulations available in the parent directory:
		sogen-deterministic
Usage:	
	Command-line arguments: [-option [value]]. . . [--option [value]]. . .
		-a, --sim-args           [arguments]      : INCLUDING THIS WILL PASS ALL SUCCEEDING ARGUMENTS TO THE SIMULATION PROGRAM.   
		-e, --exec               [path]           : The path of the executable that should run simulations of parameter sets. Default is   "../deterministic  ".   
	   	-r, --ranges-file        [filename]       : The .ranges file from which the ranges to sample from will be taken. Default is "sampleData/initial.ranges".  
	  	-o, --output-file        [filename]       : The .params file to which the output of the sampling will be written. Default is   "sampleData/initial.params". \
	  	-c, --cubes-file         [filename]       : The .cubes file to which the successful ranges will be written. Default is   "sampleData/initial.cubes  ".  
	   	-m, --segments-file      [filename]       : The .segments file to which the indicies of the IHS-chosen segments (hypercubes) can be written. Default is   "sampleData/initial.segments  ".  
	   	-d, --duplication-factor [int]            : Duplication factor. A positive, non-zero integer. More duplication means more attempts at well-distributed LHS, but also means more time+memory used.   
	   	-n, --segments-number    [int]            : Number of segments (hypercubes) that should be chosen from the parameter space. This is equivalent to the number of points you would like to sample. \n" ;
	   	-l, --lhs-runs           [int]            : Number of times you would like run the hypercube sampling. While the number of segments determines how many cubes are in the space, the number of runs determines how many times to resample from the same space.  
	   	-s, --random-seed        [int]            : Seed for random number generator, a postivite non-zero integer. If unset, a seed will be generated based on system time. Set this value to get the same output from multiple runs of the program    
	   	-p, --processes          [int]            : For parallelizing, this specifies how many processes the simulation should be spread accross. Default is 1, which is just serial simulation.   
	   	-D, --recursion-depth    [int]            : How deep the program should recurse looking for parameter sets. The range of the resulting values will be (ranges)*n^D.   
	   	-M, --mutant-threshold   [float]          : Percentage of mutants which must pass for a cube to be considered valid.   
	   	-C, --cube-threshold     [float]          : Percentage of cubes (ranges) within a space that must pass the mutant threshold to make that space valid.   	
	   	-b, --take-best                           : Including this will cause the simulation to analyze the best cube found when no cube in the space passed the mutant threshold. Disabled by default.    
	   	-i, --increase-threshold                  : Including this will increase the threshold for how many mutants have to pass (-M) at each depth of the recursive search. Disabled by default.    
	   	-w, --write-segments                      : Including this will enable writing the segment selections to a file. Automatically activated if a --segments-file is used.    
	   	-q, --quiet                               : Turn off printing messages to standard output (for sampling program only). Disabled by default.    
	   	-h, --help                                : This help menu.   
	Example: ./ihsample -r data.ranges -o data.params --segments-number 10 -p 4 -M .70 -i -random-seed 314 -C .66 --sim-args -q	

