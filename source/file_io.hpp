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

#ifndef FILE_IO_HPP
#define FILE_IO_HPP

char* read_file(char* );
int count_dims(char*, int*);
void parse(int , char* , double*, int* );
void sample_write ( int , int , int , int, int, int, double* , int*, char*);
void single_sample_write (int , int , int , int , double* , char* );
void range_write ( int append, int dim_num, int segments, int duplication, double grade, double* ranges, char *file_out_name );
void simple_write ( int append, int dim_num, int num_runs, int num_samples, double** samples, char *file_out_name );
#endif
