/*

[2018] Isaias Faria Silva 
All Rights Reserved.
 
NOTICE: All information contained here in this file is the 
        property of Isaias Faria Silva. 
        If you want to use this, please put the contact 
        rights information.
 
PURPOSE:Practical activity of Parallel Programming
		Clermont-Auvergne University - Computer Science Institute
		ZZ3 - 2018 January 

DESCRIPTION:
		This is just a "hello world!" using OpenMPI

INSTRUCTIONS:
		How to compile: gcc -fopenmp hello.c -o hello
		How to execute: ./hello NBOFTHREADS
			ps: NBOFTHREADS is optional. The defaut value is 8 threads.
			
Last update: 15 january 2018


This project is also available at Github
	>> https://goo.gl/23qq51 
	
*/


#include <stdio.h>
#include "omp.h"

int main(int argc, char **argv)
{
	int NBofTHREADS = 8;
	if(argc > 1){
		NBofTHREADS = atoi(argv[1]);   
	}

	printf("Before PARALLEL REGION : There are %d threads\n\n" , omp_get_num_threads()) ;
	//defininf the num_of_threads
	#pragma omp parallel num_threads(NBofTHREADS)	
	{
		printf("In the PARALLEL REGION : There are %d threads\n\n" , omp_get_num_threads()) ;
		printf("Hello World from TID %d!\n", omp_get_thread_num());
	}
	printf("After PARALLEL REGION : There are %d threads\n\n" ,	omp_get_num_threads()) ;
}