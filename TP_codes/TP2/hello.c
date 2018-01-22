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
		This is just a "hello world!" using MPI

INSTRUCTIONS:
		How to compile: mpicc –o hello hello.c
		How to execute: ./hello
			
			
Last update: 22 january 2018

This project is also available at Github
	>> https://goo.gl/23qq51 
	
*/



#include "mpi.h"
#include <stdio.h>
int main(int argc, char **argv){

//  // ================== QUESTIONS 1 AND 2 ======================= //
	// These questions were made successfully.                      //
	// ============================================================ //


//  // ======================= QUESTION 3 ========================= //
	// Executing with: mpirun -np 4 hello                           //
	//  Hello from proc. 0 of 4 at processor node28.rcisima.isima.fr// 
	//  Hello from proc. 2 of 4 at processor node28.rcisima.isima.fr//
 	//  Hello from proc. 1 of 4 at processor node28.rcisima.isima.fr// 
 	//  Hello from proc. 3 of 4 at processor node28.rcisima.isima.fr// 
 	// -------------------------------------------------------------//
	// Executing with: mpirun -np 4 hello                           //
	//  Hello from proc. 0 of 4 at processor node28.rcisima.isima.fr// 
	//  Hello from proc. 2 of 4 at processor node28.rcisima.isima.fr//
 	//  Hello from proc. 1 of 4 at processor node28.rcisima.isima.fr// 
 	//  Hello from proc. 3 of 4 at processor node28.rcisima.isima.fr// 
 	// =============================================================//


//aux variables
	char name[200]; // To save the processor's name
	int len;		// The size of the processor's name

	int rang, nbProcs;
	MPI_Init( &argc, &argv );
	MPI_Comm_rank( MPI_COMM_WORLD, &rang );
	MPI_Comm_size( MPI_COMM_WORLD, &nbProcs );
	MPI_Get_processor_name(name, &len);
	printf( " Hello from proc. %d of %d at processor %s \n ",rang, nbProcs,name);
	

	MPI_Finalize();
	return 0;
}