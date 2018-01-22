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
		How to execute: mpirun -np 32 -hostfile myhostfile hello
			
			
Last update: 22 january 2018

This project is also available at Github
	>> https://goo.gl/23qq51 

*/



#include "mpi.h"
#include <stdio.h>
//Aditional libs section
   // Lib to use sched_getcpu()
   #include <sched.h>


int main(int argc, char **argv){

//  // ================== QUESTIONS 1 AND 2 ======================= //
	// These questions were made successfully.                      //
	// ============================================================ //


//  // ======================= QUESTION 3 ========================= //
	// Executing with: mpirun -np 4 hello                           //
    //  Hello from proc. 0 of 4 at processor node28.rcisima.isima.fr//
    //    which id is 1                                             //
    //  Hello from proc. 2 of 4 at processor node28.rcisima.isima.fr//
    //    which id is 16                                            //
    //  Hello from proc. 3 of 4 at processor node28.rcisima.isima.fr//
    //    which id is 24                                            //
    //  Hello from proc. 1 of 4 at processor node28.rcisima.isima.fr//
    //    which id is 8                                             // 
 	// =============================================================//

//  // ====================== QUESTION 4 =========================== //
	// After had created the hostfile file, I executed with -np 32   //
	// and defining the hostfile and the output was:                 //
 	// Hello from proc. 25 of 32 at processor node28.rcisima.isima.fr//
	// which id is 11                                                //
 	// Hello from proc. 20 of 32 at processor node28.rcisima.isima.fr//
	// which id is 2                                                 //
 	// Hello from proc. 16 of 32 at processor node28.rcisima.isima.fr//
	// which id is 3                                                 //
 	// Hello from proc. 17 of 32 at processor node28.rcisima.isima.fr//
	// which id is 15                                                //
 	// Hello from proc. 21 of 32 at processor node28.rcisima.isima.fr//
	// which id is 14                                                //
 	// Hello from proc. 31 of 32 at processor node28.rcisima.isima.fr//
	// which id is 10                                                //
 	// Hello from proc. 18 of 32 at processor node28.rcisima.isima.fr//
	// which id is 1                                                 //
 	// Hello from proc. 19 of 32 at processor node28.rcisima.isima.fr//
	// which id is 24                                                //
 	// Hello from proc. 22 of 32 at processor node28.rcisima.isima.fr//
	// which id is 17                                                //
 	// Hello from proc. 23 of 32 at processor node28.rcisima.isima.fr//
	// which id is 9                                                 //
 	// Hello from proc. 27 of 32 at processor node28.rcisima.isima.fr//
	// which id is 12                                                //
 	// Hello from proc. 30 of 32 at processor node28.rcisima.isima.fr//
	// which id is 6                                                 //
 	// Hello from proc. 24 of 32 at processor node28.rcisima.isima.fr//
	// which id is 4                                                 //
 	// Hello from proc. 26 of 32 at processor node28.rcisima.isima.fr//
	// which id is 16                                                //
 	// Hello from proc. 28 of 32 at processor node28.rcisima.isima.fr//
	// which id is 0                                                 //
 	// Hello from proc. 29 of 32 at processor node28.rcisima.isima.fr//
	// which id is 8                                                 //
 	// Hello from proc. 0 of 32 at processor node26.rcisima.isima.fr //
	// which id is 7                                                 //
 	// Hello from proc. 2 of 32 at processor node26.rcisima.isima.fr //
	// which id is 1                                                 //
 	// Hello from proc. 4 of 32 at processor node26.rcisima.isima.fr //
	// which id is 20                                                //
 	// Hello from proc. 3 of 32 at processor node26.rcisima.isima.fr //
	// which id is 28                                                //
 	// Hello from proc. 5 of 32 at processor node26.rcisima.isima.fr //
	// which id is 10                                                //
 	// Hello from proc. 7 of 32 at processor node26.rcisima.isima.fr //
	// which id is 14                                                //
 	// Hello from proc. 1 of 32 at processor node26.rcisima.isima.fr //
	// which id is 24                                                //
 	// Hello from proc. 6 of 32 at processor node26.rcisima.isima.fr //
	// which id is 22                                                //
 	// Hello from proc. 8 of 32 at processor node26.rcisima.isima.fr //
	// which id is 5                                                 //
 	// Hello from proc. 9 of 32 at processor node26.rcisima.isima.fr //
	// which id is 15                                                //
 	// Hello from proc. 11 of 32 at processor node26.rcisima.isima.fr//
	// which id is 9                                                 //
 	// Hello from proc. 12 of 32 at processor node26.rcisima.isima.fr//
	// which id is 2                                                 //
 	// Hello from proc. 13 of 32 at processor node26.rcisima.isima.fr//
	// which id is 13                                                //
 	// Hello from proc. 14 of 32 at processor node26.rcisima.isima.fr//
	// which id is 0                                                 //
 	// Hello from proc. 10 of 32 at processor node26.rcisima.isima.fr//
	// which id is 6                                                 //
 	// Hello from proc. 15 of 32 at processor node26.rcisima.isima.fr//
	// which id is 8                                                 //
	// ==============================================================//

//aux variables
	char name[200]; // To save the processor's name
	int len;		// The size of the processor's name

	int rang, nbProcs;
	MPI_Init( &argc, &argv );
	MPI_Comm_rank( MPI_COMM_WORLD, &rang );
	MPI_Comm_size( MPI_COMM_WORLD, &nbProcs );
	MPI_Get_processor_name(name, &len);
	printf( " Hello from proc. %d of %d at processor %s which id is %d\n ",rang, nbProcs,name,sched_getcpu());
	

	MPI_Finalize();
	return 0;
}