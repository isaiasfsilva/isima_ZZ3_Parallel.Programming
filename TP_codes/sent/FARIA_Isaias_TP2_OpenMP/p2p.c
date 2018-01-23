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
		This is just a communication test using MPI

INSTRUCTIONS:
		How to compile: mpicc –o p2p p2p.c
		How to execute: mpirun -np 32 -hostfile myhostfile p2p
			
			
Last update: 22 january 2018

This project is also available at Github
	>> https://goo.gl/23qq51 

*/


#include "mpi.h"
#include <stdio.h>
#include <string.h>
int main(int argc, char **argv){

//  // ========================  QUESTION 5 ========================= //
	// In this question we can see the effect of the parameter SOURCE //
	// We can define q source to receive some message OR we can listen//
	// from everyone :)                                               //
	// ===============================================================//
	int rang, nbProcs, dest=0, source, etiquette = 50;
	MPI_Status statut;
	char message[100];
	MPI_Init( &argc, &argv );
	MPI_Comm_rank( MPI_COMM_WORLD, &rang );
	MPI_Comm_size( MPI_COMM_WORLD, &nbProcs );
	if ( rang != 0 ) {
		sprintf( message, "Bonjour de la part de P%d!\n" , rang	);
		MPI_Send( message, strlen(message)+1, MPI_CHAR,dest, etiquette, MPI_COMM_WORLD );
	}else
		for ( source=1; source<nbProcs; source++ ) {
			MPI_Recv( message, 100, MPI_CHAR, MPI_ANY_SOURCE,
			etiquette, MPI_COMM_WORLD, &statut );
			printf( "%s", message );
		}
	MPI_Finalize();

	return 0 ;
}