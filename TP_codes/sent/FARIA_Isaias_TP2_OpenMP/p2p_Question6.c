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
		How to compile: mpicc –o p2p_Question6 p2p_Question6.c
		How to execute: mpirun -np 32 -hostfile myhostfile p2p_Question6
			
			
Last update: 22 january 2018

This project is also available at Github
	>> https://goo.gl/23qq51 

*/

#include "mpi.h"
#include <stdio.h>
#include <string.h>
#include <string.h>
int main(int argc, char **argv){


	int rang, nbProcs, dest=0, source, etiquette = 50;
	MPI_Status statut;
	char message[300];
	char messageS[300];
	sprintf( message, " ");
	MPI_Init( &argc, &argv );
	MPI_Comm_rank( MPI_COMM_WORLD, &rang );
	MPI_Comm_size( MPI_COMM_WORLD, &nbProcs );

	//Everyone except the last thread has to read something from its succ
	if(rang<(nbProcs-1)){
		MPI_Recv( message, 300, MPI_CHAR, rang+1,etiquette, MPI_COMM_WORLD, &statut );
	}	
	//Everyone except the first thread has to send something to its pred
	if ( rang != 0 ) {		
		
		//----------------------------------------------------//
		// This part is so much important!!! Here we desloc 5 //
		// caracters to right and put in the first 5 positions//
		// my contribuition P MYID.                           //
		//----------------------------------------------------//
		sprintf( messageS, ", P%d%*.*s" , rang,5,295,message);

		//send to my predecessor
		MPI_Send( messageS,strlen(messageS)+1, MPI_CHAR,rang-1, etiquette, MPI_COMM_WORLD );
	}else{
		//Just the first thread has to print something
		message[0]=' ';
		printf( "Bonjour de la part de%s", message );
	}

	MPI_Finalize();
	return 0 ;
}
