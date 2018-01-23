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
		How to compile: mpicc –o p2p_Question6_bin p2p_Question6_bin.c -lm
		How to execute: mpirun -np 32 -hostfile myhostfile p2p_Question6_bin
			
			
Last update: 23 january 2018

This project is also available at Github
	>> https://goo.gl/23qq51 

*/

#include "mpi.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
int main(int argc, char **argv){


	int rang, nbProcs, dest=0, source, etiquette = 50;
	
	//Id of next and previous levels
	int prevID,nextID, squareNB,step;
	int myDegree=0;
	//aux variables
	int j,i;
	MPI_Status statut;
	char message[300],messageR[300];
	char messageS[300];
	sprintf( message, " ");
	MPI_Init( &argc, &argv );
	MPI_Comm_rank( MPI_COMM_WORLD, &rang );
	MPI_Comm_size( MPI_COMM_WORLD, &nbProcs );

	squareNB = sqrt(nbProcs); //NB of layers
    if(squareNB*squareNB!=nbProcs)
    	squareNB++;
	sprintf( message, "-");

// ===============================//
//      RECEIVING SECTION         //
//================================//
	for(i=0; i<=squareNB;i++){ // Find the notes that has to receive something at each layer
		step = pow(2,i+1);
		if(rang%step==0){
			nextID = rang^(1<<i);
			if(nextID > rang && nextID < nbProcs){
				//I have to receive somethinf from nextID
				MPI_Recv( messageR, 300, MPI_CHAR, nextID,etiquette, MPI_COMM_WORLD, &statut );

				//I use a little text identifier "-" to concatenate the new string 
				// with my buffer
				int idx = 0;
				while (message[idx] != '-') idx++;
				sprintf( message, "%*.*s%s" ,0,idx,message,messageR);

				//Number of received messages. This important to calc 
				//the next step to send the message
				myDegree++;

			}

		}


	}


// ===============================//
//        SENDING SECTION         //
//================================//
	//Now they have to sent the message
	if(rang==0){

		//Node 0 just print the buffer

		//before we haave to clean the fist common "," and the text identifier "-"
		message[0]=' ';
		int idx = 0;
		while (message[idx] != '-') idx++;
		message[idx]=' ';

		//print bufffer
		printf( "Bonjour de la part de%s", message );
	}else{
		//calc the node to send the buffer + my contribuition
		//identify my level
		while(squareNB>0 && (rang%(int)pow(2,squareNB))!=0){			
			squareNB--;
		}
		//calc the next neightbors		
	    prevID = rang^(1<<squareNB);

		//Put my contribuition
//----------------------------------------------------//
// This part is so much important!!! Here we desloc 5 //
// caracters to right and put in the first 5 positions//
// my contribuition P MYID.                           //
//----------------------------------------------------//
	
		sprintf( messageS, ", P%d%*.*s" , rang,1,299,message);
		
		//Send the message
		MPI_Send( messageS,strlen(messageS)+1, MPI_CHAR,prevID, etiquette, MPI_COMM_WORLD );
	}

	
	MPI_Finalize();
	return 0 ;
}
