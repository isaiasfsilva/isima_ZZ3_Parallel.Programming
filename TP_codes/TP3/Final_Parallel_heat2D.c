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
        This is a Parallell version of the TP3 - Using MPI_Scatter and MPI_Gatter

             Communication locale / globale - Différences finies appliquées à 
             l’équation de la chaleur bidimensionnelle

INSTRUCTIONS:
        How to compile: mpicc -o Final_Parallel_heat2D Final_Parallel_heat2D.c -lm
        How to execute: mpirun -np NB_Of_THREADS Final_Parallel_heat2D
            
            
Last update: 12 February 2018

This project is also available at Github
    >> https://goo.gl/23qq51 


*/


/// ---------------------------------- R E S U L T S  -------------------------------------///
/*
    For this TP I've worked in two versions. The first one (file Parallel_heat2D_WithoutSCATTER-GATTER.c)
    all the communications are based in Recv and Sendv. The data sharing was defined by me in
    execution time. When this one was working fine, I started this version using the functions
    Gather and Scatter. In every iteration all threads send the data for the thread 0 using Gather
    and the last and first line for each neighbors using the P2P functions (Recv and Send). Actually
    I used the ISend (assyncronous version) to send and the syncronous version to receive. I've
    made some tests and everything is ok. After some executions changing the size of the matrix
    and the number of threads I've the following results:

    NOTE: The number of iterations was fixed in 100. For the benchmask I didn't execute writing 
          each step on the out_put file. I did it because the writing procedure isn't parralelized
          and it blockk the execution time.

    ___________________________________________
   |               Matrix 500x500              |
    -------------------------------------------
   |Nb of threads| Cpu Time (s) | comm Time (s)|
    -------------------------------------------
   |      1      |     3.95     |     0.02     |
   |      4      |     1.14     |     0.06     |
   |      8      |     1.03     |     0.08     |
   |      12     |     0.83     |     0.09     |
   |      16     |     0.60     |     0.11     |
    -------------------------------------------


    ___________________________________________
   |              Matrix 1000x1000             |
    -------------------------------------------
   |Nb of threads| Cpu Time (s) | comm Time (s)|
    -------------------------------------------
   |      1      |     16.67    |     0.26     | 
   |      4      |      4.04    |     0.25     |
   |      8      |      2.74    |     0.28     |
   |      12     |      2.08    |     0.29     |
   |      16     |      2.49    |     0.31     |
    -------------------------------------------


    ___________________________________________
   |              Matrix 1500x1500             |
    -------------------------------------------
   |Nb of threads| Cpu Time (s) | comm Time (s)|
    -------------------------------------------
   |      1      |      19.06     |    0.37    |   
   |      4      |       5.77     |    0.39    |
   |      8      |       4.62     |    0.41    | 
   |      12     |       3.25     |    0.45    |     
   |      16     |       2.91     |    0.51    | 
    -------------------------------------------


    ___________________________________________
   |              Matrix 10000x10000           |
    -------------------------------------------
   |Nb of threads| Cpu Time (s) | comm Time (s)|
    -------------------------------------------
   |      1      |      277.02    |    4.75    |   
   |      4      |      90.77     |    4.98    |
   |      8      |      59.88     |    5.21    | 
   |      12     |      47.14     |    5.98    |     
   |      16     |      37.21     |    6.19    | 
    -------------------------------------------



*/
/// --------------------------- E N D   O F   R E S U L T S  -----------------------------///


//
// heat2D_Jacobi.c
//
// Diffusion de chaleur par la methode de Jacobi
//
// Utilisation : heat2D_Jacobi nbLigs nbCols
// les arguments sur la ligne de commande sont optionnels, 
// leurs valeurs sont =80 par defaut 
//

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define MAXVALUE 511  // la valeur maximale des donnees, utilisee dans initData()
#define NBSTEP   100    // le nombre d'iterations

/* Masque pour la mise a jour des donnees */
double MASK[] = {0.0, 0.1, 0.0, 0.1, 0.6, 0.1, 0.0, 0.1, 0.0};

/**
 * Mise a jour des lignes de start a end
 * @param :
 *      data : anciennes valeurs
 *      newData : nouvelles valeurs
 *      nbCols : nombre de donnees en Y
 */
void updateNLigs(double *old_data, double *new_data, int start, int end, int nbCols, double *masque){
    int ix, iy, ixm, iym;
    double *ptr_masque;
    
    for (ix=start; ix<=end; ix++) {                                 //lOOP FOR EACH LINE
        for (iy=1; iy<nbCols-1; iy++) {                             //LOOP FOR EACH COLUMN
            
            ptr_masque = masque;
            *(new_data+ix*nbCols+iy) = 0.0;                         //CLEAR THE CURRENT POSITION
            for (ixm=-1; ixm<=1; ixm++)                             //RUN AT THE MASK - LINE BY LINE
                for (iym=-1; iym<=1; iym++) {                       //RUN AT THE MASK - COLUMN BY COLUMN
                    *(new_data+ix*nbCols+iy) += *ptr_masque * *(old_data + (ix+ixm)*nbCols + iy+iym);
                    ptr_masque++;
                }
        }
    }
}


/**
 * Initialisation de donnees : un domaine de nbLigs x nbCols points
 */
double * initData(int nbLigs, int nbCols)
{
	double *data = NULL;
	int    i, j;

	data = (double *) malloc(nbLigs * nbCols * sizeof(double));

	if (data) {
		for (i=0; i<nbLigs; i++)
			for (j=0; j<nbCols; j++)
				*(data + i*nbCols + j) = ((double) (i * (nbLigs-i-1) * j * (nbCols-j-1)) / (pow(nbLigs/2.0,2.0) * pow(nbCols/2.0, 2.0))) * MAXVALUE;
	}

	return data;
}


/**
 * Savegarde des donnees dans un fichier
 * @param :
 *      data : pointeur de nbLigs x nbCols double
 *      fname : nom du fichier de sauvegarde
 */
void saveData(double *data, int nbLigs, int nbCols, char *fname)
{
	FILE *fp = fopen(fname, "w");
    int   i, j;

	if (fp) {
		fprintf(fp, "%d %d\n", nbLigs, nbCols);

		for (i=0; i<nbLigs; i++)
			for (j=0; j<nbCols; j++)
				fprintf(fp, "%6.1f%c", *(data+i*nbCols+j),(j==nbCols-1)?'\n':' ');
	    printf("\n");
	
		fclose(fp);
	}

	
}


/**
 * Ajout des donnees a un fichier
 * @param :
 *      data : pointeur de nbLigs x nbCols double
 *      fname : nom du fichier de sauvegarde
 */
void appendData(double *data, int nbLigs, int nbCols, char *fname)
{
    FILE *fp = fopen(fname, "a");
    int   i, j;
    
    if (fp) {
        fprintf(fp, "\n");
        
        for (i=0; i<nbLigs; i++)
            for (j=0; j<nbCols; j++)
                fprintf(fp, "%6.1f%c", *(data+i*nbCols+j),(j==nbCols-1)?'\n':' ');
        
        fclose(fp);
    }
    
}

/**
 * Echange de 2 pointers sur double
 */
void exchangeData(double **data1, double **data2)
{
    double *tmp = *data1;
    
    *data1 = *data2;
    *data2 = tmp;
}


int main(int argc, char **argv)
{
    //TimeStamp Variables 
    double PinitTime, PendTime, CommTime=0.0,CommInitTime,calcTime=0.0, calcInittime, writingTime=0.0, startWritingTime;

    //Commom variables
    int rang, nbProcs;
    int     nbLigs=500, nbCols=500;
    //originalMatrix used by thread 0
    double *startMatrix=NULL;

    //Local data for each thread
    double *data1=NULL, *data2=NULL;

    
    if (argc>=3) {
        nbLigs = atoi(argv[1]);
        nbCols = atoi(argv[2]);
    }
    

///////////////////////
// PARALLELL SECTION //
///////////////////////

    MPI_Status statut;
    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &rang );
    MPI_Comm_size( MPI_COMM_WORLD, &nbProcs );

    //Getting the initial time
    MPI_Barrier(MPI_COMM_WORLD);

 
   
    //Nb of lines per node
    int blocksize = (nbLigs-2) / nbProcs;

    //Aux variable used in the calcul of sendCounts
    int sum=0;

    int *sendcounts = malloc(sizeof(int)*nbProcs);  //Memory allocation - The block size of each thread 
    int *displs = malloc(sizeof(int)*nbProcs);      //Memory allocation - The starting point at the matrix
    int i;

    //Initialize the aux vectors
    for (i = 0; i < nbProcs; i++) {
        sendcounts[i] = (blocksize)*nbCols;
        displs[i] = sum;
        sum += sendcounts[i];
    }
    //The last thread could has more elements
    sendcounts[nbProcs-1]+=((nbLigs-2)%nbProcs)*nbCols;
    if(rang==nbProcs-1){
        blocksize+=(nbLigs-2)%nbProcs;
    }

    //Allocation of the auxilioary matrix for each  thread
    data1 = (double *) calloc(sendcounts[rang]+2*nbCols, sizeof(double));
    data2 = (double *) calloc(sendcounts[rang]+2*nbCols, sizeof(double));
   

    if(rang==0){            
// = = = = = = = = = = = = = = = = = = = = = = = = = //
// I'm the master and I have to distribute the data  //
// = = = = = = = = = = = = = = = = = = = = = = = = = //

        printf("Taille du probleme : %d x %d\n", nbLigs, nbCols);
        printf("Pour la modifier : %s nbLigs nbCols\n\n", *argv);

        //creating the first matrix
        startMatrix = initData(nbLigs, nbCols);
        PinitTime = MPI_Wtime();
    }

    MPI_Barrier(MPI_COMM_WORLD);
//Sending the information for everyone
    MPI_Scatterv(startMatrix+nbCols, sendcounts, displs, MPI_DOUBLE, data1+nbCols, sendcounts[rang], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    

//Starting the processing
    
    if (data1 && data2) {
        if(rang==0)
            saveData(startMatrix, nbLigs, nbCols, "saveData_Parallel_version");
        
        int step;
	    /* Calcul de la diffusion de la chaleur */
        for (step=0; step<  NBSTEP; step++) {
            //Share my information with my neighbor
            MPI_Request SendRequest, SendRequest2;
            
            MPI_Barrier(MPI_COMM_WORLD);
            if(rang==0){
                CommInitTime = MPI_Wtime();
            } 
            

            //Sending my data
            if(rang!=0){ 
                //Ok. I'm not the root, so I've to send my first line to everyone and receive the last line from everyone else
                 MPI_Isend(data1+nbCols, nbCols , MPI_DOUBLE,rang-1,0,MPI_COMM_WORLD,&SendRequest);
                 MPI_Recv(data1, nbCols , MPI_DOUBLE,rang-1,0,MPI_COMM_WORLD,&statut);
            }

            if(rang!=nbProcs-1 && nbProcs!=1){
                 //Ok, I'm not the last thread, so I've to send my last line and to receive the first one from everyone else
                 MPI_Recv(data1 + (blocksize+1)*nbCols, nbCols , MPI_DOUBLE,rang+1,0,MPI_COMM_WORLD,&statut);
                 MPI_Isend(data1 + (blocksize)*nbCols, nbCols , MPI_DOUBLE,rang+1,0,MPI_COMM_WORLD,&SendRequest2);
            }
            
            MPI_Barrier(MPI_COMM_WORLD);
            if(rang==0){
                CommTime +=(MPI_Wtime()-CommInitTime);
                calcInittime=MPI_Wtime();
            } 
            

            updateNLigs(data1, data2, 1, blocksize, nbCols, MASK);
            exchangeData(&data1, &data2);

            //This part gets the communication time
             MPI_Barrier(MPI_COMM_WORLD);
            if(rang==0){
                calcTime+=(MPI_Wtime() - calcInittime);
                CommInitTime = MPI_Wtime();
            }
            
            //Send my data to thread 0
            MPI_Gatherv(data1+nbCols, sendcounts[rang], MPI_DOUBLE, startMatrix+nbCols, sendcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
            MPI_Barrier(MPI_COMM_WORLD);
            if(rang==0){
                CommTime += (MPI_Wtime()-CommInitTime);
                
            } 
            
            // the thread 0 must save the data at the file
            if(rang==0){
                startWritingTime = MPI_Wtime();
                //appendData(startMatrix, nbLigs, nbCols, "saveData_Parallel_version");
                writingTime += (MPI_Wtime()-startWritingTime);
            }


        }
     }

    if (data2) free(data2);
    if (data1) free(data1);
    if (startMatrix) free(startMatrix);

    //Get the total execution time
    MPI_Barrier(MPI_COMM_WORLD);
    if(rang==0){
        PendTime = MPI_Wtime();
        printf("\tTotal CPU Time %f sec\n", (PendTime - PinitTime));
        printf("\tTotal Communication Time (considering the Gatter at each iteration)= %f sec\n", CommTime);
        printf("\tTotal Calculation Time = %f sec\n", calcTime);
        printf("\tTotal Writing Time = %f sec\n", writingTime);
           
    }

    //Finalize the MPI context
    MPI_Finalize();
   return EXIT_SUCCESS;
}
