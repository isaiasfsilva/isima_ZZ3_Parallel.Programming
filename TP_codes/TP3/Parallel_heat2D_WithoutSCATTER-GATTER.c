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
        This is a Parallell version of the TP3 - Part 1
             Communication locale / globale - Différences finies appliquées à 
             l’équation de la chaleur bidimensionnelle

INSTRUCTIONS:
        How to compile: mpicc –o Parallel_heat2D_WithoutSCATTER-GATTER Parallel_heat2D_WithoutSCATTER-GATTER.c -lm
        How to execute: mpirun -np NB_Of_THREADS Parallel_heat2D_WithoutSCATTER-GATTER
            
            
Last update: 30 january 2018

This project is also available at Github
    >> https://goo.gl/23qq51 


 //Update: 06 February 2018
 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   PLEASE READMEE   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++   EVERYTHING IS OK BUT IF WE HAVE 4 NODES AND MATRIX SIZE 10X10, WE HAVE A ERROR; ++++++++++++++++++++++++++++++++++++++++++++++++++++++++

 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   IS SMETHING WITH INDEXATION   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


//Update: 11 February 2018

++++++++++++++++++++++++++++++ I stopped this version because the subjet of this TP was changed. When I got the first version 1 month ago wasn't necessary to use Scatter/Gather ++++++++++++++++++

 //                Now i'm starting a new version using the MPI_Scatter and MPI_Gather. Check the file Final_Parallel_heat2D.c to verify the correct and final version of this TP.

*/


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


#define MAXVALUE 5  // la valeur maximale des donnees, utilisee dans initData()
#define NBSTEP   5    // le nombre d'iterations

/* Masque pour la mise a jour des donnees */
double MASK[] = {0.0, 0.5, 0.5, 0.5, 0.6, 0.5, 0.5, 0.7, 0.0};

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

    //Commom variables
    int rang, nbProcs;

    double *startMatrix=NULL;

    double *data1=NULL, *data2=NULL;
    int     nbLigs=10, nbCols=10;
    
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

    //Nb of lines per node
    int blocksize = (nbLigs-2) / nbProcs;

    //Local variables
    int startLine = 1 + blocksize*rang;
    int endLine   = startLine + blocksize;
    int step      =0;

    //Just the last one
    if(rang==nbProcs-1){
        endLine+=(nbLigs-2)%nbProcs;
    }

    data1 = (double *) calloc((endLine - startLine +2)*nbCols, sizeof(double));
    data2 = (double *) calloc((endLine - startLine +2)*nbCols, sizeof(double));

    if(rang==0){            
// = = = = = = = = = = = = = = = = = = = = = = = = = //
// I'm the master and I have to distribute the data  //
// = = = = = = = = = = = = = = = = = = = = = = = = = //

        printf("Taille du probleme : %d x %d\n", nbLigs, nbCols);
        printf("Pour la modifier : %s nbLigs nbCols\n\n", *argv);

        //creating the first matrix
        startMatrix = initData(nbLigs, nbCols);


        //DEBUG
        //printf("The original matrix is:\n"); 
        int j;
        int i;
        for (i=0; i<nbLigs; i++){
            for (j=0; j<nbCols; j++){
               // printf("%6.1f%c", *(startMatrix+i*nbCols+j),(j==nbCols-1)?'\n':' ');
            }
        }
      
        




        //Send for everyone
        //int i;
        for (i = 1; i < nbProcs-1; ++i){
            //printf("Sending to the (%d) %d lines [starting in %d]\n",i,(blocksize+2),(i*blocksize));
            MPI_Send(startMatrix + (i*blocksize)*nbCols, (blocksize+2)*nbCols , MPI_DOUBLE,i,0,MPI_COMM_WORLD);
        }


        //The last node may have more work to do, so its data slide could be bigger
        if(nbProcs!=1){
           // printf("Sending to the last guy %d lines [%d - %d]\n",(blocksize+2) + (nbLigs-2)%nbProcs,(nbProcs-1)*blocksize,(blocksize+2) + (nbLigs-2)%nbProcs);
            MPI_Send(startMatrix + (nbProcs-1)*blocksize*nbCols, ((blocksize+2) + (nbLigs-2)%nbProcs)*nbCols  , MPI_DOUBLE,i,0,MPI_COMM_WORLD);
        }

        //Copy my part            
        memcpy ( data1, startMatrix, sizeof(double)*(blocksize+2)*nbCols );

        if (startMatrix) free(startMatrix);



    }else{
        //Receiving my data slice to start ly job
        // printf("I'm %d and I'm Receiving %d lines:\n", rang,endLine-startLine+2); 
        MPI_Recv(data1,(endLine-startLine+2)*nbCols , MPI_DOUBLE,0,0,MPI_COMM_WORLD,&statut);
    }
    
//Starting the processing
    
    if (data1 && data2) {


     
        //saveData(data1, nbLigs, nbCols, "saveData");
        

	    /* Calcul de la diffusion de la chaleur */
        for (step=0; step<  2; step++) {

      
            

            updateNLigs(data1, data2, 1, endLine-startLine, nbCols, MASK);

            exchangeData(&data1, &data2);



            MPI_Request SendRequest, SendRequest2;
            //Sending my data
            if(rang!=0){

                 MPI_Isend(data1+nbCols, nbCols , MPI_DOUBLE,rang-1,0,MPI_COMM_WORLD,&SendRequest);
                 MPI_Recv(data1, nbCols , MPI_DOUBLE,rang-1,0,MPI_COMM_WORLD,&statut);
            }
            if(rang!=nbProcs-1 && nbProcs!=1){
               //  printf("I'm %d and I'm sending the line %d to my successor \n", rang,(endLine-startLine));
                 MPI_Recv(data1 + (endLine-startLine+1)*nbCols, nbCols , MPI_DOUBLE,rang+1,0,MPI_COMM_WORLD,&statut);

                 MPI_Isend(data1 + (endLine-startLine)*nbCols, nbCols , MPI_DOUBLE,rang+1,0,MPI_COMM_WORLD,&SendRequest);
            }

        }
     }
    

    //Sending everything to root
     if(rang!=0){
         MPI_Send(data1 + nbCols, (endLine-startLine)*nbCols , MPI_DOUBLE,0,0,MPI_COMM_WORLD);
     }else{
        int i;
        for(i=1; i<nbProcs-1;i++){
            MPI_Recv(startMatrix + (i*(endLine-startLine)+1)*nbCols,(endLine-startLine)*nbCols , MPI_DOUBLE,i,0,MPI_COMM_WORLD,&statut);
        }
        //Receive from the last one
        if(nbProcs!=1)
            MPI_Recv(startMatrix + ((nbProcs-1)*blocksize+1)*nbCols, ((blocksize) + (nbLigs-2)%nbProcs)*nbCols  , MPI_DOUBLE,(nbProcs-1),0,MPI_COMM_WORLD,&statut);

        memcpy ( startMatrix, data1, sizeof(double)*(blocksize+1)*nbCols );


        printf("The final matrix is:\n"); 
        int j;
        for (i=0; i<nbLigs; i++){
            for (j=0; j<nbCols; j++){
                printf("%6.1f%c", *(startMatrix+i*nbCols+j),(j==nbCols-1)?'\n':' ');
            }
        }
     }

    if (data2) free(data2);
    if (data1) free(data1);

    MPI_Finalize();
   return EXIT_SUCCESS;
}
