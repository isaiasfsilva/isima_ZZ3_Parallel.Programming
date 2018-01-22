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
		How to compile: gcc -fopenmp calcMult.c -o calcMult
		How to execute: ./calcMult NBOFTHREADS SIZEMATRIX VERBOSE
			ps: NBOFTHREADS is optional. The defaut value is 3 threads.
				SIZEMATRIX  is optional. The defaut value is 3 - matrix 3x3
				VERBOSE (0 or 1) means VERBOSE MODE. It is optional. 
					The VERBOSE mode is disabled by defaut.

Last update: 15 january 2018

This project is also available at Github
	>> https://goo.gl/23qq51 
	
*/


#include <stdio.h>
#include "omp.h"
#include <unistd.h>


//Global variables
int NBofTHREADS = 3;
int N  = 3;
int VERBOSE=0;

//Global structures

typedef struct CKS CKS;
struct CKS{   //Checksun structure
	int *CSL; //Save the line checksun
	int *CSC; //Sabe the column checksun
};

//COMMON FUNCTIONS
int isEqual(CKS* a, CKS* b){
	int i;
	for(i=0;i<N;i++){
		if(a->CSL[i]!=b->CSL[i] || a->CSC[i]!=b->CSC[i])
			return 0;
	}
	return 1;
}

//////////////////////////////////////
//     Sequential multiplication    //
//           Question 3.1           //
//////////////////////////////////////

void question3dot1(int SLEEP, CKS* cks){

//local variables
	int a[N][N], b[N][N], res[N][N];
	int i,j,k,aux=0;
	
//starting A and B
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			a[i][j]=aux;
			b[i][j]=aux;
			aux++;
		}
	}

	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			aux = 0;
	        for (k = 0; k < N; k++) {
	          aux = aux + a[i][k]*b[k][j];
	        }	 
	        res[i][j] = aux;	          
		    usleep(SLEEP);           
      	}
	}

	if(VERBOSE){
		printf("Results of 3.1 \n");
	//printing results
		printf("\nMatrix A\n\t");
		for (i = 0; i < N; i++) {
			for (j = 0; j < N; j++)
				printf("%d\t", a[i][j]);
				printf("\n\t");
		}

	    printf("\nMatrix B\n\t");
		for (i = 0; i < N; i++) {
			for (j = 0; j < N; j++)
				printf("%d\t", b[i][j]);
				printf("\n\t");
		}

		printf("\nMatrix RES\n\t");
	    for (i = 0; i < N; i++) {
			for (j = 0; j < N; j++)
				printf("%d\t", res[i][j]);

				printf("\n\t");
		}
	}

	if(cks!=NULL){ // CALC THE CHECKSUN
		for(i=0;i<N;i++){
			for(j=0;j<N;j++){
				cks->CSL[i]+=res[i][j]*(j+1);
				cks->CSC[j]+=res[i][j]*(i+1); 
			}
		}
	}
}


//////////////////////////////////////
//      Parallel multiplication     //
//           Question 3.2           //
//              -------             //
// Here I divised the tasks between //
// my threads by hands. The question//
// 3.3 will be scheduled by the     //
// library openMP                   //
//                                  //
//////////////////////////////////////

void question3dot2(int SLEEP, CKS* cks){
//local variables
	int a[N][N], b[N][N], res[N][N];
	int i,j;

//Parallel Variables
	int nbItensByThread = N/NBofTHREADS;
	
//starting A and B
	#pragma omp parallel num_threads(NBofTHREADS)
	{
		//local variables
		int myid =  omp_get_thread_num();
		int base_ = myid*nbItensByThread;
		int upper = base_+nbItensByThread;
		upper +=(myid==(NBofTHREADS-1))?(N%NBofTHREADS):0;
		int aux_ =base_*N; ;
		int i_,j_;
		for (i_ = base_; i_ < upper; i_++) {
			for(j_=0;j_<N;j_++){
				a[i_][j_]=aux_;
				b[i_][j_]=aux_;
				aux_++;
			}
		}
	}

//Multipling
	#pragma omp parallel num_threads(NBofTHREADS)
	{
		//local variables
		int myid =  omp_get_thread_num();
		int base_ = myid*nbItensByThread;		
		int upper = base_+nbItensByThread;
		upper +=(myid==(NBofTHREADS-1))?(N%NBofTHREADS):0;
		int aux_ =base_*N; ;
		int i_,j_,k_;
		for (i_ = base_; i_ <upper; i_++) {
			for(j_=0;j_<N;j_++){
				aux_ = 0;
		        for (k_ = 0; k_ < N; k_++) {
		          aux_ = aux_ + a[i_][k_]*b[k_][j_];
		        }	 
		        res[i_][j_] = aux_; 
		        usleep(SLEEP);      
			}
		}
	}

//printing results
	if(VERBOSE){
		printf("Results of 3.2 \n");
	//printing results
		printf("\nMatrix A\n\t");
		for (i = 0; i < N; i++) {
			for (j = 0; j < N; j++)
				printf("%d\t", a[i][j]);
				printf("\n\t");
		}

	    printf("\nMatrix B\n\t");
		for (i = 0; i < N; i++) {
			for (j = 0; j < N; j++)
				printf("%d\t", b[i][j]);
				printf("\n\t");
		}

		printf("\nMatrix RES\n\t");
	    for (i = 0; i < N; i++) {
			for (j = 0; j < N; j++)
				printf("%d\t", res[i][j]);

				printf("\n\t");
		}
	}

	if(cks!=NULL){ // CALC THE CHECKSUN

		#pragma omp parallel num_threads(NBofTHREADS)
		{
			//local variables
			int myid =  omp_get_thread_num();
			int base_ = myid*nbItensByThread;		
			int upper = base_+nbItensByThread;
			upper +=(myid==(NBofTHREADS-1))?(N%NBofTHREADS):0;
			int aux_ =base_*N; ;
			int i_,j_,k_;
			for (i_ = base_; i_ <upper; i_++) {
				for(j_=0;j_<N;j_++){
					cks->CSL[i_]+=res[i_][j_]*(j_+1);
				}	
			}
		}

		#pragma omp parallel num_threads(NBofTHREADS)
		{
			//local variables
			int myid =  omp_get_thread_num();
			int base_ = myid*nbItensByThread;		
			int upper = base_+nbItensByThread;
			upper +=(myid==(NBofTHREADS-1))?(N%NBofTHREADS):0;
			int aux_ =base_*N; ;
			int i_,j_,k_;
			for (j_ = base_; j_ <upper; j_++) {
				for(i_=0;i_<N;i_++){
					cks->CSC[j_]+=res[i_][j_]*(i_+1); 
				}	
			}
		}




		
	}
}

//////////////////////////////////////
//      Parallel multiplication     //
//           Question 3.3           //
//              -------             //
// Here the tasks are automatically //
// divised between the threads using//
// the static and dynamic scheduling//
//                                  //
//////////////////////////////////////

void question3dot3(int MODE_DYNAMIC, int SLEEP, CKS* cks){
//local variables
	int a[N][N], b[N][N], res[N][N];
	int i,j;
	
//starting A and B
	if(MODE_DYNAMIC){
		#pragma omp parallel for schedule(static) num_threads(NBofTHREADS)
		for (i=0;i<N;i++){
			int aux_ =N*i;
			if(VERBOSE)
				printf("Task %d at thread %d\n", i, omp_get_thread_num());
			int j_;
			for(j_=0;j_<N;j_++){
				a[i][j_]=aux_;
				b[i][j_]=aux_;
				aux_++;
			}
		}
		
	
		#pragma omp parallel for schedule(static) num_threads(NBofTHREADS)
		for (i = 0; i < N; i++) {

			if(VERBOSE)
				printf("Task %d at thread %d\n", i, omp_get_thread_num());
			int j_, aux_;
			for (j_ = 0; j_ < N; j_++) {
				aux_ = 0;
				int k_;
		        for (k_ = 0; k_ < N; k_++) {
		          aux_ = aux_ + a[i][k_]*b[k_][j_];
		        }	 
		        res[i][j_] = aux_;   
		        usleep(SLEEP);        
	      	}
		}

	}else{ //Dynamic MODE
		#pragma omp parallel for schedule(dynamic) num_threads(NBofTHREADS)
		for (i=0;i<N;i++){
			int aux_ = N*i;

			if(VERBOSE)
				printf("Task %d at thread %d\n", i, omp_get_thread_num());
			int j_;
			for(j_=0;j_<N;j_++){
				a[i][j_]=aux_;
				b[i][j_]=aux_;
				aux_++;
			}
		}

		#pragma omp parallel for schedule(dynamic) num_threads(NBofTHREADS)
		for (i = 0; i < N; i++) {
			if(VERBOSE)
				printf("Task %d at thread %d\n", i, omp_get_thread_num());
			int j_,aux_;
			for (j_ = 0; j_ < N; j_++) {
				
				aux_ = 0;
				int k_;
		        for (k_ = 0; k_ < N; k_++) {
		          aux_ = aux_ + a[i][k_]*b[k_][j_];
		        }	 
		        res[i][j_] = aux_;  
		        usleep(SLEEP);   

	      	}
		}
	}


	if(VERBOSE){
		printf("Results of 3.3 \n");
	//printing results
		printf("\nMatrix A\n\t");
		for (i = 0; i < N; i++) {
			for (j = 0; j < N; j++)
				printf("%d\t", a[i][j]);
				printf("\n\t");
		}

	    printf("\nMatrix B\n\t");
		for (i = 0; i < N; i++) {
			for (j = 0; j < N; j++)
				printf("%d\t", b[i][j]);
				printf("\n\t");
		}

		printf("\nMatrix RES\n\t");
	    for (i = 0; i < N; i++) {
			for (j = 0; j < N; j++)
				printf("%d\t", res[i][j]);

				printf("\n\t");
		}
	}

	if(cks!=NULL){ // CALC THE CHECKSUN
		if(MODE_DYNAMIC){
			#pragma omp parallel for schedule(dynamic) num_threads(NBofTHREADS)
			for (i = 0; i < N; i++) {
				int j_;
				for(j_=0;j_<N;j_++){
					cks->CSL[i]+=res[i][j_]*(j_+1);					
				}	
			}

			#pragma omp parallel for schedule(dynamic) num_threads(NBofTHREADS)
			for (j = 0; j < N; j++) {
				int i_;
				for(i_=0;i_<N;i_++){					
					cks->CSC[j]+=res[i_][j]*(i_+1); 
				}	
			}


		}else{
			#pragma omp parallel for schedule(static) num_threads(NBofTHREADS)
			for (i = 0; i < N; i++) {
				int j_;
				for(j_=0;j_<N;j_++){
					cks->CSL[i]+=res[i][j_]*(j_+1);					
				}	
			}

			#pragma omp parallel for schedule(static) num_threads(NBofTHREADS)
			for (j = 0; j < N; j++) {
				int i_;
				for(i_=0;i_<N;i_++){					
					cks->CSC[j]+=res[i_][j]*(i_+1); 
				}	
			}

		}


	}

}


int main(int argc, char **argv){

//Local variables
	double PinitTime, PendTime;
	int i;

//Reading parameters
	if(argc > 1){// SET THE NB OF THREADS
		NBofTHREADS = atoi(argv[1]);   
	}

	if(argc > 2){// SET THE SIZE OF THE MATRICES
		N = atoi(argv[2]);   
	}

	if(argc > 3){// SET EXECUTION MODE 
		VERBOSE=atoi(argv[3]);   
	}
    
    

// ============= QUESTION 3.1 ===================
	PinitTime= omp_get_wtime();
	question3dot1(0, NULL);
	PendTime= omp_get_wtime();

	printf("\nQUESTION 3.1 = %f mlsec\n", (PendTime- PinitTime ));


// ============= QUESTION 3.2 ===================
	PinitTime= omp_get_wtime();
	question3dot2(0, NULL);
	PendTime= omp_get_wtime();

	printf("\nQUESTION 3.2 (schedule defined by me :) ) = %f mlsec\n", (PendTime- PinitTime ));


// ============= QUESTION 3.3 ===================
	// ======================================================== //
	//                    Question 3.3:                         //
	//                        --//--                            //
	// To verify that both calculs are identic just execute with//
	// the VERBOSE mode. Otherwise you can verify the checksum  //
	// verification too.                                        //
	// =========================================================//


// ============= QUESTION 3.4 ===================
	PinitTime= omp_get_wtime();
	question3dot3(0,0, NULL);
	PendTime= omp_get_wtime();

	printf("\nQUESTION 3.4 (schedule STATIC) = %f mlsec\n", (PendTime- PinitTime ));

	PinitTime= omp_get_wtime();
	question3dot3(1,0, NULL);
	PendTime= omp_get_wtime();

	printf("\nQUESTION 3.4 (schedule DYNAMIC) = %f mlsec\n", (PendTime- PinitTime ));

// ============= QUESTION 3.5 ===================
	// ======================================================== //
	//                    Question 3.5:                         //
	//                        --//--                            //
	// In my program there are 2 syncronizations. The first one //
	// after initializing the matrix before starting the calc of//
	// RES=AxB. I don't think that we can reduce this number    //
	// because to calc the first value of RES we need the first //
	// ligne of A and the first column of B. So, we can think   //
	// in some strategy but the performance will be almost the  //
	// the same.												//
	// =========================================================//

// ============= QUESTION 3.6 ===================
	// ======================================================== //
	//                    Question 3.6:                         //
	//                        --//--                            //
	// We can see that the size of the matrix and the number of //
	// threads have a high impact in the execution time. Of     //
	// if we increase the size of the matrix the time to calc   //
	// AxB will be bigger. The interesting result is in the     //
	// number of threads. If we increase until 16, each time    //
	// the performance is better. BUT if we set 50 threads we'll//
	// have a lot of "context switch" at the each CPU core. So  //
	// the time increases too beceuse these switchs are to much //
	// expensive.                                               //
	// I don't know why but using the schedule defined by me the//
	// calc time is so much better than defining schedule static//
	// or dynamic. May be for this case or this kind of calcul  //
	// it makes sense. Comparing static and dynamic scheduling  //
	// the dynamic scheduling is always better than the static  //
	// mode.                                                    //
	// =========================================================//


// ============= QUESTION 3.7 ===================
	// ======================================================== //
	//                    Question 3.7:                         //
	//                        --//--                            //
	// Every test has sleep SLEEP secondes after each calc of   //
	// each element of res. Now all iterations have to sleep    //
	// 0.01 seconds. This will increase the execution time      //
	// =========================================================//

	int SLEEP_time = 10000;

	printf("\nQUESTION 3.7 - Everyone is sleeping %f seconds after each calc\n", SLEEP_time);

	PinitTime= omp_get_wtime();
	question3dot1(SLEEP_time,NULL);
	PendTime= omp_get_wtime();

	printf("\nSEQUENTIAL = %f mlsec\n", (PendTime- PinitTime ));

	PinitTime= omp_get_wtime();
	question3dot2(SLEEP_time, NULL);
	PendTime= omp_get_wtime();

	printf("\tPARALLEL (schedule defined by me :) ) = %f mlsec\n", (PendTime- PinitTime ));

	PinitTime= omp_get_wtime();
	question3dot3(0,SLEEP_time, NULL);
	PendTime= omp_get_wtime();

	printf("\tPARALLEL (schedule STATIC) = %f mlsec\n", (PendTime- PinitTime ));

	PinitTime= omp_get_wtime();
	question3dot3(1,SLEEP_time, NULL);
	PendTime= omp_get_wtime();

	printf("\tPARALLEL (schedule DYNAMIC) = %f mlsec\n", (PendTime- PinitTime ));

// ============= QUESTION 3.8 ===================
	// ======================================================== //
	//                    Question 3.8:                         //
	//                        --//--                            //
	// Actually the parallel code is so much faster than the    //
	//sequential version. Sleeping 0.001 seconds with 16 threads//
	// and the matrices 100x100 the sequential time is 100 sec  //
	// and the parallel time is 6 with scheduling by myself, 3  //
	// with static scheduling and 8.24 with dynamic scheduling. //
	// Yes, it's so much better :)                              //
	// ======================================================== //

// ============= QUESTION 4.1 ===================
	
	printf("\nQUESTION 4.1 - All tests are beeing tested again... wait a minute...\n");
	CKS seq, pa_static, pa_dynamic, pa_byhand;
	seq.CSL       = (int*) calloc(N,sizeof(int));
	seq.CSC       = (int*) calloc(N,sizeof(int));	
	pa_static.CSL = (int*) calloc(N,sizeof(int));
	pa_static.CSC = (int*) calloc(N,sizeof(int));	
	pa_dynamic.CSL= (int*) calloc(N,sizeof(int));
	pa_dynamic.CSC= (int*) calloc(N,sizeof(int));	
	pa_byhand.CSL = (int*) calloc(N,sizeof(int));
	pa_byhand.CSC = (int*) calloc(N,sizeof(int));	
	
	
	question3dot1(SLEEP_time,&seq);
	question3dot2(SLEEP_time, &pa_byhand);
	question3dot3(0,SLEEP_time, &pa_static);
	question3dot3(1,SLEEP_time, &pa_dynamic);

	if(isEqual(&seq,&pa_static) && isEqual(&pa_static, &pa_dynamic) && isEqual(&pa_dynamic,&pa_byhand)){ 
		printf("\tVALID CHECKSUN\n");
	}else{
		printf("\tINVALID CHECKSUN\n");
	}

//  // ================ == TIME TABLE == =====================  //
	//  ______________________________________________________  //
	// |Table of results with 16 threads - SLEEP=0.03s - (sec)| //
	//  ------------------------------------------------------  //
	// |  MODE  |   5x5   | 10x10  |  20x20 | 50x50 | 100x100 | //
	//  ------------------------------------------------------  //
	// |   SEQ  |   0.25  |  1.00  |  4.03  | 25.25 |  101.04 | //
	// |PAR. STA|   0.05  |  0.10  |  0.41  | 2.018 |  7.072  | //
	// |PAR. DYN|   0.05  |  0.10  |  0.41  | 2.020 |  7.073  | //
	// |PAR. BYH|   0.25  |  1.00  |  1.01  | 2.525 |  10.10  | //
	//  ------------------------------------------------------| //
	//   Table 1 - Time results considering SLEEP as 0.03s.     //
	//             All times are in seconds.                    //
	// ======================================================== //


}