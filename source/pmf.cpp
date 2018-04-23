//The name on top of each function is only of the person who wrote
//majority of that function and implemented the first version. 
//Both of us worked on the whole code eventually because of debugging, 
//integrating the code, and for error resolution.
#include "pmf.h"
#include "sched.h"
#include <iostream>
#include <omp.h>
#include <cstring>
#include <stdio.h>
#include <cstdlib>
#include <ctime>
#include <math.h>
#include <string>

/* TODO :
	Done - Create s threads
	Done - Initialize scheduler for s threads
	DONE - Parallelize this:
	DONE - pick a free block from scheduler
	DONE - do sgd on this block
	DONE - put block back in shceduler
	DONE - do for given number of iterations 
	DONE - TRIED & DISCARDED (instead do it till error becomes < eps?)
	DONE - return factors
*/

// Written by Arnav Kansal
void read_matrix(char filename[], matrix_size *ms, float **R){
	FILE * fp; 
	fp = fopen(filename, "r");
	if(!fp){
		printf("Cannot open file %s\n", filename);
		exit(1);
	}

	for(int i=0; i<ms->row_size; ++i){
		for(int j=0; j<ms->col_size; ++j){
			fscanf(fp,"%f ",&R[i][j]);
		}
	}
	fclose(fp); 
	return;
}

// Written by Arnav Kansal
void dump_matrix(char filename[], float **R, matrix_size ms){
	FILE * fp;
  fp = fopen(filename,"w");
  if(!fp){
    printf("Cannot create the file %s\n", filename);
    exit(1);
  }
  fprintf(fp,"%d %d\n ",ms.row_size, ms.col_size);
  for(int i=0; i<ms.row_size; ++i){
  	for(int j=0; j<ms.col_size; ++j){
  		fprintf(fp,"%f ",R[i][j]);
  	}
  	fprintf(fp,"\n");
  }   
 
  fclose(fp);
}

// Written by Arnav Kansal
float calc_err(decomposition *dec, float **R, matrix_size ms){
	float ret = 0.0;
	for(int i=0;i<ms.row_size;++i){
		for(int j=0;j<ms.col_size;++j){
			float error = 0.0;
			float **X = dec->X; float **Y = dec->Y;
			for(int k=0;k<dec->MX.col_size;++k){
				error += (X[i][k]*Y[j][k]);
			}
			error = R[i][j] - error;
			ret += error;
		}
	}	
	return ret;
}

// Written by Ojas Deshpande
// assume num_threads divides, MR.row_size, MR.col_size
void _launch_sched(matrix_size MR, int num_threads){
	init_sched(MR.row_size, MR.col_size);
	int sub_rows = MR.row_size / (num_threads+1);
	int sub_cols = MR.col_size / (num_threads+1);
	for(int i=0;i<num_threads+1;i++){
		for(int j=0;j<num_threads+1;j++){
			matrix_size mr = { .row_size = sub_rows, .col_size = sub_cols };
			block cur_block = {
				.x_index = i*sub_rows, 
				.y_index = j*sub_cols,
				.MB = mr, 
				.num_updates = 0, 
				.priority = (0 + rand()/(RAND_MAX + 1.0))
			};
			push_block(cur_block);
		}
	}
	return;
}	

// Written by Arnav Kansal
//each thread will do this separately
void _factorize_block(prob_params *params, float** R, matrix_size MR, \
	decomposition *dec, int iteration){
	while(iteration < params->num_iter){
		block cur_block;

		#pragma omp critical
		{
			cur_block = get_block();
		}
		//do sgd on this block
		/*
	 	*	e_ij = R_ij - XiYj
	 	*  X_i <- X_i + lr *(e_ij Yj - lambda Xi)
	 	*  Y_j <- Y_j + lr *(e_ij Xi - mu Yj)
	 	*/
	
		float iter_err = 0.0;
		//printf("%d %d\n", cur_block.MB.row_size, cur_block.MB.col_size);
		for(int i=0;i<cur_block.MB.row_size;++i){
			for(int j=0;j<cur_block.MB.col_size;++j){
				// computing e_ij
				float error_ij = 0.0;
				float **X = dec->X; float **Y = dec->Y;
				//#pragma omp parallel for
				for(int k=0;k<params->dim;++k){
					error_ij += (X[i + cur_block.x_index][k]*Y[cur_block.y_index + j][k]);
				}
				error_ij = R[i + cur_block.x_index][j + cur_block.y_index] - error_ij;
				iter_err += error_ij;
				if(fabs(error_ij) > 0.001){
					// new X,Y
					float *X_new = new float[params->dim]();
					float *Y_new = new float[params->dim]();
					//#pragma omp parallel for
					for(int k=0;k<params->dim;++k){
						X_new[k] = X[i+cur_block.x_index][k] + (params->lr * ( (error_ij*Y[j + cur_block.y_index][k]) - (params->lambda*X[i + cur_block.x_index][k]) ) );
						Y_new[k] = Y[j+cur_block.y_index][k] + (params->lr * ( (error_ij*X[i + cur_block.x_index][k]) - (params->mu*Y[j + cur_block.y_index][k]) ) );
					}
					// update X,Y
					std::memcpy(X[i + cur_block.x_index], X_new, (params->dim)*sizeof(float));
					std::memcpy(Y[j + cur_block.y_index], Y_new, (params->dim)*sizeof(float));
					//#pragma omp parallel for
					//for(int k=0;k<params->dim;++k){
					//	Y[k][j + cur_block.y_index] = Y_new[k];
					//}
					free(X_new);
					free(Y_new);
				}
			}
		}
		#ifdef DEBUG
	  		printf("error in iteration %d: %f\n", iteration, iter_err);
		#endif
		#pragma omp critical
		{
			push_block(cur_block);
		}
		iteration++;
	}
}

// Written by Ojas Deshpande
void matrix_factorize(prob_params *params, float** R,\
  decomposition *dec, matrix_size MR){
	if(dec->MX.row_size != MR.row_size || dec->MY.row_size != MR.col_size || dec->MX.col_size != dec->MY.col_size){
		std::cout << "Matrix dimensions don't match. Exiting!";
		exit(1);
	}
	omp_set_dynamic(0);
	omp_set_num_threads(params->num_threads);
	omp_init_lock(&wrlock);
	_launch_sched(MR, params->num_threads);
	#pragma omp parallel for num_threads(params->num_threads)
	for(int i=0;i<params->num_threads;i++){
		_factorize_block(params, R, MR, dec, 0);
	}
	return;
}

// Written by Arnav Kansal
void random_data(decomposition *dec, float **R, matrix_size *mr){
	// fill dec->X
	if(dec->MX.col_size != dec->MY.col_size)
		return;
	/*for(int i=0; i<dec->MX.row_size; ++i){
		for(int j=0; j<dec->MX.col_size; ++j){
			dec->X[i][j] = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
		}
	}
	// fill dec->Y
	for(int i=0; i<dec->MY.row_size; ++i){
		for(int j=0; j<dec->MY.col_size; ++j){
			dec->Y[i][j] = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
		}
	}
	// calculate R
	mr->row_size = dec->MX.row_size;
	mr->col_size = dec->MY.col_size;
	for(int i=0; i<mr->row_size; ++i){
		for(int j=0; j<mr->col_size; ++j){
			R[i][j] = 0.0;
			for(int k=0; k<dec->MX.col_size; ++k){
				R[i][j] += (dec->X[i][k] * dec->Y[k][j]);
			}
		}
	}*/
	for(int i=0; i<dec->MX.row_size; ++i){
		for(int j=0; j<dec->MX.col_size; ++j){
			dec->X[i][j] = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
		}
	}
	// fill dec->Y
	for(int i=0; i<dec->MY.row_size; ++i){
		for(int j=0; j<dec->MY.col_size; ++j){
			dec->Y[i][j] = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
		}
	}
}

// Written by Ojas Deshpande
int main(int argc, char* argv[]){
	srand (static_cast <unsigned> (time(0)));
	prob_params *params = (prob_params*)malloc(sizeof(prob_params));
	char* fileName;
	matrix_size mx_s;
	for(int i=1;i<argc;i++){
		if(!strcmp(argv[i],"-dim")){
			params->dim = atoi(argv[i+1]);
		}
		else if(!strcmp(argv[i],"-lambda")){
			params->lambda = atof(argv[i+1]);
		}
		else if(!strcmp(argv[i],"-mu")){
			params->mu = atof(argv[i+1]);
		}
		else if(!strcmp(argv[i],"-iter")){
			params->num_iter = atoi(argv[i+1]);
		}
		else if(!strcmp(argv[i],"-threads")){
			params->num_threads = atoi(argv[i+1]);
		}
		else if(!strcmp(argv[i],"-lr")){
			params->lr = atof(argv[i+1]);
		}
		else if(!strcmp(argv[i], "-filename")){
			fileName = argv[i+1];
		}
		else if(!strcmp(argv[i], "-matSize")){
			mx_s.row_size = atoi(argv[i+1]);
			mx_s.col_size = atoi(argv[i+2]);
		}
	}
	#ifdef DEBUG
		printf("dim: %d lambda: %f mu: %f iter: %d threads: %d lr: %f\n", params->dim, params->lambda, params->mu, params->num_iter, params->num_threads, params->lr);	
	#endif
	float** R = new float*[mx_s.row_size];
	for(int i=0;i<mx_s.row_size;i++) R[i] = new float[mx_s.col_size];
	read_matrix(fileName, &mx_s, R);
 	decomposition *dec = (decomposition*)malloc(sizeof(decomposition));
	dec->X = new float*[mx_s.row_size];
	for(int i=0;i<mx_s.row_size;i++) dec->X[i] = new float[params->dim];
	dec->Y = new float*[mx_s.col_size];
	for(int i=0;i<mx_s.col_size;i++) dec->Y[i] = new float[params->dim];

	dec->MX.row_size = mx_s.row_size;
	dec->MX.col_size = params->dim;
	dec->MY.row_size = mx_s.col_size;
	dec->MY.col_size = params->dim;
	
	
	random_data(dec, R, &mx_s);
	matrix_factorize(params, R, dec, mx_s);
	printf("final error: %f after %d iterations\n",calc_err(dec,R,mx_s), params->num_iter);
	dump_matrix("X.gen",dec->X,dec->MX);
	dump_matrix("Y.gen",dec->Y,dec->MY);
	return 0;
}
