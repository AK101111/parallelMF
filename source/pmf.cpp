#include "pmf.h"
#include "sched.h"
#include <iostream>
#include <omp.h>
#include <cstring>
#include <stdio.h>
#include <cstdlib>
#include <ctime>

/* TODO :
	Done - Create s threads
	Done - Initialize scheduler for s threads
	Parallelize this:
	DONE - pick a free block from scheduler
	do sgd on this block
	DONE - put block back in shceduler
	DONE - do for given number of iterations (instead do it till error becomes < eps?)
	return factors
*/

void read_matrix(char filename[], matrix_size *ms, float **R){
  FILE * fp; 
  fp = fopen(filename, "r");
  if(!fp){
    printf("Cannot open file %s\n", filename);
    exit(1);
  }

  fscanf(fp,"%d ",&(ms->row_size));
  fscanf(fp,"%d ",&(ms->col_size));

  R = new float*[ms->row_size];
	for(int i=0;i<ms->row_size;i++) R[i] = new float[ms->col_size];

  for(int i=0; i<ms->row_size; ++i){
  	for(int j=0; j<ms->col_size; ++j){
  		fscanf(fp,"%f ",&R[i][j]);
  	}
  }
 	fclose(fp); 
}

void dump_matrix(char filename[], float **R, matrix_size ms){
	FILE * fp;
  fp = fopen(filename,"w");
  if(!fp){
    printf("Cannot create the file %s\n", filename);
    exit(1);
  }
    
  for(int i=0; i<ms.row_size; ++i){
  	for(int j=0; j<ms.col_size; ++j){
  		fprintf(fp,"%f ",R[i][j]);
  	}
  	fprintf(fp,"\n");
  }   
 
  fclose(fp);
}

float calc_err(decomposition *dec, float **R, matrix_size ms){
	float ret = 0.0;
	for(int i=0;i<ms.row_size;++i){
		for(int j=0;j<ms.col_size;++j){
			// computing e_ij
			float error = 0.0;
			float **X = dec->X; float **Y = dec->Y;
			for(int k=0;k<dec->MX.col_size;++k){
				error += (X[i][k]*Y[k][j]);
			}
			error = R[i][j] - error;
			ret += error;
		}
	}	
	return ret;
}

// TODO see change in _launch_sched. Understand ojas's use of sub_rows, sub_cols.
// what is rows/cols_in_use	
int _randomizer(float **R, matrix_size MR, int **correspondance){
	return 0;
}

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

//each thread will do this separately
void _factorize_block(prob_params *params, float** R, matrix_size MR, \
	decomposition *dec, int iteration){
	if(iteration == params->num_iter){
		return;
	}
	block cur_block;
	
	#pragma omp critical 
	cur_block = get_block();
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
			for(int k=0;k<params->dim;++k){
				error_ij += (X[i + cur_block.x_index][k]*Y[k][cur_block.y_index + j]);
			}
			error_ij = R[i + cur_block.x_index][j + cur_block.y_index] - error_ij;
			iter_err += error_ij;
			// new X,Y

			float *X_new = new float[params->dim]();
			float *Y_new = new float[params->dim]();
			for(int k=0;k<params->dim;++k){
				X_new[k] = X[i+cur_block.x_index][k] + (params->lr * ( (error_ij*Y[k][j + cur_block.y_index]) - (params->lambda*X[i + cur_block.x_index][k]) ) );
				Y_new[k] = Y[k][j+cur_block.y_index] + (params->lr * ( (error_ij*X[i + cur_block.x_index][k]) - (params->mu*Y[k][j + cur_block.y_index]) ) );
			}
			// update X,Y
			for(int k=0;k<params->dim;++k){
				X[i + cur_block.x_index][k] = X_new[k];
				Y[k][j + cur_block.y_index] = Y_new[k];
			}
			//std::memcpy(X[i + cur_block.x_index], X_new, (params->dim)*sizeof(float));
			//std::memcpy(Y[j + cur_block.y_index], Y_new, (params->dim)*sizeof(float));
			free(X_new);
			free(Y_new);
		}
	}
	//#ifdef DEBUG
  //printf("error in iteration %d: %f\n", iteration, iter_err);
  //#endif

	#pragma omp critical
	push_block(cur_block);
	_factorize_block(params, R, MR, dec, iteration + 1);
}

void matrix_factorize(prob_params *params, float** R,\
  decomposition *dec, matrix_size MR){
	if(dec->MX.row_size != MR.row_size || dec->MY.col_size != MR.col_size || dec->MX.col_size != dec->MY.row_size){
		std::cout << "Matrix dimensions don't match. Exiting!";
		exit(1);
	}
	//omp_set_dynamic(0);
	//omp_set_num_threads(params->num_threads);
	_launch_sched(MR, params->num_threads);
	// parallelize this
	#pragma omp parallel num_threads(params->num_threads)
	for(int i=0;i<params->num_threads;i++){
		_factorize_block(params, R, MR, dec, 0);
	}
	return;
}

void random_data(decomposition *dec, float **R, matrix_size *mr){
	// fill dec->X
	if(dec->MX.col_size != dec->MY.row_size)
		return;
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
	}
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

int main(int argc, char* argv[]){
	srand (static_cast <unsigned> (time(0)));
	prob_params *params = (prob_params*)malloc(sizeof(prob_params));
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
	}
	#ifdef DEBUG
		printf("dim: %d lambda: %f mu: %f iter: %d threads: %d lr: %f\n", params->dim, params->lambda, params->mu, params->num_iter, params->num_threads, params->lr);	
	#endif
	
	matrix_size mx_s;
	float** R;
	read_matrix("data/10/10_10_1-1.R",&mx_s, R);
	//std::cin >> mx_s.row_size >> mx_s.col_size;
	decomposition *dec = (decomposition*)malloc(sizeof(decomposition));
	dec->X = new float*[mx_s.row_size];
	for(int i=0;i<mx_s.row_size;i++) dec->X[i] = new float[params->dim];
	dec->Y = new float*[params->dim];
	for(int i=0;i<params->dim;i++) dec->Y[i] = new float[mx_s.col_size];
	dec->MX.row_size = mx_s.row_size;
	dec->MX.col_size = params->dim;
	dec->MY.row_size = params->dim;
	dec->MY.col_size = mx_s.col_size;
	
	//random_data(dec, R, &mx_s);
	matrix_factorize(params, R, dec, mx_s);
	printf("final error: %f after %d iterations\n",calc_err(dec,R), iterations);
	//dump_matrix("1000_100.X",dec->X,dec->MX);
	//dump_matrix("100_1000.Y",dec->Y,dec->MY);

	return 0;
}