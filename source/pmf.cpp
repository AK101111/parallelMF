#include "pmf.h"
#include "sched.h"
#include <iostream>
#include <omp.h>
#include <cstring>

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

int _randomizer(float **R, matrix_size MR, unsigned int **correspondance){
	return 0;
}

// assume num_threads divides, MR.row_size, MR.col_size
void _launch_sched(matrix_size MR, unsigned int num_threads){
	init_sched();
	unsigned int sub_rows = MR.row_size / num_threads;
	unsigned int sub_cols = MR.col_size / num_threads;
	for(unsigned int i=0;i<sub_rows;i++){
		for(unsigned int j=0;j<sub_cols;j++){
			matrix_size mr = { .row_size = sub_rows, .col_size = sub_cols };
			block cur_block(i*MR.row_size, j*MR.col_size, mr, 0, 0 + rand()/(RAND_MAX + 1.0));
			push_block(cur_block);
		}
	}
	return;
}

//each thread will do this separately
void _factorize_block(prob_params *params, float** R, matrix_size MR, \
	decomposition *dec, unsigned int iteration){
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
	for(unsigned int i=0;i<cur_block.MB.row_size;++i){
		for(unsigned int j=0;j<cur_block.MB.col_size;++j){
			// computing e_ij
			float error_ij = 0.0;
			float **X = dec->X; float **Y = dec->Y;
			for(unsigned int k=0;k<params->dim;++k){
				error_ij += (X[i][k]*Y[k][j]);
			}
			error_ij = R[i][j] - error_ij;
			// new X,Y
			float *X_new = new float[params->dim]();
			float *Y_new = new float[params->dim]();
			for(unsigned int k=0;k<params->dim;++k){
				X_new[k] += (params->lr * ( (error_ij*Y[k][j]) - (params->lambda*X[i][k]) ) );
				Y_new[k] += (params->lr * ( (error_ij*X[i][k]) - (params->mu*Y[k][j]) ) );
			}
			// update X,Y
			std::memcpy(X[i], X_new, (params->dim)*sizeof(float));
			std::memcpy(Y[j], Y_new, (params->dim)*sizeof(float));
			free(X_new);
			free(Y_new);
		}
	}

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
	omp_set_dynamic(0);
	omp_set_num_threads(params->num_threads);

	_launch_sched(MR, params->num_threads);

	// parallelize this
	#pragma omp parallel num_threads(params->num_threads)
	for(unsigned int i=0;i<params->num_threads;i++){
		_factorize_block(params, R, MR, dec, 0);
	}
	return;
}

int main(int argc, char* argv[]){
	prob_params *params = (prob_params*)malloc(sizeof(prob_params));
	for(int i=1;i<argc;i++){
		if(strcmp(argv[i],"-dim")){
			params->dim = atoi(argv[i+1]);
			i += 1;
		}
		else if(strcmp(argv[i],"-lambda")){
			params->lambda = atof(argv[i+1]);
			i += 1;
		}
		else if(strcmp(argv[i],"-mu")){
			params->mu = atof(argv[i+1]);
			i += 1;
		}
		else if(strcmp(argv[i],"-iter")){
			params->num_iter = atoi(argv[i+1]);
			i += 1;
		}
		else if(strcmp(argv[i],"-threads")){
			params->num_threads = atoi(argv[i+1]);
			i += 1;
		}
		else if(strcmp(argv[i],"-lr")){
			params->lr = atof(argv[i+1]);
		}
	}
	matrix_size mx_s;
	std::cin >> mx_s.row_size >> mx_s.col_size;
	float** R = new float*[mx_s.row_size];
	for(unsigned int i=0;i<mx_s.row_size;i++) R[i] = new float[mx_s.col_size];
	decomposition *dec = (decomposition*)malloc(sizeof(dec));
	dec->X = new float*[mx_s.row_size];
	for(unsigned int i=0;i<mx_s.row_size;i++) dec->X[i] = new float[params->dim];
	dec->Y = new float*[params->dim];
	for(unsigned int i=0;i<params->dim;i++) dec->Y[i] = new float[mx_s.col_size];
	dec->MX.row_size = mx_s.row_size;
	dec->MX.col_size = params->dim;
	dec->MY.row_size = params->dim;
	dec->MY.col_size = mx_s.col_size;
	matrix_factorize(params, R, dec, mx_s);
	return 0;
}