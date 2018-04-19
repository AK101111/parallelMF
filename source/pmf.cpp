#include "pmf.h"
#include "sched.h"
#include <iostream>
#include <omp.h>

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

}

void _launch_sched(matrix_size MR, unsigned int num_threads){
	init_sched();
	int sub_rows = MR.row_size / num_threads;
	int sub_cols = MR.col_size / num_threads;
	for(int i=0;i<sub_rows;i++){
		for(int j=0;j<sub_cols;j++){
			block cur_block(i*MR.row_size, j*MR.col_size, MR, 0, 0 + rand()/(RAND_MAX + 1.0));
			push_block(cur_block);
		}
	}
	return;
}

//each thread will do this separately
void do_stuff(prob_params *params, float** R, matrix_size MR, unsigned int iteration){
	if(iteration == params->num_iter){
		return;
	}
	block cur_block;
	
	#pragma omp critical
	cur_block = get_block();
	//do sgd on this block

	#pragma omp critical
	push_block(cur_block);
	do_stuff(params, R, MR, iteration + 1);
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
	for(int i=0;i<params->num_threads;i++){
		do_stuff(params, R, MR, 0);
	}
	return;
}

int main(int argc, char* argv[]){
	prob_params *params;
	for(int i=1;i<argc;i++){
		if(argv[i] == "-dim"){
			params->dim = atoi(argv[i+1]);
			i += 1;
		}
		else if(argv[i] == "-lambda"){
			params->lambda = atoi(argv[i+1]);
			i += 1;
		}
		else if(argv[i] == "-mu"){
			params->mu = atoi(argv[i+1]);
			i += 1;
		}
		else if(argv[i] == "-iter"){
			params->num_iter = atoi(argv[i+1]);
			i += 1;
		}
		else if(argv[i] == "-threads"){
			params->num_threads = atoi(argv[i+1]);
			i += 1;
		}
	}
	matrix_size mx_s;
	std::cin >> mx_s.row_size >> mx_s.col_size;
	float** R = new float*[mx_s.row_size];
	for(int i=0;i<mx_s.row_size;i++) R[i] = new float[mx_s.col_size];
	decomposition *dec;
	dec->X = new float*[mx_s.row_size];
	for(int i=0;i<mx_s.row_size;i++) dec->X[i] = new float[params->dim];
	dec->Y = new float*[params->dim];
	for(int i=0;i<params->dim;i++) dec->Y[i] = new float[mx_s.col_size];
	dec->MX.row_size = mx_s.row_size;
	dec->MX.col_size = params->dim;
	dec->MY.row_size = params->dim;
	dec->MY.col_size = mx_s.col_size;
	matrix_factorize(params, R, dec, mx_s);
	return 0;
}