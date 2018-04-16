#include "pmf.h"
#include "sched.h"
#include <iostream>

/* TODO :
	Create s threads
	Done - Initialize scheduler for s threads
	Parallelize this:
		get a free thread
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
	block cur_block = get_block();
	//do sgd on this block

	push_block(cur_block);
	do_stuff(params, R, MR, iteration + 1);
}

void matrix_factorize(prob_params *params, float** R,\
  decomposition *dec, matrix_size MR){
	if(dec->MX.row_size != MR.row_size || dec->MY.col_size != MR.col_size || dec->MX.col_size != dec->MY.row_size){
		std::cout << "Matrix dimensions don't match. Exiting!";
		exit(1);
	}
	_launch_sched(MR, params->num_threads);

	// parallelize this
	for(int i=0;i<params->num_threads;i++){
		do_stuff(params, R, MR, 0);
	}
	return;
}

