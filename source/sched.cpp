#include "sched.h"
#include "omp.h"
#include <iostream>

priorityQ pq;
std::vector<int> rows_in_use;
std::vector<int> columns_in_use;

//void init_sched(int block_rows, int block_cols, int mat_rows, int mat_cols){
void init_sched(int mat_rows, int mat_cols){
	rows_in_use.resize(mat_rows);
	columns_in_use.resize(mat_cols);
	std::fill(rows_in_use.begin(), rows_in_use.end(), 0);
	std::fill(columns_in_use.begin(), columns_in_use.end(), 0);
	omp_init_lock(&writelock);
	return;
}

// insert blocks
void push_block(block b){
	//reset the number of updates for the block.
	b.num_updates += 1;
	// for all blocks with same number of updates, we want to select one randomly.
	b.priority = b.num_updates + rand()/(RAND_MAX + 1.0);
	omp_set_lock(&writelock);
	pq.push(b);
	omp_unset_lock(&writelock);
	rows_in_use[b.x_index] = 0;
	columns_in_use[b.y_index] = 0;
	return;
}

// TODO fix corner case. Program might loop forever/ what if pq is empty??
block get_block(){
	std::vector<block> tempList;
	// return the block with smallest number of updates, and which is also free.
	while(!pq.empty()){
		omp_set_lock(&writelock);
		block b = pq.top();
		pq.pop();
		omp_unset_lock(&writelock);
		if(rows_in_use[b.x_index] == 0 && columns_in_use[b.y_index] == 0){
			omp_set_lock(&writelock);
			for(std::vector<block>::iterator it = tempList.begin();it!=tempList.end();it++){
				pq.push(*it);
			}
			omp_unset_lock(&writelock);
			rows_in_use[b.x_index] = 1;
			columns_in_use[b.y_index] = 1;
			return b;
		}
		tempList.push_back(b);
	}
}