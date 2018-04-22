#include "sched.h"

priorityQ pq;
std::vector<int> rows_in_use;
std::vector<int> columns_in_use;

//void init_sched(int block_rows, int block_cols, int mat_rows, int mat_cols){
void init_sched(int mat_rows, int mat_cols){
	rows_in_use.resize(mat_rows);
	columns_in_use.resize(mat_cols);
	std::fill(rows_in_use.begin(), rows_in_use.end(), 0);
	std::fill(columns_in_use.begin(), columns_in_use.end(), 0);
	return;
}

// insert blocks
void push_block(block b){
	//reset the number of updates for the block.
	b.num_updates += 1;
	// for all blocks with same number of updates, we want to select one randomly.
	b.priority = b.num_updates + rand()/(RAND_MAX + 1.0);
	pq.push(b);
	rows_in_use[b.x_index] = 0;
	columns_in_use[b.y_index] = 0;
	return;
}

// TODO fix corner case. Program might loop forever/ what if pq is empty??
block get_block(){
	std::vector<block> tempList;
	// return the block with smallest number of updates, and which is also free.
	while(!pq.empty()){
		block b = pq.top();
		pq.pop();
		if(rows_in_use[b.x_index] == 0 && columns_in_use[b.y_index] == 0){
			for(std::vector<block>::iterator it = tempList.begin();it!=tempList.end();it++){
				pq.push(*it);
			}
			rows_in_use[b.x_index] = 1;
			columns_in_use[b.y_index] = 1;
			return b;
		}
		tempList.push_back(b);
	}
}