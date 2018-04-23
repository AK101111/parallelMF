#ifndef SCHED_H_
#define SCHED_H_

#include <stdlib.h>
#include <vector>
#include <queue>
#include "pmf.h"
#include "omp.h"

extern std::vector<int> rows_in_use;
extern std::vector<int> columns_in_use;

typedef struct block{
  int x_index;
  int y_index;
  matrix_size MB;
  int num_updates;
  double priority;
  //block();
  //block(int x_ind, int y_ind, matrix_size ms, int num_upd, double pr) :
  //            x_index(x_ind), y_index(y_ind), MB(ms), num_updates(num_upd), priority(pr) {}
}block;

// If something is free (1 instead of 0), return that block
// And among all free blocks, return the one with the earliest insertion time (priority)
struct lessThan{
  bool operator()(const block& lhs, const block& rhs) const{
  	return lhs.priority > rhs.priority;
  }
};

typedef std::priority_queue<block, std::vector<block>, lessThan> priorityQ;


extern priorityQ pq;

//void init_sched(int block_rows, int block_cols, int mat_rows, int mat_cols);
void init_sched(int mat_rows, int mat_cols);
block get_block();
void push_block(block b);
extern omp_lock_t wrlock;

#endif // SCHED_H_