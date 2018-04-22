#ifndef PMF_H_
#define PMF_H_

#include "stddef.h"

typedef struct matrix_size{
  int row_size;
  int col_size;
} matrix_size;

typedef struct decomposition{
	float **X;
  matrix_size MX;
  float **Y;
  matrix_size MY;
}decomposition;

typedef struct prob_params{
  float lambda;
  float mu;
  int num_iter;
  float iter;
  int dim;
  int num_threads;
  float lr;
}prob_params;

int _randomizer(float **R, matrix_size MR, int **correspondance);
void _launch_sched(matrix_size MR, int num_threads);
void _factorize_block(prob_params *params, float** R, matrix_size MR, \
  decomposition *dec, int iteration);
void matrix_factorize(prob_params *params, float** R,\
  decomposition *dec, matrix_size MR);


#endif // PMF_H_