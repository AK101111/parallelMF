#ifndef PMF_H_
#define PMF_H_

typedef struct matrix_size{
  size_t row_size;
  size_t col_size;
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
  float num_iter;
  float iter;
  size_t dim;
  unsigned int num_threads;
}prob_params;

int _randomizer(float **R, matrix_size MR, unsigned int **correspondance);
int _launch_sched(matrix_size MR, unsigned int num_threads);
int matrix_factorize(prob_params *params, float** R,\
  decomposition *dec, matrix_size MR);

#endif // PMF_H_