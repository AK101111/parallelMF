# parallelMF
Parallel implementation of Matrix Factorization using OpenMP.
Given matrix R [m*n] and k>0, solves 
`R = XY`
where X [m*k] and Y[k*n]

---

## Requirements

This code requires openmp based c++ compiler.
We have tested the code on `crunchy` machines with gcc-6.4.0.

---

## Build


```bash
make
```

sample run:

```
./pmf -dim <dimension> -lambda <lambdaX> -mu <lambdaY> -iter <num-iter> -threads <threads for openmp> -lr <learning rate> -filename <datafile> -matSize <matrix size->m> <matrix size->n>
```

```
./pmf -dim 10 -lambda 0.01 -mu 0.01 -iter 10000 -threads 49 -lr 0.001 -filename data/1000/1000_1000_10-1.R -matSize 1000 1000
```

---

## Directory structure
The source code is contained in the `source/` dir. The matrix factorization code lives in `pmf.h` and the scheduler has been defined in the `sched.h`.

---

