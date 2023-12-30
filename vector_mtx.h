#ifndef VECTOR_MTX_H
#define VECTOR_MTX_H

double *vector_malloc(int nmax); //allocates memory for 1d array
double **mtx_malloc(int mmax, int nmax); //allocates memory for 2d array
void mtx_free(double **mtx, int mmax);

#endif
