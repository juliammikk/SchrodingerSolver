#include <stdio.h>
#include <stdlib.h>
#include "vector_mtx.h"

// allocate memory space for a 1D array
double *vector_malloc(int nmax)
{
    double *pt;
    int n;
    pt = (double *)malloc(sizeof(double)*nmax);

    // initialize all entries to 0
    for(n=0; n<nmax; n++) pt[n] = 0.0;
    return pt;
}

// allocate memory space for a 2D array
double **mtx_malloc(int mmax, int nmax)
{
    double **pt;
    int m, n;
    pt = (double **)malloc(sizeof(double *)*mmax);
    for(m=0; m<mmax; m++){
        pt[m] = (double *)malloc(sizeof(double)*nmax);
    }

    // initialize all entries to 0
    for(m=0; m<mmax; m++){
        for(n=0; n<nmax; n++) {pt[m][n] = 0.0;}
    }
    return pt;
}

// free the memory space allocated by mtx_malloc
void mtx_free(double **mtx, int mmax)
{
    int m;
    for(m=0; m<mmax; m++) {
        free(mtx[m]); }
    free(mtx);
}