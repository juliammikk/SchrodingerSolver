#ifndef NUMEROV_PARAMS_H
#define NUMEROV_PARAMS_H
#include "params.h"

//define new type Func_1D to be a pointer to function that take in two args and return a double
typedef double (*Func_1D)(double, DynamicVars*);

typedef struct numerov_params{
    double x_f; //max x
    double x_i; //min x
    double y_0; //y(x_i)
    double y_1; //y(x_f)
    int nmax; //number of sampling points
    double h; //step size
    Func_1D NumerovFunc_F; //y"=Fy

} NumerovParams;
#endif