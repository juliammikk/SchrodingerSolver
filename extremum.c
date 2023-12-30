#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "derivatives.h"
#include "solve.h"
#include "extremum.h"

typedef double (*FuncPT)(double);

FuncPT original_func;
double Extremum_DF(double x);

double Extremum_getExtremum(FuncPT func, double x_init, double *curvature){
    double tol, df, extremum_x;
    int count;
    original_func = func;

    count = 0;
    tol = 1.0e-10;
    extremum_x = Solve_Newton(0.0, Extremum_DF, x_init, tol, &count);
    *curvature = Derivative_SecondD(extremum_x, func);
    return extremum_x;
}

double Extremum_DF(double x){
    //double f;
    return Derivative_FirstD(x, original_func); //is x here correct
}