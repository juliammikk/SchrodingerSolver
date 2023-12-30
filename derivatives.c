#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "derivatives.h"

double Derivative_FirstD(double x, double (*func)(double)){
    double df, h;
    h = 1.05e-5;

    if(x!=0.0) h = h*x;

    df = ((*func)(x+h) - (*func)(x-h))/(2.0*h);
    return df;
}

double Derivative_SecondD(double x, double (*func)(double)){
    double dff, h;
    h = 1.05e-5;

    if(x!=0.0) h = h*x;

    dff = ((*func)(x+h) + (*func)(x-h) - 2*(*func)(x))/(h*h);
    return dff;
}