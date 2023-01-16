//
// Created by megan on 14/01/2023.
//

// File: solve.c
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "solve.h"

// Numerical derivative
// This function is visible only within this file.
double Solve_Get_Df(double (*func)(double), double x);
// We could implement this numerical derivative as a general utility function
// and in fact we will do that shortly. We are including this here to make
// Solve function self-contained. But that is not strictly necessry.





// Solve f(x) = nu using bisect method
// Note that count is passed by reference (that is, the pointer is passed)
// because we want to change the value of it inside Solve_Bisect
double Solve_Bisect(double nu, double (*func)(double), double x_min, double x_max, double tol, int *count){
    double x_mid, f_max, f_min, f_mid, err;
    int count_max = 1000; // Large enough.
    *count += 1; // add 1 whenever BisectSolve is called

    // Warn and exit
    if(*count > count_max)
    {
        fprintf(stderr, "Solve_Bisect: Done %d iterations without convergence.\n",
                 count_max);
        fprintf(stderr, "Exiting.\n");
        exit(0);
    }


    f_max = (*func)(x_max) - nu; // Calculate f_max = f(x_max) - nu
    f_min = (*func)(x_min) - nu; //Calculate f_min = f(x_min) - nu

    if(f_max*f_min > 0.0) // we can’t find a solution within the range
    {
        //Warn and exit
        fprintf(stderr, "Solve_Bisect: Solution cannot be found in given range.\n",
                count_max);
        fprintf(stderr, "Exiting.\n");
        exit(0);


    }

    x_mid = (x_min + x_max)/2.0;
    f_mid = (*func)(f_mid) - nu;

    // Calculate the error
    if(nu != 0.0) err = fabs(f_mid/nu);
    else err = fabs(f_mid);

    // If err < tol, we have a solution and the calculation ends.
    if(err < tol) { return x_mid; }
    if(f_mid*f_max < 0.0) // the solution is between x_mid and x_max
    {
        // Call Solve_Bisect with the range (x_mid, x_max)
        return Solve_Bisect(nu, func, x_mid, x_max, tol, count);
    }
    else if(f_min*f_mid < 0.0) // the solution is between x_min and x_mid
    {
        return Solve_Bisect(nu, func, x_mid, x_mid, tol, count);
    }
    else // one of the factors is zero, return said factor
    {
        if(f_mid == 0.0) return x_mid;
        else if(f_max == 0.0) return x_max;
        else return x_min;
    }
}



// First, let's find the derivative at a point
// This uses f’(x) = (f(x+h) - f(x-h))/(2h) + O(h^2)
double Solve_Get_Df(double (*func)(double), double x_old)
{
    double h, df;
    if(x_old != 0.0) { h = x_old*1.0E-5; }
    else {h = 1.0E-5; }
    // This is how one accesses the function passed to this function via the function pointer (*func)
    df = (*func)(x_old+h) - (*func)(x_old-h);
    df /= 2.0*h;
    return df;
}



// solves nu = func(x) by newton’s method
// using x_{n+1} = x_n + (nu -f(x_n))/f’(x_n)
double Solve_Newton(double nu, double (*func)(double), double x_0, double tol, int *count)
{
    double x_old, x_new, err, df, h;
    int count_max;
    count_max = 1000;
    x_old = x_0; // Initial value
    do { // Do the following while the condition in the while(...)
// below is satisfied
        df = Solve_Get_Df(func, x_old); // Get the derivative
        if(fabs(df) < tol) // Derivative is too small
        {
            //Warn and exit
            fprintf(stderr, "Derivative is too small.\n",
                    count_max);
            fprintf(stderr, "Exiting.\n");
            exit(0);
        }
// TD: Implement x_new = x_old + (nu - f(x_old))/f’(x_old);
// TD: Calculate err = |((x_new-x_old)/x_old)|;
// Absolute value function is fabs(x)
// TD: Set x_old = x_new
// TD: Add 1 to *count
// This is how you access a variable that is passed by reference.
        if(*count == count_max) // Too many iterations
        {
            //Warn and exit
            fprintf(stderr, "Too many iterations.\n",
                    count_max);
            fprintf(stderr, "Exiting.\n");
            exit(0);

        }
    } while(err > tol);
    return x_new;
}





