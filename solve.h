//
// Created by megan on 14/01/2023.
//
#ifndef SCHRODINGERSOLVER_SOLVE_H
#define SCHRODINGERSOLVER_SOLVE_H


// Solve_Bisect solves f(x) = nu using the bisect search method
double Solve_Bisect(double nu, double (*func)(double), double x_min, double x_max, double tol, int *count);

// Solve_Newton solves f(x) = nu using Newton’s method
double Solve_Newton(double nu, double (*func)(double), double x_0, double tol, int *count);

#endif