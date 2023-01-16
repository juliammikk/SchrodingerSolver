#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "solve.h"

double f_solve(double x);
int main(void)
{
    double x_max, x_min, x, tol;
    int count;
    printf("Solve_Bisect\n");
    fprintf(stdout, "Solve_Bisect\n");
    count = 0;
    tol = 1.0e-10;
    x = Solve_Bisect(0.0, f_solve, 0.1, 4.0, tol, &count);
    fprintf(stdout, "count = %d\n", count);
    fprintf(stdout, "x = %e\n", x);
    fprintf(stdout, "Solve_Newton\n");
    count = 0;
    tol = 1.0e-10;
    x = Solve_Newton(0.0, f_solve, 4.0, tol, &count);
    fprintf(stdout, "count = %d\n", count);
    fprintf(stdout, "x = %e\n", x);
    return 1;
}

double f_solve(double x)
{
    return sin(x);

}
