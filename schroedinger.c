// File schroedinger.c
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "numerov.h"
#include "numerov_params.h"
#include "radial_eq_functions.h"
#include "params.h"
#include "schroedinger.h"
#include "solve.h"
#include "init.h"


NumerovParams *COM_N_PARAMS_F;
NumerovParams *COM_N_PARAMS_B;
DynamicVars *COM_D_PARAMS;
double *COM_Y_F; // Forward wavefunction
double *COM_Y_B; // Backward wavefunction

// These functions are needed only within this file
double Schroedinger_GetDf_nmax(double *y, NumerovParams *N_Params);
void Schroedinger_InitializeCommonVar(void);
double Schroedinger_GetError(void);
void Schroedinger_EvolveForward(void);
void Schroedinger_EvolveBackward(void);
void Schroedinger_CalcRunScales(double Et);
double Schroedinger_CalcRc(double Eb, double r_init);
double Schroedinger_GetBoundStateError(double Et);
void Schroedinger_PlotData(double Et_min, double Et_max);

void Schroedinger_GetBoundState
(DynamicVars *D_Params, NumerovParams *N_Params_f,
NumerovParams *N_Params_b, double *yf, double *yb)
{
    double Et_min, Et_max, tol, err;
    int count;
    double x, y;

    COM_N_PARAMS_F = N_Params_f;
    COM_N_PARAMS_B = N_Params_b;
    COM_D_PARAMS = D_Params;
    COM_Y_F = yf;
    COM_Y_B = yb;

    Et_min = D_Params->Et_min;
    Et_max = D_Params->Et_max;
    printf("et_min is: %f\n", Et_min);
    printf("et_max is: %f\n", Et_max);
    // To plot the error err = u’_I/u_I - u’_II/u_II at x = xc
    //
    Schroedinger_PlotData(Et_min, Et_max);

    // Schroedinger_GetBoundStateError returns
    // err = u’_I/u_I - u’_II/u_II at x = xc
    count = 0;
    tol = 1.0e-6;
    Solve_Bisect(0.0, Schroedinger_GetBoundStateError, Et_min, Et_max, tol, &count);
    fprintf(stderr, "count = %d\n", count);

    return;
}// Schroedinger_GetBoundState

// To plot
// err = u’_I/u_I - u’_II/u_II at x = xc
// between Et_min and Et_max
void Schroedinger_PlotData(double Et_min, double Et_max)
{
    FILE *output;
    int n, nmax;
    double dEt, Et, err;
    // Open a file to record this data and put it to output
    // Choose the file name according to the naming convention
    output = fopen("bound_state_err.dat", "w");
    // Set nmax = 1000;
    nmax = 1000;
    // Set dEt = (Et_max - Et_min)/(nmax)
    dEt = (Et_max - Et_min)/(nmax);
    for(n=0; n<=nmax; n++)
    {
        // Set Et = Et_min + n*dEt
        Et = Et_min + n*dEt;
        err = Schroedinger_GetBoundStateError(Et);
        // Print Et and err to output
        fprintf(output, "%e %e\n", Et, err);
    }
    // Don’t forget to close the file
    fclose(output);
    return;
}// Schroedinger_Plot_Data


// This function returns
// err = u’_I/u_I - u’_II/u_II at x = xc
double Schroedinger_GetBoundStateError(double Et)
{
    double err;
    Schroedinger_CalcRunScales(Et);
    Schroedinger_InitializeCommonVar();
    Schroedinger_EvolveForward();
    Schroedinger_EvolveBackward();
    err = Schroedinger_GetError();
    return err;
}// Schroedinger_GetBoundStateError


// Get the classical turning point r_c
// and use it to set xc and xf
//
void Schroedinger_CalcRunScales(double Et)
{
    double r_init, r_min, r_c, Ea, Eb;
    Ea = PARAM_DATA.Ea; // Energy scale
    Eb = Et*Ea; // Energy corresponding to Et
    // Get COM_D_PARAMS->kb = sqrt(2.0*Eb*mass)
    // getting mass from PARAM_DATA
    COM_D_PARAMS->kb = sqrt(2.0*Eb*(PARAM_DATA.mass));
    // Set Et in COM_D_Params to be Et
    COM_D_PARAMS->Et = Et;
    // Set Eb in COM_D_Params to be Eb
    COM_D_PARAMS->Eb = Eb;
    // Set r_min to be r0 from PARAM_DATA
    r_min = PARAM_DATA.r0;
    // Set r_init to be 1.1*r_min
    r_init = 1.1*r_min;
    // This calculates the classical turning point
    // by solving -Eb = V(r) + ell(ell+1)/(2 mass r^2)
    r_c = Schroedinger_CalcRc(Eb, r_init);
    // Set rc in COM_D_PARAMS to be r_c
    COM_D_PARAMS->rc = r_c;
    // Set xc in COM_D_PARAMS to be r_c*ka
    // getting ka from PARAM_DATA
    COM_D_PARAMS->xc = r_c*(PARAM_DATA.ka);
    // The wavefunction behaves like exp(-(kb/ka)*x)
    // where kb = sqrt(2*mass*Eb)
    // We take x_f = 20*(ka/kb)
    // Set COM_D_PARAMS->xf to be 20 times ka/kb.
    // getting ka and kb from appropriate parameter structures
    COM_D_PARAMS->xf = 20*(PARAM_DATA.ka)/(COM_D_PARAMS->kb);
    // Set COM_D_PARAMS->rf to be xf/ka
    // getting ka and xf from appropriate parameter structures
    COM_D_PARAMS->rf = (COM_D_PARAMS->xf)/(PARAM_DATA.ka);
    return;
}// Schroedinger_CalcRunScales


// Solve V(r) + ell(ell+1)/(2mu r^2) = -Eb
// This is given
double Schroedinger_CalcRc(double Eb, double r_init)
{
    double f, E_min, r_c, tol;
    int count;
    // solve -Eb = V_eff(r)
    tol = 1.0e-8;
    count = 0;
    r_c = Solve_Newton(-Eb, RadialEqFunctions_Veff, r_init, tol, &count);
    return r_c;
}// Get the turning point larger than Rmin



// Initialize all other necessary variables
void Schroedinger_InitializeCommonVar(void)
{
    double kb, rf, h;
    // For the forward evolution, x_i = 0 and x_f = xc
    // y_0 = y[0] = 0 and y_1 = y[1] can be any small number
    // We set it to 0.1
    // Set COM_N_PARAMS_F->x_i to be 0.0
    COM_N_PARAMS_F->x_i = 0.0;
    // Set COM_N_PARAMS_F->x_f to be COM_D_PARAMS->xc
    COM_N_PARAMS_F->x_f = COM_D_PARAMS->xc;
    // Set COM_N_PARAMS_F->y_0 to be 0.0
    COM_N_PARAMS_F->y_0 = 0.0;
    // Set COM_N_PARAMS_F->y_1 to be 0.1
    COM_N_PARAMS_F->y_1 = 0.1;
    // Set COM_N_PARAMS_F->h to be (x_f - x_i)/(COM_N_PARAMS_F->nmax);
    COM_N_PARAMS_F->h = (COM_N_PARAMS_F->x_f - COM_N_PARAMS_F->x_i)/(COM_N_PARAMS_F->nmax);

    // Backward evolution params
    // The wavefunction behaves like exp(-(kb/ka)*x)
    // where kb = sqrt(2*mass*Eb)
    // We take x_f = 20*(ka/kb)
    // so that y[0] = exp(-20) = 2E-9
    // and y[1] = exp(-(kb/ka)*(x_f-h))
    // y_0 = y[0] should be a small number and
    // y_1 = y[1] should be a small number > y[0]
    // The x range is xc < x < xf
    // or 0 < x’ < xf-xc
    // Set COM_N_PARAMS_B->x_i to be 0.0;
    COM_N_PARAMS_B->x_i = 0.0;
    // Set COM_N_PARAMS_B->x_f to be COM_D_PARAMS->xf - COM_D_PARAMS->xc;
    COM_N_PARAMS_B->x_f = COM_D_PARAMS->xf - COM_D_PARAMS->xc;
    // Set COM_N_PARAMS_B->h to be (x_f - x_i)/(COM_N_PARAMS_B->nmax);
    COM_N_PARAMS_B->h = (COM_N_PARAMS_B->x_f - COM_N_PARAMS_B->x_i)/(COM_N_PARAMS_B->nmax);

    kb = COM_D_PARAMS->kb;
    rf = COM_D_PARAMS->rf;
    h = COM_N_PARAMS_B->h;
    // Set COM_N_PARAMS_B->y_0 to be exp(-kb*rf);
    COM_N_PARAMS_B->y_0 = exp(-kb*rf);
    // SEt COM_N_PARAMS_B->y_1 to be exp(-kb*(rf-h));
    COM_N_PARAMS_B->y_1 = exp(-kb*(rf-h));
    return;
}// Schroedinger_Initialize


// This is given
void Schroedinger_EvolveForward(void)
{
    double f, df;
    int nmax;
    double *yf;
    NumerovParams *N_Params_f;
    DynamicVars *D_Params_f;
    yf = COM_Y_F;
    N_Params_f = COM_N_PARAMS_F;
    D_Params_f = COM_D_PARAMS;
    nmax = N_Params_f->nmax;
    N_Params_f->NumerovFunc_F = RadialEqFunctions_F_Forward;
    Numerov_Advance(yf, N_Params_f, D_Params_f);
    return;
}// Schroedinger_EvolveForward



// This is given
void Schroedinger_EvolveBackward(void)
{
    double f, df;
    int nmax;
    double *yb;
    NumerovParams *N_Params_b;
    DynamicVars *D_Params_b;
    yb = COM_Y_B;
    N_Params_b = COM_N_PARAMS_B;
    D_Params_b = COM_D_PARAMS;
    nmax = N_Params_b->nmax;
    N_Params_b->NumerovFunc_F = RadialEqFunctions_F_Backward;
    Numerov_Advance(yb, N_Params_b, D_Params_b);
    return;
}// Schroedinger_EvolveBackward


// This implements the numerical derivative at the end of the list.
// That is, given y[n-2], y[n-1], y[n], get the slope at y[n].
// Given these 3 points, the slope at x_n is given by
// y’[n] = (3y[n] - 4y[n-1] + y[n-2])/(2h)
//
double Schroedinger_GetDf_nmax(double *y, NumerovParams *N_Params)
{
    double df, h;
    int nmax;
    nmax = N_Params->nmax;
    h = N_Params->h;
    // TD: Implement
    // df = (3y[n] - 4y[n-1] + y[n-2])/(2h)
    // with n = nmax
    df = (3.0*y[nmax] - 4.0*y[nmax-1] + y[nmax-2])/(2.0*h);
    return df;
}// Schroedinger_GetDf_nmax



double Schroedinger_GetError(void)
{
    double *yf, *yb;
    NumerovParams *N_Params_f, *N_Params_b;
    double df, df_f, df_b;
    N_Params_f = COM_N_PARAMS_F;
    N_Params_b = COM_N_PARAMS_B;
    yf = COM_Y_F;
    yb = COM_Y_B;
    // Get the derivative at nmax using yf and N_Params_f
    df = Schroedinger_GetDf_nmax(yf, N_Params_f);
    // Devide the derivative by yf[nmax] and put it in df_f
    df_f = df/yf[N_Params_f->nmax];
    // Get the derivative at nmax using yb and N_Params_b
    df = -Schroedinger_GetDf_nmax(yb, N_Params_b);
    // Devide the derivative by yb[nmax] and put it in df_b
    df_b = df/yb[N_Params_b->nmax];
    df = df_f - df_b;
    return df;
}// Schroedinger_GetError
