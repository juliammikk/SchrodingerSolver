#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "radial_eq_functions.h"
#include "numerov_params.h"
#include "params.h"
#define alpha_EM (1.0/137.0) // The fine structure constant


// This is in the unit of 1/fm or 1/nm
double RadialEqFunctions_V(double r)
{
    double f;
    double A, R_A, a, R0, V0;

    if(PARAM_DATA.nucA == 0.0){
        f = -alpha_EM/r;
    } else {
        // Nuclear Woods-Saxon potential
        V0 = 50.0/hbarc; // MeV to 1/fm
        a = 0.7;
        R_A = PARAM_DATA.r0;
        f = -V0/(1.0 + exp((r-R_A)/a));
    }
    return f;
}

//*** No user serviceable parts from here on ***//


// Veff(r) = V(r) + ell(ell+1)/(2 m r^2)
double RadialEqFunctions_Veff(double r)
{
    double f, ell, mass;
    ell = (double) PARAM_DATA.ell; // PARAM_DATA.ell is an integer

    mass = PARAM_DATA.mass;

    f = RadialEqFunctions_V(r) + ell*(ell+1)/(2*mass*r*r);

    return f;
}


// F in y’’ = Fy for u_I
// This is in x = ka*r
double RadialEqFunctions_F_Forward(double x, DynamicVars *Dyn_Vars)
{
    double x0, ka, r, f, g, Ea, Et, ell, eps;

    ell = (double) PARAM_DATA.ell;
    x0 = PARAM_DATA.x0;
    ka = PARAM_DATA.ka;
    Ea = PARAM_DATA.Ea;
    Et = Dyn_Vars->Et;

    // Small number to prevent x = 0
    eps = 1.0e-15;
    x += eps;
    r = x/ka;

    f = ell*(ell+1)/pow(x,2) + RadialEqFunctions_V(r)/Ea + Et;

    return f;
}// RadialEqFunctions_F_Forward



// F in Y’’ = FY Backward evolution
// Here we use y = x_f - x
double RadialEqFunctions_F_Backward(double y, DynamicVars *Dyn_Vars)
{
    double x0, ka, r, f, g, ell, Ea, Et, x;

    ell = (double) PARAM_DATA.ell;
    x0 = PARAM_DATA.x0;
    ka = PARAM_DATA.ka;
    Ea = PARAM_DATA.Ea;
    Et = Dyn_Vars->Et;

    x = Dyn_Vars->xf - y;
    r = x/ka;

    f = ell*(ell+1)/pow(x,2) + RadialEqFunctions_V(r)/Ea + Et;
    return f;
}// RadialEqFunctions_F_Backward
