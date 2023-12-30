#ifndef RADIAL_EQ_FUNCTIONS_H
#define RADIAL_EQ_FUNCTIONS_H
#include "params.h" // Because we use DynamicVars below

// Potential energy in r
double RadialEqFunctions_V(double r);

// V_eff = V + ell(ell+1)/(2m r^2)
double RadialEqFunctions_Veff(double r);

// The F function in the RHS of the differential equation y’’ = Fy
// This is for the evolution forward from x = x_i
double RadialEqFunctions_F_Forward(double x, DynamicVars *Dyn_Vars);

// This is for the evolution backward from x = x_f
double RadialEqFunctions_F_Backward(double x, DynamicVars *Dyn_Vars);

#endif