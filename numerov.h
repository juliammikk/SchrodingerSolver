#ifndef NUMEROV_H
#define NUMEROV_H
#include "numerov_params.h"

void Numerov_Advance(double *y, NumerovParams* Num_Params, DynamicVars* Dyn_Vars);

#endif
