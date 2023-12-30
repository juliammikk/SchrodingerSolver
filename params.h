#ifndef PARAMS_H
#define PARAMS_H

#define hbarc (197.3)

//collection of constants needed in calculation
typedef struct params{
    double mass;
    double Ea; //energy scale
    double ka; //momentum scale
    double r0; //length scale
    double x0; //r0*k
    int ell;

    char *mass_unit;
    char *length_unit;

    double Et_min; //min energy eigenvalue
    double Et_max;

    double nucA;
    double nucZ;

} Params;

//collection of variables used in calculation but that do not remain constant
typedef struct dynamic_vars{
    double Eb; //absolute value of bound energy
    double kb; //sqrt(2*mass*Eb)
    double rc; //turning point radius
    double Et; //Eb/Ea
    double xc; //ka*rc
    double rf; //last point
    double xf; //ka*rf

    double Et_min;
    double Et_max;
} DynamicVars;

extern Params PARAM_DATA; //extern means public
#endif