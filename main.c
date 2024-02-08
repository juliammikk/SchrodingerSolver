#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "params.h"
#include "numerov_params.h"
#include "init.h"
#include <errno.h>
#include "schroedinger.h"
#include "vector_mtx.h"


// Functions needed only in this file
// Read in Params parameters
void ReadIn_Params(char *input_file, DynamicVars *Dyn_Vars);


void ReadIn_Num_Params
        (char *input_file_name,
         NumerovParams *Num_Params_f, NumerovParams *Num_Params_b);


// Record parameters
void Record_Params(NumerovParams Num_Params_f, NumerovParams Num_Params_b);
// Record results
void Record_Results(DynamicVars Dyn_Vars,
                    NumerovParams Num_Params_f, NumerovParams Num_Params_b,
                    double *yf, double *yb);



Params PARAM_DATA;


// This program is to be invoked as
// executable input_file1 input_file2
// the word "input_file1" and "input_file2" are then put into "argv" below
// argv[0] is the name of the executable (e.g. schroedinger)
// argv[1] is the name of the first input file (e.g. input_coulomb)
// argv[2] is the name of the second input file (e.g. input_n_params)
int main(int argc, char **argv)
{
    DynamicVars Dyn_Vars; // These parameters are calculated
    NumerovParams Num_Params_f; // For the forward evolution of u_I
    NumerovParams Num_Params_b; // For the backward evolution of u_II

    double *yf, *yb; // yf contains u_I, yb contains u_II

    ReadIn_Params(argv[1],&Dyn_Vars); // Reads in the initial data to first input file in arguments

    // Reads in the initial data from the second input file
    ReadIn_Num_Params(argv[2], &Num_Params_f, &Num_Params_b);


    Init_CalcScales(); // Get the energy and length scales to prepare for solving the differential eq

    // Record the parameters
    Record_Params(Num_Params_f, Num_Params_b);

    // Allocate memory for the forward wavefunction yf
    yf = vector_malloc(Num_Params_f.nmax+1);
    // Allocate memory for the backward wavefunction yb
    yb = vector_malloc(Num_Params_b.nmax+1);

    Schroedinger_GetBoundState(&Dyn_Vars, &Num_Params_f, &Num_Params_b, yf, yb);

    Record_Results(Dyn_Vars, Num_Params_f, Num_Params_b, yf, yb);

    return 0;
}// main



void ReadIn_Params(char *input_file, DynamicVars *Dyn_Vars)
{
    FILE *input;
    double x;
    int ix;
    char *mass_unit;
    input = fopen(input_file, "r"); // Open the input file to "r"ead
    if (input == NULL) {
        perror("Error opening input file");
        // handle the error, such as exiting the program or returning an error code
    }


   double mass;
   fscanf(input, "%le", &mass);
   mass /= hbarc;
   PARAM_DATA.mass = mass;

// Read in the mass unit
    mass_unit = (char *) malloc(sizeof(char)*10);
// First allocate enough memory to hold it.

// From the second line, read in a line and put it in PARAM_DATA.mass_unit
    fscanf(input, "%s", mass_unit);
    PARAM_DATA.mass_unit = mass_unit;

    if(strcmp(mass_unit, "eV")==0) {PARAM_DATA.length_unit = "nm";}
    else if(strcmp(mass_unit, "MeV")==0) {PARAM_DATA.length_unit = "fm";}
    else {
        fprintf(stderr, "ReadIn_Params: %s is an unknown unit.\n", mass_unit);
        fprintf(stderr, "Known units are eV and MeV.\n");
        fprintf(stderr, "Exiting.\n");
        exit(0);
    }

    //read in orbital angular momentum
    fscanf(input, "%d", &ix);
    PARAM_DATA.ell = ix;

    //Read in the atomic mass A
    fscanf(input, "%le", &x);
    PARAM_DATA.nucA = x;

    //Read in the atomic mass Z
    fscanf(input, "%le", &x);
    PARAM_DATA.nucZ = x;

    fscanf(input, "%le", &x);
    Dyn_Vars->Et_min = x;

    fscanf(input, "%le", &x);
    Dyn_Vars->Et_max= x;

    fclose(input); // Always close an opened file
}// ReadIn_Params


void Record_Params(NumerovParams Num_Params_f, NumerovParams Num_Params_b){
    double x;
    int i;
    FILE *output;

    output = fopen("recorded_params.dat", "w");

    fprintf(output, "mass = %le %s\n", PARAM_DATA.mass*hbarc, PARAM_DATA.mass_unit);
    fprintf(output, "r0 = %le %s\n", PARAM_DATA.r0, PARAM_DATA.length_unit);
    fprintf(output, "Ea = %le %s\n", PARAM_DATA.Ea*hbarc, PARAM_DATA.mass_unit);
    fprintf(output, "ka = %le\n", PARAM_DATA.ka*hbarc);
    fprintf(output, "ell = %d\n", PARAM_DATA.ell);
    fprintf(output, "x0 = %e\n", PARAM_DATA.x0);
    fprintf(output, "nucA = %e\n", PARAM_DATA.nucA);
    fprintf(output, "nucZ = %e\n", PARAM_DATA.nucZ);
    fprintf(output, "nmax (forward) = %d\n", Num_Params_f.nmax);
    fprintf(output, "nmax (backward) = %d\n", Num_Params_b.nmax);

    fclose(output);

    return;
}// Record_Params


// Read in NumerovParams data
void ReadIn_Num_Params(char *input_file_name, NumerovParams *Num_Params_f, NumerovParams *Num_Params_b)
{
    FILE *input_file;
    double x;
    int ix;

    input_file = fopen(input_file_name, "r");

    fscanf(input_file, "%d", &ix);
    Num_Params_f->nmax = ix;

    fscanf(input_file, "%d", &ix);
    Num_Params_b->nmax = ix;

    fclose(input_file);

    return;
}// ReadIn_Num_Params


void Record_Results(DynamicVars Dyn_Vars,
                    NumerovParams Num_Params_f, NumerovParams Num_Params_b,
                    double *yf, double *yb)
{
    double Et;
    FILE *output;
    int n;
    double x;


    output = fopen("schroedinger_results.dat","w");
    Et = Dyn_Vars.Et;
    fprintf(output, "Et = %e\n", Et);
    fprintf(output, "Eb = %e %s\n", Et*PARAM_DATA.Ea*hbarc, PARAM_DATA.mass_unit);
    fclose(output);


// Record the forward going solution
    output = fopen("forward_solution.dat","w");
    for(n=0; n<=Num_Params_f.nmax; n++)
    {
        x = Num_Params_f.x_i + n*Num_Params_f.h; // Dimensionless
        x /= PARAM_DATA.ka; // 1/energy
        x *= hbarc; // length dimension
        fprintf(output, "%e %e\n", x, yf[n]/yf[Num_Params_f.nmax]);
    }
    fclose(output);
// Record the backward going solution
    output = fopen("backward_solution.dat","w");
    for(n=0; n<=Num_Params_b.nmax; n++)
    {
        x = Num_Params_b.x_f - (Num_Params_b.x_i + n*Num_Params_b.h); // Dimensionless
        //x /= PARAM_DATA.ka; // 1/energy ???
        x *= hbarc; // length dimension
        fprintf(output, "%e %e\n", x, yb[n]/yb[Num_Params_b.nmax]);
    }
    fclose(output);
}// Record_Results
