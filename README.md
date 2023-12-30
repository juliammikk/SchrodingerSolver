# Schrodinger Solver






To run either:

1) Set the program arguments in your preferred IDE to input_data then run the main method.
(I've found that to make this work you also need to move the two input files into the cmake_build_debug folder)

2) Run from terminal using
gcc -o SchrodingerSolver main.c params.h numerov_params.h numerov.c vector_mtx.c init.c radial_eq_functions.c solve.c derivatives.c extremum.c -lm

and then
.\SchrodingerSolver.exe input_coulomb input_n_params

The parameters will be printed to the file "recorded_params.dat" in the cmake-build-debug directory.



