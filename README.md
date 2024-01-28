# Schrodinger Solver

A program to numerically solve the radial part of the Schrödinger equation given an arbitrary confining potential V(r).

<img src="http://www.sciweavers.org/tex2img.php?eq=%28-%5Cfrac%7B1%7D%7B2%5Cmu%7D%5Cfrac%7Bd%5E%7B2%7D%7D%7Bdr%5E%7B2%7D%7D%2B%5Cfrac%7Bl%28l%2B1%29%7D%7B2%5Cmu%20r%5E%7B2%7D%7D%2BV%28r%29%29u_%7Bn_%7Br%7Dl%7D%28r%29%20%3D%20E_%7Bn_%7Br%7Dl%7D%20u_%7Bn_%7Br%7Dl%7D%28r%29&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0" align="center" border="0" alt="(-\frac{1}{2\mu}\frac{d^{2}}{dr^{2}}+\frac{l(l+1)}{2\mu r^{2}}+V(r))u_{n_{r}l}(r) = E_{n_{r}l} u_{n_{r}l}(r)" width="375" height="50" />


Finds both the bound state wavefunctions and the energy eigenvalues, using either the Shooting method or the Numerov method.




To run either:

1) Set the program arguments in your preferred IDE to input_data then run the main method.
(I've found that to make this work you also need to move the two input files into the cmake_build_debug folder)

2) Run from terminal using
gcc -o SchrodingerSolver main.c params.h numerov_params.h numerov.c vector_mtx.c init.c radial_eq_functions.c solve.c derivatives.c extremum.c -lm

and then
.\SchrodingerSolver.exe input_coulomb input_n_params

The parameters will be printed to the file "recorded_params.dat" in the cmake-build-debug directory.



