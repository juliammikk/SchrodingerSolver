cmake_minimum_required(VERSION 3.20)
project(SchrodingerSolver C)

set(CMAKE_C_STANDARD 99)

add_executable(SchrodingerSolver solve.c derivatives.h derivatives.c extremum.c extremum.h numerov_params.h params.h vector_mtx.c vector_mtx.h numerov.c numerov.h main.c init.c init.h radial_eq_functions.c radial_eq_functions.h schroedinger.c schroedinger.h)
target_link_libraries(SchrodingerSolver m)
