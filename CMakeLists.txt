cmake_minimum_required(VERSION 3.20)
project(SchrodingerSolver C)

set(CMAKE_C_STANDARD 99)

add_executable(SchrodingerSolver main.c solve.h solve.c)
