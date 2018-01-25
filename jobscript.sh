#!/bin/sh
#mpirun -np <number_of_processes>  <executable file's absolute path>
mpirun -np 1 ~/a.out
mpirun -np 2 ~/a.out
mpirun -np 4 ~/a.out
mpirun -np 8 ~/a.out
mpirun -np 16 ~/a.out
mpirun -np 32 ~/a.out
mpirun -np 64 ~/a.out

