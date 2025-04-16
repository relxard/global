#!/bin/sh

rm *.mod

#----------------------------
# simple single processor

  gfortran -ogblmdl38 -O2 -fconvert=big-endian -frecord-marker=4 -ffree-form gblmdl38.f

#----------------------------
# multiple processor; before running: export OMP_NUM_THREADS={number}

# gfortran -ogblmdl38 -O2 -fopenmp -fconvert=big-endian -frecord-marker=4 -ffree-form gblmdl38.f
