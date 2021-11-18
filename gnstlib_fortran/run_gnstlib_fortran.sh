#!/bin/sh

CXXDIR="../gnstlib"
SOURCE_DIR="$CXXDIR/src"
INCLUDE_DIR="$CXXDIR/include"

# compile gnstlib c++
echo "Compiling GNSTLIB C++ and C API..."
g++ -c -O3 -fopenmp -lquadmath -W -fpic -std=c++11 -fext-numeric-literals \
gnstlib_capi.cpp $SOURCE_DIR/*.cpp -I $INCLUDE_DIR

# compile fortran module and C api
echo "Compiling Fortran module..."
gfortran -Wall -Wextra -c gnstlib_mod.f90 -o gnstlib_mod.o

# compile test
gfortran -Wall -Wextra -c test.f90 -o test.o

# build
echo "GNSTLIB Fortran test"
gfortran -fopenmp *.o -lquadmath -o gnstlib_fortran_test -lstdc++

# remove objects
rm *.o *.mod

echo "Done!"