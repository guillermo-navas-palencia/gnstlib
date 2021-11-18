#!/bin/sh

# Generate Python bindings for GNSTLIB C++ library

# switch to anaconda 3.5 (optional)
export PATH=/home/$USER/anaconda3/bin:$PATH

# project parameters
PROJECT="gnstlib"
PROJECT_ROOT="../$PROJECT"
ANACONDA="/home/$USER/anaconda3/include/python3.*m/"
NUMPY="/home/$USER/anaconda3/lib/python3.*/site-packages/numpy/core/include/"

SOURCE_DIR="$PROJECT_ROOT/src"
INCLUDE_DIR="$PROJECT_ROOT/include"

# swig
echo "generating swig wrapper..."
swig -c++ -python gnstlib.i

#compile
echo "compiling gnstlib library..."
g++ -c -O3  -fopenmp -lquadmath -W -fpic -std=gnu++11 \
    $SOURCE_DIR/*.cpp gnstlib_wrap.cxx -I $ANACONDA -I $INCLUDE_DIR -I $NUMPY

#build library
echo "building python module..."
g++ -shared  *.o -lgomp -lquadmath -o  _gnstlib.so

#move objects to build
mkdir -p build
mv *.o gnstlib_wrap.cxx build/