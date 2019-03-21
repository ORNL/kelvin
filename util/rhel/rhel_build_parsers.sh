#!/bin/sh

# $1 = build dir
# $2 = install dir

# Install parsers
cd $1
echo $PWD
echo $1 $2
git clone https://github.com/jayjaybillings/parsers
mkdir parsers/build
cd parsers/build/
cmake .. -DCMAKE_INSTALL_PREFIX=$2 -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++
make -j5 VERBOSE=1
make test
make install
make clean

