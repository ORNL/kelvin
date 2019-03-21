#!/bin/sh

# $1 = build dir
# $2 = install dir

# Install Kelvin
cd $1
cmake ../ -DPARSERS_DIR=$2 -DMFEM_DIR=$2 -DCMAKE_INSTALL_PREFIX=$2 -DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc -G"Eclipse CDT4 - Unix Makefiles" -DCMAKE_ECLIPSE_VERSION=4.5 -DCMAKE_BUILD_TYPE=Release
make -j5
make test
make install

