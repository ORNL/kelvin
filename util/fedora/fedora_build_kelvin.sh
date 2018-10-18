#!/bin/sh

# $1 = build dir
# $2 = install dir

# Install Kelvin
cd $1
cmake ../ -DPARSERS_DIR=$2 -DMFEM_DIR=$2 -DCMAKE_INSTALL_PREFIX=$2
make -j4
make test
make install

