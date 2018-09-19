#!/bin/sh

# Install Kelvin
echo $PWD
mkdir Kelvin/build
cd Kelvin/build/
cmake ../ -DPARSERS_DIR=/usr/local -DMFEM_DIR=/usr/local
make -j4
make test
cd $kelvin_basedir

