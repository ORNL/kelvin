#!/bin/sh

# Install Kelvin
mkdir Kelvin/build
cd Kelvin/build/
cmake ../ -DPARSERS_DIR=/usr/local -DMFEM_DIR=/usr/local
make -j4
make test

