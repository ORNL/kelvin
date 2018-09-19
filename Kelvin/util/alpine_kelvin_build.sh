#!/bin/sh

# Install dependencies
apk add gcc g++ cmake git boost boost-dev doxygen make

# Install parsers
cd home
git clone https://github.com/jayjaybillings/parsers
mkdir parsers/build
cd parsers/build/
cmake ..
make
make test
make install
make clean

# Install mfem
cd /home
git clone https://github.com/mfem/mfem
mkdir mfem/build
cd mfem/build/
cmake ..
make
make install
make clean

# Install Kelvin
mkdir /home/furnace-prototype/Kelvin/build
cd /home/furnace-prototype/Kelvin/build/
cmake ../ -DPARSERS_DIR=/usr/local -DMFEM_DIR=/usr/local
make -j4
make test

