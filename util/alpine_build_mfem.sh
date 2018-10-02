#!/bin/sh

# Install mfem
echo $PWD
cd /home
git clone https://github.com/mfem/mfem
mkdir mfem/build
cd mfem/build/
cmake ..
make
make install
make clean
cd $kelvin_basedir

