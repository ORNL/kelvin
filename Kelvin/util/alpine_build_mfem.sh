#!/bin/sh

# Install mfem
cd /home
git clone https://github.com/mfem/mfem
mkdir mfem/build
cd mfem/build/
cmake ..
make
make install
make clean

