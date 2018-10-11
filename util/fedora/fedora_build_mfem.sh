#!/bin/sh

# $1 = build dir
# $2 = install dir

# Install mfem
cd $1
git clone https://github.com/mfem/mfem
mkdir mfem/build
cd mfem/build/
cmake .. -DCMAKE_INSTALL_PREFIX=$2
make -j4
make install
make clean

