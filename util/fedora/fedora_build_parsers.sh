#!/bin/sh

# $1 = build dir
# $2 = install dir

# Install parsers
cd $1
echo $PWD
git clone https://github.com/jayjaybillings/parsers
mkdir parsers/build
cd parsers/build/
cmake .. -DCMAKE_INSTALL_PREFIX=$2
make -j3
make test
make install
make clean

