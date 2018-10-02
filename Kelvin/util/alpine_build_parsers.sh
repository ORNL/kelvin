#!/bin/sh

# Install parsers
echo $PWD
ls
cd /home
echo $PWD
ls
git clone https://github.com/jayjaybillings/parsers
mkdir parsers/build
cd parsers/build/
cmake ..
make
make test
make install
make clean
cd $kelvin_basedir

