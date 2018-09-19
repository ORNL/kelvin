#!/bin/sh

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

