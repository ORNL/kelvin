#!/bin/sh

# Install Kelvin
echo $PWD
mkdir Kelvin/build
cd Kelvin/build/
cmake ../ -DPARSERS_DIR=/usr/local -DMFEM_DIR=/usr/local -DCMAKE_BUILD_TYPE=Release
make -j4
make test
# Exit with the status of the test command
exit $?

