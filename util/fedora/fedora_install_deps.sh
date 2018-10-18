#!/bin/sh

# Install dependencies
echo $PWD
dnf install gcc cmake git boost boost-test doxygen make
cd $kelvin_baseDir
