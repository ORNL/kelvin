#!/bin/sh

# Note that GCC is not installed by this script because RHEL uses an old version that is not supported by Kelvin.

# Install dependencies
echo $PWD
yum install cmake git boost boost-test doxygen make
cd $kelvin_baseDir
