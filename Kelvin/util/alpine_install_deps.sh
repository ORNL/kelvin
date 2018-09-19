#!/bin/sh

# Install dependencies
echo $PWD
apk add gcc g++ cmake git boost boost-dev doxygen make
cd $kelvin_basedir
