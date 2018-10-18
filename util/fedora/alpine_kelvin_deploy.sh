#!/bin/sh

# Delete unneeded dependencies
apk del gcc g++ cmake git boost boost-dev doxygen make

# Reinstall the C++ runtime library
apk add libstdc++

# Install Kelvin
cd Kelvin/build
make install

# Clean up source directories
cd $kelvin_basedir
cd ..
rm -rf furnace-prototype 
