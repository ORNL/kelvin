#!/bin/sh

# Delete unneeded dependencies
apk del gcc g++ cmake git boost boost-dev doxygen make

# Reinstall the C++ runtime library
apk add libstdc++

# Install Kelvin
cd /home/furnace-prototype/Kelvin/build
make install

# Clean up source directories
cd /home
rm -rf furnace-prototype mfem parsers
