#!/bin/sh

./alpine_before_script.sh
./alpine_install_deps.sh
./alpine_build_parsers.sh
./alpine_build_mfem.sh
./alpine_build_kelvin.sh

