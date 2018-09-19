#!/bin/sh

./Kelvin/util/alpine_before_script.sh
./Kelvin/util/alpine_install_deps.sh
./Kelvin/util/alpine_build_parsers.sh
./Kelvin/util/alpine_build_mfem.sh
./Kelvin/util/alpine_build_kelvin.sh

