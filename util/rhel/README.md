These scripts will help you build Kelvin on Red Hat Enterprise Linux versions 7.6 and greater.

Prerequisites
==
You must install a modern C++ compiler in order for Kelvin to compile. You can do this by compiling from scratch or installing gcc 7 or 8 from the developer tools. Once you have installed an appropriate compiler, and checked it by executing ```gcc --version```, then you may execute the build scripts as usual by executing ```sh ./rhel_buildAll.sh```.

Notes
==
The scripts attempt to overwrite the cc and cxx links in /usr/bin  with gcc and g++, which should be pointing to modern versions if the compiler has been updated as described above. Each script passes ```-DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++``` to the CMAKE script. If a user needs to use a different compiler, they must update each script.
