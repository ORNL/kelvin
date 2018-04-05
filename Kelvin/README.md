HOWTO Work with Furnace
=

Prerequisites
==

Kelvin requires ORNL's Fire utility library which you can download as source from https://github.com/jayjaybillings/fire. Kelvin also requires the Modified Finite Element Method (MFEM) library, which you can download and build from http://mfem.org.

Fire
===

Fire can be built with all standard options. The "make install" step in Fire's build should be executed with the option -DCMAKE_INSTALL_PREFIX=<your_install_path> configured. See the Fire build instructions for more details ("Installation Step" in Fire's README.md).

Build
==

Kelvin uses CMake to compiles and build the source code into a library and executables. The entire project can be configured and build with a few easy steps, starting from the source directory:

```bash
$ mkdir build/
$ cd build/
$ cmake ../ -DPARSERS_DIR=$PARSERS_ROOT_DIR -DMFEM_DIR=$MFEM_ROOT_DIR
$ make
```

Where $PARSERS_ROOT_DIR and $MFEM_ROOT_DIR are the locations of the Parsers and MFEM installation directories. Both directories are the installation directories, not the build or source directories. The following is an example for a user name bob:

```bash
$ cmake ../ -DPARSERS_DIR=/home/bob/parsers-install -DMFEM_DIR=/home/bob/mfem-install
```

For developers who are working with Eclipse, CMake will auto-generate the Eclipse project files with the following configuration options and setup a debug build with the following configuration:

```bash
$ cmake ../ -DCMAKE_BUILD_TYPE=Debug -G"Eclipse CDT4 - Unix Makefiles" -DCMAKE_ECLIPSE_VERSION=4.5 -DPARSERS_DIR=/home/bob/parsers-install -DMFEM_DIR=/home/bob/mfem-install 
```

Input
==

Kelvin uses a simple INI file format for input which consists of blocks, keys and values. See input.ini for a working example of the format.

Boundary Conditions
===

Dirichlet/Essential Boundary Conditions are specified using the conventions of MFEM. Meshes should number all boundary elements sequentially by side such that all boundary elements on side 1 are on number 1, 2 for those on side 2, etc.

By default the thermal boundary condition on each side is set using the surfaceTemperature element of the [thermal] block. Heat fluxes can be specified on a side-by-side basis. Initial temperature values for all interior nodes is set using the initialTemperature element.

Meshes
===

Mesh option block

```
[mesh]
file= # The name of the file in the local file system
order= # The order of the finite elements in the file 
```

Kelvin supports any meshes supported by MFEM. Note that when a NETGEN neutral format mesh is used it is necessary to add the word "NETGEN" as the first line in the mesh file if it is not already available.