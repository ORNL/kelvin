HOWTO Work with Furnace
=

Prerequisites
==

Furnace requires ORNL's Fire utility library which you can download as source from https://github.com/jayjaybillings/fire. Furnace also requires the Modified Finite Element Method (MFEM) library, which you can download and build from http://mfem.org.

Fire
===

Fire can be built with all standard options. The "make install" step in Fire's build should be executed with the option -DCMAKE_INSTALL_PREFIX=<your_install_path> configured. See the Fire build instructions for more details ("Installation Step" in Fire's README.md).

Build
==

Modify the Makefile to point to the proper MFEM and Fire directories. Then simply type "make" at the prompt:

```bash
$ make
```

Input
==

Furnace uses a simple INI file format for input which consists of blocks, keys and values. See input.ini for a working example of the format.

Boundary Conditions
===

Dirichlet/Essential Boundary Conditions are specified using the conventions of MFEM. Meshes should number all boundary elements sequentially by side such that all boundary elements on side 1 are on number 1, 2 for those on side 2, etc.

Meshes
===

Mesh option block

```
[mesh]
file= # The name of the file in the local file system
order= # The order of the finite elements in the file 
```

Furnace supports any meshes supported by MFEM. Note that when a NETGEN neutral format mesh is used it is necessary to add the word "NETGEN" as the first line in the mesh file if it is not already available.