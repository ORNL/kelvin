PMGen - Particle Mesh Generator
=

This is a simple utility script for generating the particle list and background mesh for a given input mesh.

PMGen uses other classes in Kelvin and MFEM to generate the particle mesh.

Usage
==

Run `pmgen -h` for a full description of the options available in PMGen. 

Input 
==

Any mesh supported by MFEM will work.

Output
==

PMGen will generate a comma separated list of particle positions in a comma separate variables (CSV) file called particles.csv by default. 

