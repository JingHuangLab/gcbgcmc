# gcbgcmc
A C++ implementation of grid cavity bias GCMC interface with CHARMM. Note the CHARMM itself has a built-in GCMC module which is much faster and more robust. This implementation aims to carry out GCMC with the Drude force field so that it can be more flexible and take advantage of CHARMM's versatility. To use this program, by default, it is thought you are familiar with CHARMM.

# Compilation
**Prerequisites**:   
- **MersenneTwister.h**: a random number generator that can be found on GitHub.  
- a C++ compiler support **C++11** standard.  
- cmake ( version >= **3.17** works fine, lower version may also work.)

To compile the code:
```bash
mkdir build
cd build
cmake ../src -DCMAKE_BUILD_TYPE=Release
cmake --build .
```
The executable binary file **mygcbgcmc.exe** will be generated.
# Usage
To use mygcbgcmc, the molecular dynamics (MD) engine CHARMM need to be properly installed, so that mygcbgcmc can call CHARMM to calculate energy, force and propagate MD simulations.
To run the GCMC tasks, you need to prepare X files:
 - input file for mygcbgcmc (control the GCMC job, e.g. steps, GCMC region...)
 - input template for CHARMM (control how force and energy are calculated and how MD are carried out)
 - CHARMM coordinate file (define the system coordinate)
 - CHARMM psf file (define the system topology)
Then run GCMC by:
```bash
mygcbgcmc.exe -i input -o log
```
Details about the input file are provided in the example/ and the corresponding README.md.  

# Reference

 1.  [D.J. Adams (1975): Grand canonical ensemble Monte Carlo for a Lennard-Jones fluid, Molecular Physics: An International Journal at the Interface Between Chemistry and Physics, 29:1, 307-311](http://dx.doi.org/10.1080/00268977500100221)
 2.  [Mihaly Mezei (1980) A cavity-biased (T, V, µ) Monte Carlo method for the computer simulation of fluids, Molecular Physics: An International Journal at the Interface Between Chemistry and Physics, 40:4, 901-906](http://dx.doi.org/10.1080/00268978000101971)
 3.  [Mihaly Mezei (1987) Grand-canonical ensemble Monte Carlo study of dense liquid, Molecular Physics: An International Journal at the Interface Between Chemistry and Physics, 61:3, 565-582](http://dx.doi.org/10.1080/00268978700101321)

# Contact

Issues are welcomed for discussion and suggestions.
