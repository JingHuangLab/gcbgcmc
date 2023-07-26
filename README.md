# mygcbgcmc
A C++ implemention of grid cavity bias GCMC interface with CHARMM.

# Compilation
**Prerequests**:   
- **MersenneTwister.h**: a random number generator.  
- a C++ compiler support **C++11** standard  
- cmake version >= **3.17** (works fine, lower version may also work)

To compile the code:
```bash
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build .
```
# Usage

# Reference

 1.  [D.J. Adams (1975): Grand canonical ensemble Monte Carlo for a Lennard-Jones fluid, Molecular Physics: An International Journal at the Interface Between Chemistry and Physics, 29:1, 307-311](http://dx.doi.org/10.1080/00268977500100221)
 2.  [Mihaly Mezei (1980) A cavity-biased (T, V, µ) Monte Carlo method for the computer simulation of fluids, Molecular Physics: An International Journal at the Interface Between Chemistry and Physics, 40:4, 901-906](http://dx.doi.org/10.1080/00268978000101971)
 3.  [Mihaly Mezei (1987) Grand-canonical ensemble Monte Carlo study of dense liquid, Molecular Physics: An International Journal at the Interface Between Chemistry and Physics, 61:3, 565-582](http://dx.doi.org/10.1080/00268978700101321)

# Contact
