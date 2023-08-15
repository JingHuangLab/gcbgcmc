/**
* This file is part of **gcbgcmc**.
* 
* gcbgcmc is free software: you can redistribute it and/or modify it under
* the terms of the GNU Lesser General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
* 
* gcbgcmc is distributed with the hope that it will be useful,
* but **WITHOUT ANY WARRANTY**; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
* 
* See the GNU Lesser General Public License for more details.
*/

#ifndef MGGCBGCMCMOVE_H
#define MGGCBGCMCMOVE_H

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include"const.hpp"
#include"myrand.hpp"
using namespace std;


//
void translation(double xmin, double ymin, double zmin, 
                double xlen, double ylen, double zlen,
                int natoms, double (*new_position)[3], MTRand &mtrand);

//
void rotate_mol(int natoms, double (*position)[3], const double (*sposition)[3],
                const double rotation_point[3],MTRand &mtrand);

//void getCOM(int natoms, double (&com)[3], double *mass, double (*position)[3]);

void getCOG(int natoms, double (&gom)[3], double (*position)[3]);

//double rotate();


#endif