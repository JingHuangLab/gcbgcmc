/**
* This file is part of **mygcbgcmc**.
* 
* mygcbgcmc is free software: you can redistribute it and/or modify it under
* the terms of the GNU Lesser General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
* 
* mygcbgcmc is distributed with the hope that it will be useful,
* but **WITHOUT ANY WARRANTY**; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
* 
* See the GNU Lesser General Public License for more details.
*/

#ifndef MYGCBGCMCGCMCMOL_H
#define MYGCBGCMCGCMCMOL_H

#include<iostream>
#include<string>
#include"myrand.hpp"
using namespace std;

#define NATOMS 3

class GCMCMOLE {
    
    public:
    int natoms,nrealatoms;
    double mu; //chemical potential

    string resn; //resname

    // pay attention to the following comment out !
    // since the arguments are pointer so no need to do dynamical memory management
    const double (*smolecule)[3];// = new double[natoms][3];
    double *atomMassList;// = new double[natoms];
    string *atomNameList;// = new string[natoms];
    string *atomTypeList;// = new string[natoms];

    //double COM[3]={0.0, 0.0, 0.0};
    double GOM[3]={0.0, 0.0, 0.0}; // center of mass and geometry


    GCMCMOLE(int natoms, int nrealatoms, double aMu, const double (*molecule)[3],
                string aresname, string* typelist, string* namelist);

    void new_position(double xmin, double ymin, double zmin, 
                        double xlen, double ylen, double zlen,
                        int natoms, double (*position)[3],MTRand &mtrand);

};

#endif