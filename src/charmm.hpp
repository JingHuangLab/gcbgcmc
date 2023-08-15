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


#ifndef MYGCBGCMCCHARMM_H
#define MYGCBGCMCCHARMM_H


#include<iostream>
#include<string>
#include<vector>
#include<fstream>
#include<iomanip>
#include"gcmcmole.hpp"

using namespace std;

//CharmmSystem
// hold a class for charmm system
class CharmmSystem{
    
    public:
        
        int natoms;
        vector<int> resid;
        vector<int> resnumber;
        vector<string> resname;
        vector<string> atomtype;
        vector<string> segid;
        vector<double> xpos;
        vector<double> ypos;
        vector<double> zpos;

        void clear_all();
};


class Charmm {

    public:

    string charmmEXE;
    string inpTemplate;
    //double dens;//density
    double systemEnergy;

    CharmmSystem crd;//
    GCMCMOLE gcmcmol;//

    vector<int> atomIndex;//for the method getAtomIndex
    vector<int> residList;//for the method findResname

    vector<double> tmpx,tmpy,tmpz;
    vector<double> backupx,backupy,backupz;

    int nReal;

    Charmm(string crdfile, string psffile, GCMCMOLE agcmcmol, double density,
            string intputTemplate, string charmmexec);

    //
    inline void clearTMPxyz(){
        if (!tmpx.empty()){
            tmpx.clear();
            tmpy.clear();
            tmpz.clear();
        }
    }
    // 
    void readCrd(string crdfile);
    //
    void update();
    //
    void restore();
    //
    void getAtomIndex(int resid, string resname,int natoms);
    //
    void writeCrd();
    //
    void writeCrd(int natoms,double (*position)[3],string resname,
                        string* atomtype);
    //
    void findResname(string resname);
    //
    bool ngcmcmols(string resname);
    //
    void get_gcmc_atom_position();
    //
    void getRealAtomPosition();
    //
    void getNRealAtoms();
    //
    double insertion(int natoms, double (*insertPos)[3],
                        string resname,string* atomtype);
    //
    void restoreInsertion();
    //
    double deletion(int natoms, int resid, string resname);
    //
    void restore_deletion();
    //
    double displacement(int natoms, int iresid, double (*newPos)[3]);
    //
    void restoreDisplacement();
    //
    double runMD();
    //
    double inline readEnergy(string filename);
    //
    double energy(int jobtype=0);
    //
    void check_output();

};


#endif