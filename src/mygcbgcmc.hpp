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

#ifndef CPPMYGCBGCMC_H
#define CPPMYGCBGCMC_H

#include<iostream>
#include<string>
#include<cmath>
#include"charmm.hpp"
#include"grid.hpp"
#include"const.hpp"
#include"grid.hpp"
#include"myrand.hpp"


using namespace std;

class GCBGCMC {
    public:
        MTRand mtrand;
        Charmm mysystem;
        Grid mygrid;
        int nstep=0;
        int ndisplacement=0;
        int displacement_accept=0;
        int ndeletion=0;
        int deletion_accept=0;
        int ninsertion=0;
        int insertion_accept=0;

        double beta;
        double Bfactor;
        
        double old_energy,new_energy;

        //template pointer for the insertion and deletion position
        double (*insertPos)[3];
        double (*deletePos)[3];


        // constructor
        GCBGCMC(Charmm &system, Grid &mygrid, double temperature);

        // 
        bool acceptance(double dE, double cbfactor);
        
        //
        void statastic() const;
        
        //displacement
        void displacement();
        
        //deletion
        void deletion();
        
        //insertion
        void insertion();
        
        //GCMC main loop
        void run(int nsteps);

        //GCMC/MD main loop
        void run(int nsteps, int nmc);


        //~GCBGCMC();
        
        
};

// inherit from GCBGCMC to do special care for Drude
class DGCBGCMC: GCBGCMC{
    public:

        // use the same constructor as the GCBGCMC
        using GCBGCMC::GCBGCMC;

        //deletion
        void ddeletion();
        
        //insertion
        void dinsertion();

        //GCMC main loop
        void run(int nsteps);

        //GCMC/MD main loop
        void run(int nsteps, int nmc);


};

#endif