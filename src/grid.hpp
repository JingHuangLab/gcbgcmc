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

#ifndef MYGCBGCMCGRID_H
#define MYGCBGCMCGRID_H

#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<cmath>
#include<omp.h>
using namespace std;

/**
 * @brief:
*/
class Grid{

    public:
        double xmin,ymin,zmin,xmax,ymax,zmax,dx0,dy0,dz0,r_cut,r_cut2;
        double dxcav,dycav,dzcav,hdxcav,hdycav,hdzcav;
        double xlen,ylen,zlen,hxlen,hylen,hzlen;
        double V;
        int nx,ny,nz,ncells;
        int xsdepth,ysdepth,zsdepth;
        int ncavities=0,tmp_ncavities=0;
        //
        vector<double> xcell_center;
        vector<double> ycell_center;
        vector<double> zcell_center;

        vector<int> occupancy,tmp_occupancy;
        vector<int> changedCells;
        vector<int> searchList;// 1D contain 3D
        //vector<int> dsa;

        vector<int> cavities;
        vector<int> tmp_cavities;
        
        vector<int> cavityList;

        //bool *cavities

        double cellBoundary[3]={0.0};


        

        //
        Grid(double x_min, double y_min, double z_min, 
            double x_max, double y_max, double z_max, 
            double dx, double dy, double dz, double r);
        //
        void initialCavityList(int natoms, vector<double> &xpos, 
                    vector<double> &ypos, vector<double> &zpos);
        //
        void getCellBoundary(int icell);
        //
        inline int cellPos2Index(int i, int j, int k);
        //
        bool ifOccupied(double &x0, double &y0, double &z0, 
                    double &x1, double &y1, double &z1,
                    bool pbcx=true, bool pbcy=true, bool pbcz=true) const;
        //
        void directionIndex(int index, int nd, int sd, vector<int> &dsa);
        //
        void updateArea(double *changedPosition);
        //
        void updateInsertion(int natoms, double (*changedPosition)[3], int n=0);
        //
        void updateDeletion(int natoms, double (*changedPosition)[3], int n=0);
        //
        void updateDisplacement(int natoms, vector<double> &xpos, 
                    vector<double> &ypos, vector<double> &zpos);
        //
        void updateDelOccupancy();
        //
        void updateInsOccupancy();
        //
        void getCavities();
        //
        void restoreDeletion();
        //
        double calProbcavities() const;

        void writeOccupancy(string fname,int n);

        //~Grid();
    
    
};

#endif
