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

#include"gcmcmove.hpp"



/**
 * @brief Translation
 * @param xmin the mini along x axis
 * @param ymin the mini along y axis
 * @param zmin the mini along z axis
 * @param xlen the box length along x axis
 * @param ylen the box length along y axis
 * @param zlen the box length along z axis
 * @param mtrand a MTRand class used for random number
 * 
*/
void translation(double xmin, double ymin, double zmin, 
                double xlen, double ylen, double zlen,
                int natoms, double (*new_position)[3], MTRand &mtrand){
    //random point first
    double randomx=xmin+myrand(mtrand)*xlen;
    double randomy=ymin+myrand(mtrand)*ylen;
    double randomz=zmin+myrand(mtrand)*zlen;
    
    for(int i=0;i!=natoms;i++){ // the nth atoms
        new_position[i][0]+=randomx;
        new_position[i][1]+=randomy;
        new_position[i][2]+=randomz;
    }
}

/**
 * @brief Randomly rotate a molecule according to a point
 * 
 * @param natoms number of atoms in the molecule
 * @param position position of the return value
 * @param sposition standary position
 * @param rotationPoint rotation according to which
 * @param mtrand a MTRand class to generate random number 
*/
void rotate_mol(int natoms, double (*position)[3], 
                const double (*sposition)[3],
                const double rotation_point[3], MTRand &mtrand){
    
    //prepare a rotation vecotr
    double v[3]={0.0, 0.0, 0.0};
    while(true){
        v[0]=2*myrand(mtrand)-1;
        v[1]=2*myrand(mtrand)-1;
        v[2]=v[0]*v[0]+v[1]*v[1];
        if(v[2]<1){break;}
    }

    v[0] = 2*v[0]*sqrt(1 - v[2]);
    v[1] = 2*v[1]*sqrt(1 - v[2]);
    v[2] = 1 - 2*v[2];

    double u[3][3]={{1.0, 0.0, 0.0},
                    {0.0, 1.0, 0.0},
                    {0.0, 0.0, 1.0}};

    double t[3][3]={{0.0, 0.0, 0.0},
                    {0.0, 0.0, 0.0},
                    {0.0, 0.0, 0.0}};

    t[2][1]= v[0];
    t[1][2]=-v[0];
    t[0][2]= v[1];
    t[2][0]=-v[1];
    t[1][0]= v[2];
    t[0][1]=-v[2];       
   
    double phi=2*myrand(mtrand)*PI;
    
    double val1=sin(phi);
    double val2=1-cos(phi);

    //u=val1*t+u+dot(t,t)*val2;
    for(int i=0;i!=3;i++){
        for(int j=0;j!=3;j++){
            u[i][j]=val1*t[i][j]+u[i][j]+ val2*(t[i][0]*t[0][j]+t[i][1]*t[1][j]+t[i][2]*t[2][j]);
            //u[i][1]=val1*t[i][1]+u[i][1]+val2*(t[i][0]*t[0][1]+t[i][1]*t[1][1]+t[i][2]*t[2][1]);
            //u[i][2]=val1*t[i][2]+u[i][2]+val2*(t[i][0]*t[0][2]+t[i][1]*t[1][2]+t[i][2]*t[2][2]);
        }
    }

    //apply rotation matrix u on the molecule 
    double xv,yv,zv; //
    for(int i=0;i!=natoms;i++){
        xv=sposition[i][0]-rotation_point[0];
        yv=sposition[i][1]-rotation_point[1];
        zv=sposition[i][2]-rotation_point[2];

        position[i][0]=rotation_point[0]+u[0][0]*xv+u[0][1]*yv+u[0][2]*zv;
        position[i][1]=rotation_point[1]+u[1][0]*xv+u[1][1]*yv+u[1][2]*zv;
        position[i][2]=rotation_point[2]+u[2][0]*xv+u[2][1]*yv+u[2][2]*zv;
    }


}

//void getCOM(int natoms, double (&com)[3], double *mass, double (*position)[3]){}

/**
 * @brief Get the geometry of a molecule
 * @param natoms number of atoms
 * @param gom return pointer
 * @param position coordinates of the molecules
*/
void getCOG(int natoms, double (&gom)[3], double (*position)[3]){
    
    for(int atom=0;atom!=natoms;atom++){
            gom[0]+=position[atom][0];
            gom[1]+=position[atom][1];
            gom[2]+=position[atom][2];
    }

    for(int i=0;i!=3;i++){
        gom[i]/=natoms;
    }
}


#if DEBUGMOVE
int main(){

    MTRand mtrand;
    double tip3p[3][3]={{ 0.03030173, 0.160511300, -0.354791880},
                        {-0.31833426, 0.561022395,  0.441622365},
                        {0.288032544,-0.721533685, -0.086830479}};
    
    double newtip3p[3][3]={0};    

    double gom[3]={0.0, 0.0, 0.0};

    rotate_mol(3, newtip3p, tip3p, gom, mtrand);

    for(int i=0;i!=3;i++){
        cout<<newtip3p[i][0]<<", "<<newtip3p[i][1]<<", "<<newtip3p[i][2]<< endl;
    }
}
#endif
