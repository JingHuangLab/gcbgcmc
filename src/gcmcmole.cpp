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

#include"gcmcmole.hpp"
#include"gcmcmove.hpp"

/**
 * @brief construct the GCMCMOLE class
 * @param natoms number of atoms in the GCMC molecule
 * @param nrealatoms number of real atoms in the GCMC molecule, 
 *      e.g. exclude lonepair, drude
 * @param mu the chemical potential of each GCMC molecule
 * @param molecule the standard coordinates of the GCMC molecule 
 *      with geometry at [0.0, 0.0, 0.0]
 * @param resname the residue name of the GCMC molecule
*/

GCMCMOLE::GCMCMOLE(int Natoms, int nrealatoms, double aMu, const double (*molecule)[3],
                string resname, string* typelist,
                string* namelist):natoms(Natoms),nrealatoms(nrealatoms),smolecule(molecule){

                mu=aMu;
                
                //smolecule=molecule;
                resn=resname;
                //atomMassList=masslist;
                atomNameList=namelist;
                atomTypeList=typelist;
                
                }

/**
 * @brief Generate a new position of the GCMC molecule
*/
void GCMCMOLE::new_position(double xmin, double ymin, double zmin, 
                        double xlen, double ylen, double zlen,
                        int natoms, double (*position)[3], MTRand &mtrand){

    //rotate fist 
    rotate_mol(natoms, position, smolecule, GOM, mtrand);
    //cout<<endl;
    //for(int j=0;j!=3;j++){
    //        cout<<position[j][0]<<","<<position[j][1]<<","<<position[j][0]<<endl;
    //}


    //translate
    translation(xmin, ymin, zmin, 
                xlen, ylen, zlen, natoms, position, mtrand);

}


#if DEBUGGCMCMOL
int main(){
    int natoms=3,nRealAtoms=3;
    double mu=-5.80;
    //string linebuffer;
    string resname="TIP3";
    double masslist[3]={16.0, 1.0, 1.0};
    string typelist[3]={"OT", "HT","HT"};
    string namelist[3]={"OH2", "H1", "H2"};
    double stip3p[3][3]={{0.03030173,  0.1605113,  -0.35479188},
                        {-0.31833426,  0.561022395, 0.441622365},
                        { 0.288032544,-0.721533685, -0.086830479}};
    
    GCMCMOLE a(natoms, nRealAtoms, mu,stip3p,resname,
                masslist,typelist,namelist);
    

        double b[3][3]={0};


        
        a.new_position(0.0, 0.0, 0.0, 5.0, 5.0, 5.0, 3, b);
        
        for(int i=0;i!=3;i++){
                
                cout<< b[i][0] <<"\t"<< b[i][1]<<"\t"<<b[i][2]<<endl;

        }
        cout <<"destruct"<<endl;
}
#endif