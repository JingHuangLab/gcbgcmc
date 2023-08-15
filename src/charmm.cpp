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

#include"charmm.hpp"

/**
 * @brief a data structure
 */
void CharmmSystem::clear_all(){
    resid.clear();
    resnumber.clear();
    resname.clear();
    atomtype.clear();
    segid.clear();
    xpos.clear();
    ypos.clear();
    zpos.clear();
}

/**
 * @biref: the charmm class 
 * 
 * 
*/
Charmm::Charmm(string crdfile, string psffile, GCMCMOLE agcmcmol, double density,
                string intputTemplate, string charmmexec):gcmcmol(agcmcmol){
    //cout<< "In Charmm constructor"<<endl;
    charmmEXE=charmmexec;
    
    //dens=density;
    
    inpTemplate=intputTemplate;
    //cout << "BEFORE reading crd"<<endl;
    
    readCrd(crdfile);
    //cout << "after reading crd"<<endl;
    
    string cmd="cp "+inpTemplate+" ./CHARMM.INP";
    system(cmd.c_str());
    
    cmd="cp "+crdfile+" ./accept.crd";
    system(cmd.c_str());

    cmd="cp "+psffile+" ./accept.psf";
    system(cmd.c_str());
}
/**
 * @brief read charmm crd file to get information about the system
 *         should use the extended format
*/
void Charmm::readCrd(string crdfile){
    ifstream fp;
    string lbuffer;
    string resnametmp;
    if(!crd.resname.empty()){crd.clear_all();}//make sure empty
    fp.open(crdfile,ios::in);
    while (!fp.eof()){
        getline(fp,lbuffer);
        if (lbuffer[0]!='*'){
            crd.natoms=stoi(lbuffer.substr(0,10));
        break;
        }
    }
    for(int i=0;i!=crd.natoms;i++){
        getline(fp,lbuffer);
        crd.resnumber.push_back(stoi(lbuffer.substr(10,10)));
        resnametmp=lbuffer.substr(22,8);//for remvoe the space
        crd.resname.push_back(resnametmp.erase(resnametmp.find_last_not_of(" ")+1));
        crd.atomtype.push_back(lbuffer.substr(32,8));
        crd.xpos.push_back(stod(lbuffer.substr(40,20)));
        crd.ypos.push_back(stod(lbuffer.substr(60,20)));
        crd.zpos.push_back(stod(lbuffer.substr(80,20)));
        crd.segid.push_back(lbuffer.substr(102,8));
        crd.resid.push_back(stoi(lbuffer.substr(112,8)));
    }
    fp.close();
}

/**
* @brief update system after accept
*/
void Charmm::update(){
    system("cp try.crd accept.crd");
    system("cp try.psf accept.psf");
    readCrd("accept.crd");

    //if accept clear the tmp for displacement
    backupx.clear(); 
    backupy.clear();
    backupz.clear();

}

/**
* do nothing
*/
void Charmm::restore(){
    return;
}

/*
 *
 */
void Charmm::getAtomIndex(int resid, string resname,int natoms){
    atomIndex.clear();
    int j =0;
    for(int i=0;i!=crd.natoms;i++){//
        if(crd.resid[i]==resid&&crd.resname[i]==resname){
            atomIndex.push_back(i);
            j+=1;
            if(j==natoms){break;}
        }
    }
}
/**
 * @brief Write out crd
 * formate of charmm crd file in extend format:
 * 
 * TITLE
 * NATOM (I10)
 * ATOMNO RESNO   RES  TYPE  X     Y     Z   SEGID RESID Weighting
 *   I10   I10 2X A8 2X A8       3F20.10     2X A8 2X A8 F20.10
 * 
 * 
 */

void Charmm::writeCrd(){
    ofstream fp;
    fp.open("try.crd",ios::out);
    fp<<"* Generate by cppmygcbGCMC\n";
    fp<<"* Have fun with it\n";
    fp<<"*\n";
    fp<<setw(10)<< right <<crd.natoms <<"  "<<"EXT\n";
    for(int i=0;i!=crd.natoms;i++){
        fp<< setw(10) << right << i+1 \
          << setw(10) << right << crd.resnumber[i] << "  " \
          << setw(8) << left << crd.resname[i] << "  " \
          << setw(8) << right << crd.atomtype[i] \
          << setw(20) << fixed << setprecision(10) <<right<< crd.xpos[i] \
          << setw(20) << fixed << setprecision(10) <<right<< crd.ypos[i] \
          << setw(20) << fixed << setprecision(10) <<right<< crd.zpos[i] << "  " \
          << setw(8)  << left << crd.segid[i] <<"  " \
          << setw(8)  << left  << crd.resid[i] \
          << setw(20) << right << setprecision(10) << "0.0000000000" <<endl;
    }

    fp.close();
}

/**
 * @brief Overload of the writeCrd
 * formate of charmm crd file in old format:
 * 
 * @param 
 */

void Charmm::writeCrd(int natoms, double (*position)[3],string resname,
                        string* atomtype){
    ofstream fp;
    fp.open("try.crd",ios::out);
    fp<<"* Generate by CPPmygcbGCMC\n";
    fp<<"* Have fun with it\n";
    fp<<"*\n";
    fp<<setw(10)<< right << crd.natoms+natoms <<"  "<<"EXT\n";
    for(int i=0;i!=crd.natoms;i++){
        fp<< setw(10) << right << i+1 \
          << setw(10) << right << crd.resnumber[i] << "  " \
          << setw(8) << left << crd.resname[i] << "  " \
          << setw(8) << right << crd.atomtype[i] \
          << setw(20) << fixed << setprecision(10) <<right<< crd.xpos[i] \
          << setw(20) << fixed << setprecision(10) <<right<< crd.ypos[i] \
          << setw(20) << fixed << setprecision(10) <<right<< crd.zpos[i] << "  " \
          << setw(8)  << left << crd.segid[i] <<"  " \
          << setw(8)  << left  << crd.resid[i] \
          << setw(20) << right << setprecision(10) << "0.0000000000" <<endl;
    }
    //new moelcules
    //note the width of atomtype
    for(int i=0;i!=natoms;i++){
        fp<< setw(10) << right << crd.natoms+i+1 \
        << setw(10) << right << crd.resnumber[crd.natoms-1]+1 << "  " \
        << setw(8) << left << resname << "  " \
        << setw(8) << left << atomtype[i] \
        << setw(20) << fixed << setprecision(10) <<right <<position[i][0] \
        << setw(20) << fixed << setprecision(10) <<right << position[i][1] \
        << setw(20) << fixed << setprecision(10) <<right << position[i][2] << "  " \
        << setw(8)  << left << "NEW" <<"  " \
        << setw(8)  << left  << 1 \
        << setw(20) << right << setprecision(5) << "0.0000000000" <<endl;
    }

    fp.close();
}


/**
 *@brief Give the resname update the resid list
 */
void Charmm::findResname(string resname){
    residList.clear();
    for(int i=0;i!=crd.natoms;i++){

        if(crd.resname[i]==resname){
            residList.push_back(crd.resid[i]);
        }
    }
    
}

/*
 * if there is GCMC molecules
 */
bool Charmm::ngcmcmols(string resname){
    int n=0;
    for(int i=crd.natoms;i>=0;i--){
        if(crd.resname[i]==resname){
            n+=1;
            if(n/gcmcmol.nrealatoms>1){
                return true;
            }
        }
    }
    return false;
}
//

/*
 * only get the gcmc atom and position
 * except the Drude and 
 */
//void Charmm::get_gcmc_atom_position(){
//    if( !tmpx.empty() ){
//        tmpx.clear();
//        tmpy.clear();
//        tmpz.clear();
//    }
//    for(int i=0;i!=crd.natoms;i++){
//
//    }
//
//}
//

/**
 *@brief  Return the number of real atoms and the position
 *          nreal and the vector
 */
void Charmm::getRealAtomPosition(){
    nReal=0;   
    clearTMPxyz();
    for(int i=0;i!=crd.natoms;i++){
        if(crd.atomtype[i].substr(0,1)=="D"||crd.atomtype[i].substr(0,2)=="LP"){
            continue;
        }
        nReal+=1;
        tmpx.push_back(crd.xpos[i]);
        tmpy.push_back(crd.ypos[i]);
        tmpz.push_back(crd.zpos[i]);
    }
}

/**
 * @brief update the number of real atoms
 */
void Charmm::getNRealAtoms(){
    nReal=0;
    for(int i=0;i!=crd.natoms;i++){
        if(crd.atomtype[i].substr(0,1)=="D"||crd.atomtype[i].substr(0,2)=="LP"){
            continue;
        }
        nReal+=1;
    }
}

/**
 * @brief Write out crd file and call charmm to return energy
 * @param insertPos the insertPosiont
 * @return energy
*/
double Charmm::insertion(int natoms, double (*insertPos)[3],
                            string resname,string* atomtype){

    //write crd and call charmm
    writeCrd(natoms, insertPos, resname, atomtype);

    systemEnergy=energy(2);

    return systemEnergy;
}
/*
 * do nothing
 */
void Charmm::restoreInsertion(){return;}

/**
 * @brief 
 *
 */
double Charmm::deletion(int natoms, int resid, string resname){
    string seleCmd="RESID "+to_string(resid)+" .and. resname "+resname;
    string cmd="sed -e s/SELECMD/\""+seleCmd+"\"/g "+inpTemplate+" > CHARMM.INP";
    system(cmd.c_str());
    getAtomIndex(resid,resname,gcmcmol.natoms);
    clearTMPxyz();
    for(int i=0;i!=natoms;i++){
        tmpx.push_back(crd.xpos[atomIndex[i]]);
        tmpy.push_back(crd.ypos[atomIndex[i]]);
        tmpz.push_back(crd.zpos[atomIndex[i]]);
    }
    return energy(3);
}

/**
 * @brief nothing to do because the crd is not changed!
 */
void Charmm::restore_deletion(){return;}
//

/**
 * @brief Do displacement in GCMC
 * @param natoms number of atoms in GCMC molecule
 * @param iresid the random selected resid
 * @param newPos the new position of the select resid
 * @return energy
*/
double Charmm::displacement(int natoms, int iresid, double (*newPos)[3]){
    
    //random choice a resid 
    getAtomIndex(iresid, gcmcmol.resn, gcmcmol.natoms);
    //backup the original and replace
    
    int j;
    for(int i=0;i!=natoms;i++){
        j=atomIndex[i];
        backupx.push_back(crd.xpos[j]);
        backupy.push_back(crd.ypos[j]);
        backupz.push_back(crd.zpos[j]);

        crd.xpos[j]=newPos[i][0];
        crd.ypos[j]=newPos[i][1];
        crd.zpos[j]=newPos[i][2];

    }
    //write out 
    writeCrd();

    return energy(1);
}
/**
 *@brief restore the atom position after reject displacement 
 *
 */
void Charmm::restoreDisplacement(){
    int j;
    for(int i=0;i!=gcmcmol.natoms;i++){
        j=atomIndex[i];

        crd.xpos[j]=backupx[i];
        crd.ypos[j]=backupy[i];
        crd.zpos[j]=backupz[i];
    }
    
    //if reject then clear the template 
    backupx.clear(); 
    backupy.clear();
    backupz.clear();
}
/**
 * @brief read energy from a file
*/
double inline Charmm::readEnergy(string filename){
    double energy;
    //read energy
    ifstream infile;
    infile.open("energy.dat");
    infile >>energy;
    infile.close();
    return energy;
}

/**
 * @brief call charmm and read the energy
 * @param jobtype   0: single point energy
 *                  1: displacement
 *                  2: insertion
 *                  3: deletion
 *                  5: run MD simulations
 * @return system energy
 */
double Charmm::energy(int jobtype){
    string exe_cmd;
    exe_cmd=charmmEXE+" jobtype="+to_string(jobtype)+" -i CHARMM.INP -o CHARMM.OUT";
    system(exe_cmd.c_str());

    //check 
    check_output();

    //double energy;
    ////read energy
    //ifstream infile;
    //infile.open("energy.dat");
    //infile >>energy;
    //infile.close();

    return readEnergy("energy.dat");
}

/**
 * @brief call charmm to run MD simulations
 * 
 * @return system energy
*/
double Charmm::runMD(){

    system("cp accept.crd beforemd.crd");
    string call_charmm_cmd=charmmEXE+" jobtype=5 -i CHARMM.INP -o CHARMM.OUT";
    system(call_charmm_cmd.c_str());
    //check
    check_output();

    //system("cp try.crd accept.crd"); we will do this in the update()
    
    return readEnergy("energy.dat");
}

/**
*@brief Test if charmm run properly
*/
void Charmm::check_output(){
    if(!system("grep --silent ABNORMAL CHARMM.OUT")){
        cout<<"Error! charmm does not run properly!!!"<<endl;
        exit(1);
    }
}


#if DEBUGCHARMM
int main(){
    
    double dens=0.03342;
    string mytemplate="energy.inp";
    Charmm mycharmm("onewater.crd", "onewater.psf", dens, mytemplate,"charmm");

    double energy=mycharmm.energy();
    cout << energy << endl;

}
#endif
