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


#include"mygcbgcmc.hpp"
#include"myrand.hpp"


//use the constructor initializer list to initialize the three class
GCBGCMC::GCBGCMC(Charmm &asystem, Grid &agrid, double temperature):
                    mysystem(asystem), mygrid(agrid){


    beta=1.0/(kB*temperature);
    Bfactor=exp(mysystem.gcmcmol.mu);//here use the input mu as B factor

    //cout<<"beta= "<<beta<<endl;
    
    mygrid.initialCavityList(mysystem.crd.natoms,
                                mysystem.crd.xpos,
                                mysystem.crd.ypos,
                                mysystem.crd.zpos);
    cout<<" #> Initialization complete"<<endl;
    
    old_energy=mysystem.energy();

    insertPos=new double [mysystem.gcmcmol.natoms][3];
    deletePos=new double [mysystem.gcmcmol.natoms][3];


    system("mkdir -p ./gcmc_traj");

}

/**
 * @brief the Metropolis accptance method
 * @param dE the inverse energy difference defined as E_old-E_new which contain the chemical potential
 * @param cbfactor the cavity biased factor due
 * @return true(accept) or false(reject)
*/
bool GCBGCMC::acceptance(double dE, double cbfactor){
    string backup;
    //backup="cp try.crd ./gcmc_traj/"+to_string(nstep)+".crd";
    //system(backup.c_str());
    //backup="cp try.psf ./gcmc_traj/"+to_string(nstep)+".psf";
    //system(backup.c_str());
    
    cout << "dE: "<<dE<<", cbfactor: "<<cbfactor;
    /*if (dE>0){
        cout<< " Acc"<<endl;

        return true;
    }*/ 
    
    if (cbfactor*exp(beta*dE)> myrand(mtrand)){
        cout<< " Acc"<<endl;
        return true;
    }
    else{
        cout<< " Rej"<<endl;
        return false;
    }

}

/**
 * @brief print out the statistic information during the simulations
 * 
*/
void GCBGCMC::statastic() const{
    cout<<"------------------------------"<<endl;
    cout<<"Dis: "<<ndisplacement<<", accept: "<<displacement_accept<<endl;
    cout<<"Ins: "<<ninsertion   <<", accept: "<<insertion_accept   <<endl;
    cout<<"Del: "<<ndeletion    <<", accept: "<<deletion_accept    <<endl;
    cout<<"------------------------------"<<endl;

}

/**
 * @brief displacement
 * 
*/
void GCBGCMC::displacement(){

    //cout<<"Begin displacement"<<endl;
    mysystem.findResname(mysystem.gcmcmol.resn);//a list of number
    int iresid, icell;
    iresid=choice(mysystem.residList, mtrand);
    mygrid.getCavities();//first update ncavites

    if(mygrid.ncavities){
        //cout<< "ncavities: "<<mygrid.ncavities<<endl;
        //mygrid.getCavities();
        //cout<<"size of cavityList: "<<mygrid.cavityList.size()<<endl;
        
        icell=choice(mygrid.cavityList, mtrand);

        mygrid.getCellBoundary(icell); 
        mysystem.gcmcmol.new_position(mygrid.cellBoundary[0],
                    mygrid.cellBoundary[1],mygrid.cellBoundary[2],
                    mygrid.dxcav, mygrid.dycav, mygrid.dzcav,
                    mysystem.gcmcmol.natoms, insertPos, mtrand);
        //cout<<endl;
        //for(int i =0; i!=3;i++){
        //cout<< insertPos[i][0]<<","<<insertPos[i][1]<<","<<insertPos[i][2]<<endl;
        //}

        new_energy=mysystem.displacement(mysystem.gcmcmol.natoms, iresid, insertPos);

        for(int i=0;i!=mysystem.gcmcmol.natoms;i++){
            deletePos[i][0]=mysystem.backupx[i];
            deletePos[i][1]=mysystem.backupy[i];
            deletePos[i][2]=mysystem.backupz[i];
        } 
    }
    else{
        mysystem.gcmcmol.new_position(mygrid.xmin,mygrid.ymin, mygrid.zmin,
                    mygrid.xlen, mygrid.ylen, mygrid.zlen,
                    mysystem.gcmcmol.natoms, insertPos,mtrand);
                    
        new_energy=mysystem.displacement(mysystem.gcmcmol.natoms, iresid, insertPos);
        
        for(int i=0;i!=mysystem.gcmcmol.natoms;i++){
            deletePos[i][0]=mysystem.backupx[i];
            deletePos[i][1]=mysystem.backupy[i];
            deletePos[i][2]=mysystem.backupz[i];
        } 
    }

    double cbfactor=1.0;
    if(acceptance(old_energy-new_energy,cbfactor)){
        displacement_accept+=1;
        mysystem.update();
        // must update insertion first or negative occupancy will occur
        //mygrid.updateInsertion(mysystem.gcmcmol.natoms, insertPos,nstep);
        //mygrid.updateDeletion(mysystem.gcmcmol.natoms, deletePos,nstep);
        mygrid.updateDisplacement(mysystem.crd.natoms,
                                mysystem.crd.xpos,
                                mysystem.crd.ypos,
                                mysystem.crd.zpos);

        old_energy=new_energy;
    }
    else{
        mysystem.restoreDisplacement();       
    }
    //cout<<"End displacement"<<endl;

}

/**
 * @brief do deletion instance
 * 
 * 
*/
void GCBGCMC::deletion(){
    //cout<<"Begin deletion"<<endl;

    if(mysystem.ngcmcmols(mysystem.gcmcmol.resn)){// if there is GCMC molecules
        
        mysystem.findResname(mysystem.gcmcmol.resn);//a list of number 
        int iresid;
        iresid=choice(mysystem.residList, mtrand);

        new_energy=mysystem.deletion(mysystem.gcmcmol.natoms,iresid,mysystem.gcmcmol.resn);
        //cout<<"In deletion call"<<endl;
        for(int i=0;i<mysystem.gcmcmol.natoms;i++){
            deletePos[i][0]=mysystem.tmpx[i];
            deletePos[i][1]=mysystem.tmpy[i];
            deletePos[i][2]=mysystem.tmpz[i];
            //cout<<deletePos[i][0]<<", "<<deletePos[i][1]<<", "<<deletePos[i][2]<<endl;
        } 

        
        //update grid before acceptance
        mygrid.updateDeletion(mysystem.gcmcmol.natoms, deletePos, nstep);

        mysystem.getNRealAtoms();
        double cbfactor=mysystem.nReal/(Bfactor*mygrid.calProbcavities());
        //double cbfactor=mysystem.nReal/(mysystem.dens*mygrid.V*mygrid.calProbcavities());
        //double cbfactor=mysystem.nReal/(3*0.23820*0.23820*0.23820*mygrid.V*mygrid.calProbcavities());


        //if(acceptance(old_energy-new_energy-mysystem.gcmcmol.mu, cbfactor)){
        if(acceptance(old_energy-new_energy, cbfactor)){

            deletion_accept+=1;
            mysystem.update();
            old_energy=new_energy;
        }else{
            mysystem.restore_deletion();
            mygrid.restoreDeletion();

        }
    }
    //cout<<"End deletion"<<endl;

}
/**
 * @brief Insertion
 * 
 * 
*/
void GCBGCMC::insertion(){

    // if cavity exists
    //cout<<"Begin Insertion"<<endl;
    mygrid.getCavities();
    double cbfactor=1.0;
    if(mygrid.ncavities){
        mygrid.getCavities();
        int icell=choice(mygrid.cavityList, mtrand);

        mygrid.getCellBoundary(icell);
        mysystem.gcmcmol.new_position(mygrid.cellBoundary[0],
                    mygrid.cellBoundary[1],mygrid.cellBoundary[2],
                    mygrid.dxcav, mygrid.dycav, mygrid.dzcav,
                    mysystem.gcmcmol.natoms, insertPos, mtrand);

        //cout<<"In insertion call"<<endl;

        //for(int j=0;j!=3;j++){
        //    cout<<insertPos[j][0]<<","<<insertPos[j][1]<<","<<insertPos[j][2]<<endl;
        //}
        
        new_energy=mysystem.insertion(mysystem.gcmcmol.natoms, insertPos,
                                        mysystem.gcmcmol.resn,mysystem.gcmcmol.atomNameList);
        mysystem.getNRealAtoms();
        
        cbfactor=Bfactor*mygrid.calProbcavities()/(
                        mysystem.nReal+mysystem.gcmcmol.nrealatoms);
        //cbfactor=mysystem.dens*mygrid.V*mygrid.calProbcavities()/(
        //                mysystem.nReal+mysystem.gcmcmol.nrealatoms);
        //cbfactor=3*0.23820*0.23820*0.23820*mygrid.V*mygrid.calProbcavities()/(
        //                mysystem.nReal+mysystem.gcmcmol.nrealatoms);
    }
    else{

        mysystem.gcmcmol.new_position(mygrid.xmin, mygrid.ymin, mygrid.zmin,
                    mygrid.xlen, mygrid.ylen, mygrid.zlen,
                    mysystem.gcmcmol.natoms, insertPos,mtrand);
        
        new_energy=mysystem.insertion(mysystem.gcmcmol.natoms, insertPos,
                                        mysystem.gcmcmol.resn,mysystem.gcmcmol.atomNameList);
    }

    //for(int i=0;i!=3;i++){
    //    cout<<insertPos[i][0]<<", "<<insertPos[i][1]<<", "<<insertPos[i][2]<<endl;
    //}
    
    double dE=old_energy-new_energy;//+mysystem.gcmcmol.mu; //times -1
    
    if(acceptance(dE,cbfactor)){
        insertion_accept+=1;
        mysystem.update();
        
        //In case charmm deal with the PBC and get 
        //a image coordination other than the insertPos, 
        //we instead use the coordination from charmm output
        //to updat the grid occupancy 
        int atomindex;
        //cout<< "Real used"<<endl;
        for(int atom=0;atom!=mysystem.gcmcmol.natoms;atom++){
            atomindex=mysystem.crd.natoms-mysystem.gcmcmol.natoms+atom;
            insertPos[atom][0]=mysystem.crd.xpos[atomindex];
            insertPos[atom][1]=mysystem.crd.ypos[atomindex];
            insertPos[atom][2]=mysystem.crd.zpos[atomindex];
            //cout<<insertPos[atom][0]<<","<<insertPos[atom][1]<<","<<insertPos[atom][2]<<endl;
        }
        mygrid.updateInsertion(mysystem.gcmcmol.natoms, insertPos, nstep);

        old_energy=new_energy;
    }else{
        mysystem.restoreInsertion();
    }
    //cout<<"End Insertion"<<endl;
}

/**
 * @brief do GCMC
 * @param nsteps number of steps to perform GCMC
*/
void GCBGCMC::run(int nsteps){
    cout<<" #> GCMC RUN START\n"<< endl;
    int movetype;
    for(int i=0; i<nsteps; i++){
        nstep=i+1;
        cout<<"Step: "<<nstep;
        movetype=randmove(mtrand);
        //movetype=0;

        switch (movetype){
            
        case 0://displacement
            cout<<" Dis ";
            ndisplacement+=1;
            displacement();
            break;
        case 1://insertion
            cout<<" Ins ";
            ninsertion+=1;
            insertion();
            break;
        case 2://deletion
            cout<<" Del ";
            ndeletion+=1;
            deletion();
            break;
        default:
            cout<<"Unknown move type: "<<movetype<<endl;
            exit(1);
        }
    }
    statastic();
}
        
/**
 * @brief do GCMC/MD
 * @param nsteps number of steps to perform GCMC
 * @param nmc perform MD simulations after each nmc GCMC steps, 
 *              how many MD steps are defined in the charmm template file 
*/
void GCBGCMC::run(int nsteps, int nmc){
    cout<<" #> GCMC/MD RUN START\n"<< endl;
    int movetype;
    for(int i=0; i<nsteps; i++){
        nstep=i+1;
        cout<<"Step: "<<nstep;
        movetype=randmove(mtrand);
        //movetype=0;

        switch (movetype){
            
        case 0://displacement
            cout<<" Dis ";
            ndisplacement+=1;
            displacement();
            break;
        case 1://insertion
            cout<<" Ins ";
            ninsertion+=1;
            insertion();
            break;
        case 2://deletion
            cout<<" Del ";
            ndeletion+=1;
            deletion();
            break;
        default:
            cout<<"Unknown move type: "<<movetype<<endl;
            exit(1);
        }

        //run MD
        if(nstep%nmc==0){
            //copy the accept as trajectory of the simulation
            string backup;
            backup="cp accept.crd ./gcmc_traj/"+to_string(nstep)+".crd";
            system(backup.c_str());
            backup="cp accept.psf ./gcmc_traj/"+to_string(nstep)+".psf";
            system(backup.c_str());
            
            cout<< " #> "<< nstep<< " steps, perform MD"<<endl;
            new_energy=mysystem.runMD();
            mysystem.update();

            mygrid.initialCavityList(mysystem.crd.natoms,
                                mysystem.crd.xpos,
                                mysystem.crd.ypos,
                                mysystem.crd.zpos);

        }
    }
    statastic();
}
//GCBGCMC::~GCBGCMC(){
//    cout<<"GCBGCMC destructor: "<<endl;
//
//    //delete [] insertPos;
//    //delete [] deletePos;
//
//    cout<<"finish GCBGCMC destructor: "<<endl;
//}

/**
 *@brief update grid for every particles after acceptance 
*/
void DGCBGCMC::dinsertion(){

    mygrid.getCavities();
    double cbfactor=1.0;
    if(mygrid.ncavities){
        mygrid.getCavities();
        int icell=choice(mygrid.cavityList, mtrand);

        mygrid.getCellBoundary(icell);
        mysystem.gcmcmol.new_position(mygrid.cellBoundary[0],
                    mygrid.cellBoundary[1],mygrid.cellBoundary[2],
                    mygrid.dxcav, mygrid.dycav, mygrid.dzcav,
                    mysystem.gcmcmol.natoms, insertPos, mtrand);
        
        new_energy=mysystem.insertion(mysystem.gcmcmol.natoms, insertPos,
                                        mysystem.gcmcmol.resn,mysystem.gcmcmol.atomNameList);
        mysystem.getNRealAtoms();
        
        cbfactor=Bfactor*mygrid.calProbcavities()/(
                        mysystem.nReal+mysystem.gcmcmol.nrealatoms);
        //cbfactor=mysystem.dens*mygrid.V*mygrid.calProbcavities()/(
        //                mysystem.nReal+mysystem.gcmcmol.nrealatoms);
        //cbfactor=3*0.23820*0.23820*0.23820*mygrid.V*mygrid.calProbcavities()/(
        //                mysystem.nReal+mysystem.gcmcmol.nrealatoms);
    }
    else{

        mysystem.gcmcmol.new_position(mygrid.xmin, mygrid.ymin, mygrid.zmin,
                    mygrid.xlen, mygrid.ylen, mygrid.zlen,
                    mysystem.gcmcmol.natoms, insertPos,mtrand);
        
        new_energy=mysystem.insertion(mysystem.gcmcmol.natoms, insertPos,
                                        mysystem.gcmcmol.resn,mysystem.gcmcmol.atomNameList);
    }

    
    double dE=old_energy-new_energy;//+mysystem.gcmcmol.mu; //times -1
    
    if(acceptance(dE,cbfactor)){
        insertion_accept+=1;
        mysystem.update();
        mygrid.initialCavityList(mysystem.crd.natoms,
                                mysystem.crd.xpos,
                                mysystem.crd.ypos,
                                mysystem.crd.zpos);
        
        old_energy=new_energy;
    }else{
        mysystem.restoreInsertion();
    }

}

/**
 * @brief update grid after accept, the first update is enough for acceptance
*/
void DGCBGCMC::ddeletion(){

    if(mysystem.ngcmcmols(mysystem.gcmcmol.resn)){// if there is GCMC molecules
        
        mysystem.findResname(mysystem.gcmcmol.resn);//a list of number 
        int iresid;
        iresid=choice(mysystem.residList, mtrand);

        new_energy=mysystem.deletion(mysystem.gcmcmol.natoms,iresid,mysystem.gcmcmol.resn);
        //cout<<"In deletion call"<<endl;
        for(int i=0;i<mysystem.gcmcmol.natoms;i++){
            deletePos[i][0]=mysystem.tmpx[i];
            deletePos[i][1]=mysystem.tmpy[i];
            deletePos[i][2]=mysystem.tmpz[i];
            //cout<<deletePos[i][0]<<", "<<deletePos[i][1]<<", "<<deletePos[i][2]<<endl;
        } 

        
        //update grid before acceptance
        mygrid.updateDeletion(mysystem.gcmcmol.natoms, deletePos, nstep);

        mysystem.getNRealAtoms();
        double cbfactor=mysystem.nReal/(Bfactor*mygrid.calProbcavities());
        //double cbfactor=mysystem.nReal/(mysystem.dens*mygrid.V*mygrid.calProbcavities());
        //double cbfactor=mysystem.nReal/(3*0.23820*0.23820*0.23820*mygrid.V*mygrid.calProbcavities());


        //if(acceptance(old_energy-new_energy-mysystem.gcmcmol.mu, cbfactor)){
        if(acceptance(old_energy-new_energy, cbfactor)){

            deletion_accept+=1;
            mysystem.update();
            mygrid.initialCavityList(mysystem.crd.natoms,
                                mysystem.crd.xpos,
                                mysystem.crd.ypos,
                                mysystem.crd.zpos);
            
            old_energy=new_energy;
        }else{
            mysystem.restore_deletion();
            mygrid.restoreDeletion();

        }
    }
}

/**
 * @brief do GCMC
 * @param nsteps number of steps to perform GCMC
*/
void DGCBGCMC::run(int nsteps){
    cout<<" #> GCMC RUN START\n"<< endl;
    int movetype;
    for(int i=0; i<nsteps; i++){
        nstep=i+1;
        cout<<"Step: "<<nstep;
        movetype=randmove(mtrand);
        //movetype=0;

        switch (movetype){
            
        case 0://displacement
            cout<<" Dis ";
            ndisplacement+=1;
            displacement();
            break;
        case 1://insertion
            cout<<" Ins ";
            ninsertion+=1;
            dinsertion();
            break;
        case 2://deletion
            cout<<" Del ";
            ndeletion+=1;
            ddeletion();
            break;
        default:
            cout<<"Unknown move type: "<<movetype<<endl;
            exit(1);
        }
    }
    statastic();
}
        
/**
 * @brief do GCMC/MD
 * @param nsteps number of steps to perform GCMC
 * @param nmc perform MD simulations after each nmc GCMC steps, 
 *              how many MD steps are defined in the charmm template file 
*/
void DGCBGCMC::run(int nsteps, int nmc){
    cout<<" #> GCMC/MD RUN START\n"<< endl;
    int movetype;
    for(int i=0; i<nsteps; i++){
        nstep=i+1;
        cout<<"Step: "<<nstep;
        movetype=randmove(mtrand);

        switch (movetype){
            
        case 0://displacement
            cout<<" Dis ";
            ndisplacement+=1;
            displacement();
            break;
        case 1://insertion
            cout<<" Ins ";
            ninsertion+=1;
            dinsertion();
            break;
        case 2://deletion
            cout<<" Del ";
            ndeletion+=1;
            ddeletion();
            break;
        default:
            cout<<"Unknown move type: "<<movetype<<endl;
            exit(1);
        }

        //run MD
        if(nstep%nmc==0){
            //copy the accept as trajectory of the simulation
            string backup;
            backup="cp accept.crd ./gcmc_traj/"+to_string(nstep)+".crd";
            system(backup.c_str());
            backup="cp accept.psf ./gcmc_traj/"+to_string(nstep)+".psf";
            system(backup.c_str());
            
            cout<< " #> "<< nstep<< " steps, perform MD"<<endl;
            new_energy=mysystem.runMD();
            mysystem.update();

            mygrid.initialCavityList(mysystem.crd.natoms,
                                mysystem.crd.xpos,
                                mysystem.crd.ypos,
                                mysystem.crd.zpos);

        }
    }
    statastic();
}

#if DEBUGMAIN
int main(){

}
#endif
