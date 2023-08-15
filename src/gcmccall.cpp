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


#include"gcmccall.hpp"
#include"mygcbgcmc.hpp"

using namespace std;
/**
 * @brief read parameter and call gcmc main
 * 
*/
void gcmcio(string parafile){

    int natoms, nRealAtoms;
    int nsteps, nmcmd;
    double mu, rou;
    double xmin, ymin, zmin, xmax, ymax, zmax, dx, dy, dz,rcut;
    double temperature;
    bool drudeff=false;// if the system have drude particles

    //string lbuffer;
    string resname;

    string crdfile;
    string psffile;
    string charmm_template;
    string charmm_execute;

    string linebuffer;
    ifstream fp;

    cout<<" ################################################" <<endl;
    cout<<" #>                 Read inputs                <#" <<endl;
    cout<<" ################################################" <<endl;

    fp.open(parafile,ios::in);
    cout<< " #> Read input file: "<<parafile<<endl;

    fp>>natoms;
    cout<<" #> natoms: "<< natoms <<endl;
    double (*SMOLECULE)[3]=new double[natoms][3];
    string *namelist=new string[natoms];
    string *typelist=new string[natoms];
    fp>>nRealAtoms;
    cout<<" #> nReals: "<< nRealAtoms<<endl;
    fp>>mu;
    //cout<<" #> chemical  potential: "<<mu<<endl;
    cout<<" #> B factor: "<<mu<<endl;

    cout<<" #> read standard molecule"<<endl;
    getline(fp,linebuffer);//need this
    
    for(int i=0;i!=natoms;i++){
        getline(fp,linebuffer);
        SMOLECULE[i][0]=stod( linebuffer.substr( 0,19) );
        SMOLECULE[i][1]=stod( linebuffer.substr(20,39) );
        SMOLECULE[i][2]=stod( linebuffer.substr(40,59) );
        cout<<" #> ["<<SMOLECULE[i][0]<<", "<<SMOLECULE[i][1]<<", "<<SMOLECULE[i][2]<<" ]"<<endl;
    }
    fp>>resname;
    cout<<" #> resname: "<<resname<<endl;
    getline(fp,linebuffer);//need this
    cout<<" #> ATOM name: ";
    for(int i=0;i!=natoms;i++){
        getline(fp,linebuffer);
        //change drude according to the atom name
        if(linebuffer.substr(0,1)=="D"){drudeff=true;}
        namelist[i]=linebuffer.substr(0,4);
        cout<<namelist[i]<<", ";
    }
    cout<<endl;

    cout<<" #> ATOM type: ";
    for(int i=0;i!=natoms;i++){
        getline(fp,typelist[i]);
        //typelist[i]=linebuffer;
        cout<<typelist[i]<<", ";
    }
    cout<<endl;

    getline(fp,crdfile);
    cout<<" #> crd file: "<<crdfile<<endl;
    getline(fp,psffile);
    cout<<" #> psf file: "<<psffile<<endl;

    getline(fp,linebuffer);
    rou=stod(linebuffer);
    cout<<" #> density is: "<< rou<<", But will not use"<<endl;

    getline(fp,charmm_template);
    cout<< " #> charmm template file: "<< charmm_template<<endl;

    getline(fp,charmm_execute);
    cout<< " #> charmm execute file: "<< charmm_execute<<endl;
    //Box information
    getline(fp,linebuffer);
    xmin=stod(linebuffer.substr( 0,9));
    ymin=stod(linebuffer.substr(10,19));
    zmin=stod(linebuffer.substr(20,29));

    getline(fp,linebuffer);
    xmax=stod(linebuffer.substr( 0,9));
    ymax=stod(linebuffer.substr(10,19));
    zmax=stod(linebuffer.substr(20,29));

    getline(fp,linebuffer);
    dx=stod(linebuffer.substr( 0,9));
    dy=stod(linebuffer.substr(10,19));
    dz=stod(linebuffer.substr(20,29));
    
    //cout<<" #> Xmin: "<<xmin<<" Ymin: "<<ymin<<" Zmin: "<<zmin<<endl;;
    //cout<<" #> Xmax: "<<xmax<<" Ymax: "<<ymax<<" Zmax: "<<zmax<<endl;;
    //cout<<" #> Dx:   "<< dx <<" Dy:   "<< dy <<" Dz  : "<<dz<<endl;


    fp>>rcut;
    cout<< " #> r_cut: "<< rcut <<endl;
    fp>>temperature;
    cout<< " #> Temperature: "<< temperature <<endl;
    fp>>nsteps;
    cout<< " #> Run GCMC steps: "<< nsteps <<endl;
    fp>>nmcmd;
    cout<< " #> Run MD each "<<nmcmd<<"  GCMC steps" <<endl;

    fp.close();
    
    cout<<" #> Reading complete" <<endl;

    //GCMCMOLE
    //const double stip3p[3][3]={{ 0.03030173,    0.160511300, -0.35479188},
    //                            {-0.31833426,  0.561022395, 0.441622365},
    //                            { 0.288032544,-0.721533685, -0.086830479}};

    GCMCMOLE gcmcMolecule(natoms, nRealAtoms, mu, SMOLECULE, resname, 
                    typelist, namelist);
    //cout<< "After initialize tip3p gcmcmole"<<endl;
    //charmm class
    //cout<< "Initialize Charmm"<<endl;

    Charmm mysystem(crdfile,psffile,gcmcMolecule,rou,charmm_template,charmm_execute);

    //cout<<"After initialize Charmm class" <<endl;
    //Grid class
    Grid mygrid(xmin, ymin, zmin, xmax, ymax, zmax, dx, dy, dz, rcut);
    //cout<<"After initialize Grid class" <<endl;
    
    
    //GCBGCMC mygcmc(mysystem, mygrid, temperature);
    
    if(drudeff){
        cout<<" #> Find Drude particles!"<<endl;
        cout<<" #> Calling Drude GCBGCMC!"<<endl;
        DGCBGCMC mygcmc(mysystem, mygrid, temperature);
            if (nmcmd==0){
                mygcmc.run(nsteps);
            }else{
                mygcmc.run(nsteps,nmcmd);
            }        
    }else{
        // no drude version
        GCBGCMC mygcmc(mysystem, mygrid, temperature);
            if (nmcmd==0){
                mygcmc.run(nsteps);
            }else{
                mygcmc.run(nsteps,nmcmd);
            }
    }

    //cout<<"After initialize GCBGCMC class" <<endl;

    //cout<<"the end"<<endl;
}
