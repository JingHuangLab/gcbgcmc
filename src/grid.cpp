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

#include "grid.hpp"

/**
 * @brief:
 * 
 * 
 * 
*/
Grid::Grid(double x_min, double y_min, double z_min, 
            double x_max, double y_max, double z_max, 
            double dx, double dy, double dz, double r){
    //
    cout<<endl; 
    cout<<" ################################################" <<endl;
    cout<<" #>                    Grid                    <#" <<endl;
    cout<<" #>         OpenMP parallel supported!         <#" <<endl;
    cout<<" ################################################" <<endl;
    #pragma omp parallel
    {
        if (omp_get_thread_num()==0){
            cout<<" #> Use "<<omp_get_num_threads()<<" threads"<<endl;
        }
    }
    cout<<"\n #> Initialize the grid" <<endl;
    
    // 
    xmin=x_min; xmax=x_max;
    ymin=y_min; ymax=y_max;
    zmin=z_min; zmax=z_max;
    cout<<" #> Grid boundary X: "<< xmin <<" --- "<<xmax<< endl;
    cout<<" #> Grid boundary Y: "<< ymin <<" --- "<<ymax<< endl;
    cout<<" #> Grid boundary Z: "<< zmin <<" --- "<<zmax<< endl;
    //
    xlen=xmax-xmin;
    ylen=ymax-ymin;
    zlen=zmax-zmin;
    hxlen=0.5*xlen;
    hylen=0.5*ylen;
    hzlen=0.5*zlen;
    //
    V=xlen*ylen*zlen;
    nx=int(xlen/dx);
    ny=int(ylen/dy);
    nz=int(zlen/dz);
    cout<<" #> Grid on X: "<< nx << endl;
    cout<<" #> Grid on Y: "<< ny << endl;
    cout<<" #> Grid on Z: "<< nz << endl;
    //
    ncells=nx*ny*nz;
    cout<<" #> Number of cells is "<<ncells<<endl;
    
    dxcav=xlen/nx;
    dycav=ylen/ny;
    dzcav=zlen/nz;
    hdxcav=0.5*dxcav;
    hdycav=0.5*dycav;
    hdzcav=0.5*dzcav;
    cout<<" #> Grid length on X: "<< dxcav << endl;
    cout<<" #> Grid length on Y: "<< dycav << endl;
    cout<<" #> Grid length on Z: "<< dzcav << endl;

    double gridx=xmin+hdxcav;
    for(int i=0;i<nx;i++){
        double gridy=ymin+hdycav;
        for(int j=0;j<ny;j++){
            double gridz=zmin+hdzcav;
            for(int k=0;k<nz;k++){
                xcell_center.push_back(gridx);
                ycell_center.push_back(gridy);
                zcell_center.push_back(gridz);
                gridz+=dzcav;
            }
        gridy+=dycav;
        }
    gridx+=dxcav;
    }
    //> debug
    //ofstream fp;
    //fp.open("gridcenter.txt",ios::out);
    //for(int i=0;i!=ncells;i++){
    //    fp<<xcell_center[i]<<", "<<ycell_center[i]<<", "<<zcell_center[i]<<endl;
    //}
    //fp.close();

    // user defined cut off
    r_cut=r;
    r_cut2=r_cut*r_cut;

    //
    xsdepth=int(r/dxcav)+1;
    ysdepth=int(r/dycav)+1;
    zsdepth=int(r/dzcav)+1;
    cout<<" #> Grid update depth on X: "<<xsdepth<<endl;
    cout<<" #> Grid update depth on Y: "<<ysdepth<<endl;
    cout<<" #> Grid update depth on Z: "<<zsdepth<<endl;

}

/**
 * @brief Occupied is defined as distance less (eqal) than rcut/
 *          If the cell is occupied by the atom
 *          If one direction is larger that rcut then atom cannot occupanied
 *          this 
 *        Sometimes, we only want to apply pbc to some of the directions
 *           
 * 
*/
bool Grid::ifOccupied(double &x0, double &y0, double &z0, 
                    double &x1, double &y1, double &z1,
                    bool pbcx, bool pbcy, bool pbcz) const{
    
    double dx = fabs(x0 - x1);
    if(dx > hxlen){dx -= xlen;}
    //> when no pbc on x,
    //> commet out the up line, uncomment the following ling
    if(dx >= r_cut){return false;}

    double dy = fabs(y0 - y1);
    if(dy > hylen){dy -= ylen;}
    //> when no pbc on y, 
    //> commet out the up line, uncomment the following ling
    if(dy >= r_cut){return false;}

    double dz = fabs(z0 - z1);
    if(dz > hzlen){dz -= zlen;}
    //> when no pbc on z, 
    //> commet out the up line, uncomment the following ling
    if(dz >= r_cut){return false;}

    if(dx*dx + dy*dy + dz*dz > r_cut2){
        return false;
    }

    return true;

}

/**
 * @brief set attribute cellBoundary according to the cell index
 * 
 * @param icell cell index
*/
void Grid::getCellBoundary(int icell){
    //only the lower point is important
    cellBoundary[0]=xcell_center[icell]-hdxcav;
    cellBoundary[1]=ycell_center[icell]-hdycav;
    cellBoundary[2]=zcell_center[icell]-hdzcav;
}

/**
 *@brief return the cell index given given the cell position
 */
inline int Grid::cellPos2Index(int i, int j, int k){
    return i*ny*nz+j*nz+k; 
}

/**
 * @brief Initialize and update the occupancy and ncavities
*/
void Grid::initialCavityList(int natoms, vector<double> &xpos, 
                    vector<double> &ypos, vector<double> &zpos){
    
    //initialize the occupancy with 0 and ccavities with false
    ncavities=0;
    occupancy=vector<int>(ncells,0);
    cavities=vector<int>(ncells,0);
    #pragma omp parallel for
    for(int icell=0;icell<ncells;icell++){
        //cout<<"x0: "<<xpos[i]<<"y0: "<<ypos[i]<<"z0: "<<zpos[i]<<endl;
        // openmp maybe work here
        for(int atom=0;atom<natoms;atom++){
            if(ifOccupied(xpos[atom], ypos[atom], zpos[atom],
                xcell_center[icell], ycell_center[icell], zcell_center[icell])){
                
                occupancy[icell]+=1;
            }
        }
    }
    
    for(int icell=0;icell<ncells;icell++){
        
        if(occupancy[icell]==0){
            cavities[icell]=1;
            ncavities+=1;
        }
    }

    cout<<" #> Availiable cavities: "<<ncavities<<"/"<<ncells<<endl;
    cout<<" #> Grid update complete!\n"<<endl;

}

/**
 * @brief calculate the rearch area index range along a direction with PBC
 * @param index the cell index along a direction, starts from 0
 * @param nd one of nx,ny,nz
 * @param sd search depth along that direction
 * @param das a vector container to put the results
 * problem: will contain multi times
*/
void Grid::directionIndex(int index, int nd, int sd, vector<int> &dsa){

    //dsa.clear();
    if(nd<=2*sd+1){//return all the length
        for(int i=0;i!=nd;i++){
            dsa.push_back(i);
        }
        return;
    }

    if(index<sd){ //go beyond the lower
        for(int i=0;i!=index+sd+1;i++){
            dsa.push_back(i);
        }
        for(int i=nd-(sd-index);i!=nd;i++){
            dsa.push_back(i);
        }
    }
    else if(index > (nd-sd-1)){ //go beyond the upper
        for(int i=index-sd;i!=nd+1;i++){
            dsa.push_back(i);
        }
        for(int i=0;i!=sd-(nd-index)+1;i++){
            dsa.push_back(i);
        }
    }
    else{//normal
        for(int i=index-sd;i!=index+sd+1;i++){
            dsa.push_back(i);
        }
    } 

}

/**
 * @brief Given the Nx,Ny,Nz of the center cell
 *          Note doing this for an atom each time
 * 
 * 
 */
void Grid::updateArea(double *changedPosition){
    int alongx,alongy,alongz;
    vector<int> xsa,ysa,zsa;
    //since the type of alongx is int,
    //so the results is int automatically
    alongx=(changedPosition[0]-xmin)/dxcav;
    alongy=(changedPosition[1]-ymin)/dycav;
    alongz=(changedPosition[2]-zmin)/dzcav;

    directionIndex(alongx,nx,xsdepth,xsa);//x direction
    //xsa.insert(xsa.end(),dsa.begin(),dsa.end());
        
    directionIndex(alongy,ny,ysdepth,ysa);//y direction
    //ysa.insert(ysa.end(),dsa.begin(),dsa.end());
        
    directionIndex(alongz,nz,zsdepth,zsa);//z direction
    //zsa.insert(zsa.end(),dsa.begin(),dsa.end());

    //get the update Area
    searchList.clear();
    for(vector<int>::iterator xiter=xsa.begin();
        xiter!=xsa.end();xiter++){
        for(vector<int>::iterator yiter=ysa.begin();
            yiter!=ysa.end();yiter++){
            for(vector<int>::iterator ziter=zsa.begin();
                ziter!=zsa.end();ziter++){
                    searchList.push_back(*xiter);
                    searchList.push_back(*yiter);
                    searchList.push_back(*ziter);
                }
        }
    }

}

/**
 * @brief update grid cavity information after insertion
 * @param natoms number of atoms
 * @param changedPosition the changed coordinates
*/
void Grid::updateInsertion(int natoms, double (*changedPosition)[3],int n){
    //int icell;
    //first get the search area
    //cout<<"run update insertion"<<endl;
    //writeOccupancy("Occinsb",n);
    //cout<<"Insert Pos: "<<endl;
    #pragma omp parallel for
    for(int icell=0;icell<ncells;icell++){
        //updateArea(changedPosition[atom]);
        //changedCells.clear();
        //for(vector<int>::iterator iter=searchList.begin();
        //    iter!=searchList.end();iter+=3){
            
        //    icell=cellPos2Index(*iter, *iter+1, *iter+2);
        //    changedCells.push_back(icell);
        //cout<<changedPosition[atom][0]<<", "<<changedPosition[atom][1]<<", "<<changedPosition[atom][2]<<endl;
        // openmp maybe work here
        for(int atom=0;atom<natoms;atom++){
            if(ifOccupied(changedPosition[atom][0], 
                          changedPosition[atom][1], 
                          changedPosition[atom][2], 
                          xcell_center[icell], 
                          ycell_center[icell], 
                          zcell_center[icell])){
                    
                occupancy[icell]+=1;
            }
        }
    }
    //cout<< "update insertion call write occupancy"<<endl;
    //writeOccupancy("Occinsa",n);
    
    updateInsOccupancy();
    //getCavities();

}


/**
 * @brief update after attemption but before real deletion
 * 
 * @param natoms number of atoms
 * @param changedPosition the changed coordinates
 * 
*/
void Grid::updateDeletion(int natoms, double (*changedPosition)[3], int n){
    //first backup the information for restore
    tmp_ncavities=ncavities;
    tmp_cavities=cavities;
    tmp_occupancy=occupancy;
    //writeOccupancy("Occdelb",n);
    //get the search area
    //int icell;
    //changedCells.clear();
    #pragma omp parallel for
    for(int icell=0;icell<ncells;icell++){
        //cout<<changedPosition[atom][0]<<", "<<changedPosition[atom][1]<<", "<<changedPosition[atom][2]<<endl;
        //openmp maybe work here
        for(int atom=0;atom<natoms;atom++){
            if(  ifOccupied(changedPosition[atom][0],
                            changedPosition[atom][1], 
                            changedPosition[atom][2], 
                            xcell_center[icell], 
                            ycell_center[icell], 
                            zcell_center[icell]) ){

                    occupancy[icell]-=1;
                    if (occupancy[icell]<0){
                        cout<< "find negative occupancy"<<endl;
                        cout<<"coor: "<<changedPosition[atom][0]<<", "<<changedPosition[atom][1]<<", "<<changedPosition[atom][2]<<endl;
                        cout<<"cell: "<<icell<<": "<<xcell_center[icell]<<", "<<ycell_center[icell]<<","<<zcell_center[icell]<<endl;
                        exit(1);
                    }
            }
        //updateArea(changedPosition[atom]);
        //for(vector<int>::iterator iter=searchList.begin();
        //    iter!=searchList.end();iter+=3){

        //    icell=cellPos2Index(*iter, *(iter+1), *(iter+2));
        //    changedCells.push_back(icell);
        //cout<<"icell loop"<<endl;
        
            //changedCells.push_back(icell);
            //cout<<icell<<": "<< changedPosition[atom][0]<<", "<< changedPosition[atom][1]<<", "<< changedPosition[atom][2]<<", "<< xcell_center[icell]<<", "<< ycell_center[icell]<<", "<< zcell_center[icell]<<":";
            

                    
                //if(!occupancy[icell]==0){//!!!maybe problem why there will be negnative???
              //      cout<< "Occ"<<endl;
                //}
            //}else{
            //    cout<<endl;
        
        }
    }

    //writeOccupancy("Occdela",n);

    updateDelOccupancy();
    //getCavities();    
}

/**
 * @brief update the Occupancy after doing Displacement,
 *          same as the initialize
*/
void Grid::updateDisplacement(int natoms, vector<double> &xpos, 
                    vector<double> &ypos, vector<double> &zpos){
    
    occupancy=vector<int>(ncells,0);
    #pragma omp parallel for
    for(int icell=0;icell<ncells;icell++){
        //cout<<"x0: "<<xpos[i]<<"y0: "<<ypos[i]<<"z0: "<<zpos[i]<<endl;
        //openmp may work here
        for(int atom=0;atom<natoms;atom++){
            if(ifOccupied(xpos[atom], ypos[atom], zpos[atom],
                xcell_center[icell], ycell_center[icell], zcell_center[icell])){
                
                occupancy[icell]+=1;
            }
        }
    }
    getCavities();

}


/**
 * @brief deletion will only decrease the occupancy and transfer 
 *          none cavities into cavities
 * 
*/
void Grid::updateDelOccupancy(){
    //for(vector<int>::iterator iter=changedCells.begin();
    //    iter!=changedCells.end();iter++){
        //omp may work here
    #pragma omp parallel for
    for(int icell=0;icell<ncells;icell++){
        if(occupancy[icell]==0){
            //if(!cavities[icell]){
            cavities[icell]=1;
                //ncavities+=1;
                //}
        }
    }
}

/**
 * @brief insertion will only increase the occupancy and transfer 
 *          cavities into none cavities
 * 
*/
void Grid::updateInsOccupancy(){

    //for(vector<int>::iterator iter=changedCells.begin();
    //    iter!=changedCells.end();iter++){
        //omp may work here
    #pragma omp parallel for
    for(int icell=0;icell<ncells;icell++){
        if(occupancy[icell]>0){
            //if(cavities[icell]){
                cavities[icell]=0;
                //ncavities-=1;
            //}
        }
    }
}

/**
 * @brief update the cavityList
 * 
*/
void Grid::getCavities(){

    cavityList.clear();
    ncavities=0;
    
    for(int i=0;i!=ncells;i++){
        if(cavities[i]){
            ncavities+=1;
            cavityList.push_back(i);
        }
    }
}

/**
 *@brief Restore the grid information before deletion
*/
void Grid::restoreDeletion(){
    ncavities=tmp_ncavities;
    cavities=tmp_cavities;
    occupancy=tmp_occupancy;
    return;
    
}

/**
 *@brief return the probability to find a cavity
*/
double Grid::calProbcavities() const{
    return double(ncavities)/ncells;
}

/**
 * @brief Destructor
 *
*/
//rid::~Grid(){
//   cout<<"Grid destructor: "<<endl;
//   //delete [] cavities;
//   //delete [] tmp_cavities;
//   cout<<"finish Grid destructor: "<<endl;
//

/**
 * @brief only for debug
 */
void Grid::writeOccupancy(string fname,int n){
    ofstream fp;
    string filename=fname+to_string(n)+".dat";
    fp.open(filename,ios::out);
    for(int i=0;i<ncells;i++){
        fp<<occupancy[i]<<endl;
    }
    fp.close();
}

#if DEBUGGRID
int main(){
    Grid mygrid(-10.0, -5.0, -5.0, 10.0, 5.0, 5.0, 0.25, 0.25, 0.25, 2.5);
    mygrid.getCellBoundary(46);
    cout<< "Boundary test:"<<endl;
    for(int i=0;i!=3;i++){
        cout << "["<<mygrid.cellBoundary[i][0]<<", "<<mygrid.cellBoundary[i][1]<<"]"<<endl;
    }
    cout<< "cell Index test:" <<mygrid.cellPos2Index(16,13,14)<<endl;

}
#endif
