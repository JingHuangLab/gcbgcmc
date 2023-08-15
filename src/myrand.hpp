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


#ifndef CPPMYGCBGCMCMYRAND_H //cpp my gcb gcmc
#define CPPMYGCBGCMCMYRAND_H


/**
 * @brief:generate random number for the gcmc
 * 
 * 
*/
#include<iostream>
#include<string>
#include<cstdlib>
#include<ctime>
#include<vector>
#include"MersenneTwister.h"
using namespace std;



/**
 * @brief Return a random number between 0 and 1 using the MTRand 
 * @param mtrand MTRand class
 * @return a random number between 0 and 1
 * 
*/
inline double myrand(MTRand &mtrand){return mtrand.rand(); }

/**
 * @brief Return a value from the give vector
 * 
 * @param candidateList a vector
 * 
 * @return a int value in the give vector
*/
inline int choice(vector<int> &candidateList, MTRand &mtrand){
    // at least one gcmcmole so it cannot be empty
    if(candidateList.empty()){
        cout<< "rand: choice is empty!!!"<<endl;
        exit(1);
    }

    //int vectorsize=candidateList.size()-1;
    //cout<<"in choice"<<endl;
    //cout<< vectorsize<<endl;
    int rand_resid=mtrand.randInt(candidateList.size()-1);
    //cout<< rand_resid<<endl;
    //int chose=candidateList[rand_resid];
    //cout<<chose <<endl;
    return candidateList[rand_resid];
}

/**
 * @brief random number among 0, 1, 2
 * 
 * @return a random int value among 0, 1, 2 
*/
inline int randmove(MTRand &mtrand){return mtrand.randInt(2);}

#if DEBUG
int main(){
    cout<< myrand()   <<endl;
    cout<< randmove()   <<endl;

    vector<int> alist;
    for(int i=15;i!=50;i++){
        alist.push_back(i);
    }
    cout << choice(alist) << endl;
    return 0;
}
#endif

#endif