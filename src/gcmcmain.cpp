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

#include <iostream>
#include <string>
#include "gcmccall.hpp"

using namespace std;

/**
 *@brief the main program of mg grid cavity biased GCMC
 * use -i to declare the input parameter file
 */
int main(int argc, char *argv[]){
    
    string inpname;
    //if the arguement is right
    if(argc!=3||argv[1]!=string("-i")){
        cout<<"Usage:"<<endl;
        cout<<"mygcbgcmc.exe -i inputfile"<<endl;
        exit(1);
    }
    //get the inputfile name and call 
    inpname=argv[2];

    //read parameter and run
    gcmcio(inpname);

}
