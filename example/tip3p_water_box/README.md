# Example1: TIP3P water box

This directory contains the files needs to perform GCMC/MD of a small water box. This example contains a water box of $20~\times~10~\times~10~\mathrm{\AA}^3$ with 72 water molecules as the initial states.
```bash
├── energy.inp: CHARMM input template  
├── example.inp: input file for mygcbgcmc  
├── initial.crd: CHARMM coordinate file  
├── initial.psf: CHARMM psf file  
└── toppar_water.str: CHARMM topology file (rtf and prm)
```
Among these files, the CHARMM-related files need no more detailed introductions. The program works by calling CHARMM to calculate the energy and force as well as propagate MD simulations. To achieve this a parameter named **jobtype** is passed to **energy.inp** to control the job type.


## The input file for mygcbgcmc
Next, let's go through the input file for **mygcbgcmc** line by line:
```bash
3      
3    
-4.98
        0.0303017300        0.1605113000       -0.3547918800
       -0.3183342600        0.5610223950        0.4416223650
        0.2880325440       -0.7215336850       -0.0868304790
TIP3 
OH2 
H1  
H2  
OT 
HT
HT
initial.crd
initial.psf
0.10026
energy.inp
charmm
     -10.0      -5.0      -5.0
      10.0       5.0       5.0
     0.250     0.250     0.250
2.5
298.15
1000
500

#do anything you want

```
---
The line 1,2 are the number of particles in the TIP3P water molecule.  
The line 3 is the Adams $B$ parameter for carrying out GCMC/MD.  
The line 4-6 are the coordinates of optimized TIP3P water molecules.  
The line 7 is the residue name of TIP3P water molecule in CHARMM.  
The line 8,9,10 are the atom name of TIP3P water molecule in CHARMM force field.  
The line 11,12,13 are the atom type of TIP3P water in CHARMM force field.  
The line 14,15 are the file name for CHARMM coordinate file and psf file of the simulating system.  
The line 16 is a float without any meaning.
The line 17 is the CHARMM input template file.  
The line 18 is the path to call CHARMM.  
The lines 19,20 defines the region where GCMC movements are happening.  
The line 21 defines the grid size for grid cavity bias.  
The line 22 defines the cut off to judge if a grid is occupied.
The line 23 defines the temperature for GCMC/MD.
The line 24 defines the total GCMC steps.
The line 25 defines MD simulation after how many steps GCMC movements and also to save a frame.

---
To run the GCMC/MD simulations:
```bash
mygcbgcmc.exe -i example.inp
```
which will output the information on the screen.
To redirect the output:
```bash
mygcbgcmc.exe -i example.inp > log
```
The output will be written into a log file.