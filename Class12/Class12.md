---
title: "Class12"
author: "Matthew Tomaneng"
date: "5/10/2018"
output: 
  html_document: 
    keep_md: yes
---


Section 1

```r
#Load bio3d and download protein we're working with
#Already have it downloaded so just put it in the project enviornment
library(bio3d)
```
Next we will tranlate the file into something R can read

```r
hiv <- read.pdb("1hsg.pdb")
hiv
```

```
## 
##  Call:  read.pdb(file = "1hsg.pdb")
## 
##    Total Models#: 1
##      Total Atoms#: 1686,  XYZs#: 5058  Chains#: 2  (values: A B)
## 
##      Protein Atoms#: 1514  (residues/Calpha atoms#: 198)
##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
## 
##      Non-protein/nucleic Atoms#: 172  (residues: 128)
##      Non-protein/nucleic resid values: [ HOH (127), MK1 (1) ]
## 
##    Protein sequence:
##       PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYD
##       QILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNFPQITLWQRPLVTIKIGGQLKE
##       ALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTP
##       VNIIGRNLLTQIGCTLNF
## 
## + attr: atom, xyz, seqres, helix, sheet,
##         calpha, remark, call
```
#Q1: What are the two non protein components
#A: Water and MK1 (ligand)

Section 1.2: Prepare initial protein and ligand input files

```r
#Trim data to components we can use isolated from each other
prot <- trim.pdb(hiv, "protein")
lig <- trim.pdb(hiv, "ligand")
#Now we are going to make these their own PDB files
write.pdb(prot, file="1hsg_protein.pdb")
write.pdb(lig, "1hsg_ligand.pdb")
```

Section 1.3: Using AutoDockTools to setup protein docking input

```r
#Most PDB files do not have charge and atom type so we need to add it using ADT
```

#Q2: Can you see the active site
#A: No, too convoluted with the structure and water atoms

```r
#Used ADT to add charges and hydrogen atoms to the PDB file and is saved as a pdbqt file
#We then made a grid box enclosing the active site
#Centered at (16,25,4) with dimmensions (30,30,30)
```

Section: 1.4 Prepare the ligand

```r
#Now doing the same thing done above but with the ligand file
```

Section: 1.5 Prepare a docking configuration file

```r
#Made a fonfig file for the dimmensions of our docking
#It will give the computer the info it needs to dock where we want it to
```

Section 2 Docking ligands into HIV-1 protease

```r
#Downloaded autdock vita to run the docking simulation
```
Section 2.2 Docking indinavir into HIV-1 protease

```r
# ~/Downloads/autodock_vina_1_1_2_mac/bin/vina --config config.txt --log
#log.tx
##This will run our docking simulation with the config file we made earlier
```

Section 2.3 Visualizing Results in R

```r
library(bio3d)
res <- read.pdb("all.pdbqt", multi=TRUE)
write.pdb(res, "results.pdb")
#Go back to VMD to load the results
```

#Q4: Qualitatively how does the generated look to the actual
#A: Pretty darn close

```r
#Quantitative Analysis of the results
# res <- read.pdb("all.pdbqt", multi=TRUE)
ori <- read.pdb("1hsg_ligand.pdbqt")
rmsd(ori, res)
```

```
##  [1]  0.697  4.195 11.146 10.606 10.852 10.945 10.945  3.844  5.473  4.092
## [11] 10.404  5.574  3.448 11.396  6.126  3.848  8.237 11.196 10.981 11.950
```
#Q5: Quantitatively how good are the docks? Is the crystal binding mode reproduced within 1Å, RMSD for all atoms?
#A:Results are way above 1Angstrom so the results must not be too close

#Q6How would you determine the RMSD for heavy atoms only (i.e. non hydrogen atoms)?
#HINT: The atom.select() function will be of help here along with the selection string “noh” for no
#hydrogens.

```r
inds.res <- atom.select(res, "noh")
inds.res
```

```
## 
##  Call:  atom.select.pdb(pdb = res, string = "noh")
## 
##    Atom Indices#: 45  ($atom)
##    XYZ  Indices#: 135  ($xyz)
## 
## + attr: atom, xyz, call
```

```r
res$xyz[inds.res$xyz]
```

```
##   [1]  9.095  9.296 16.405 16.058 17.351 15.648 15.179 10.478 12.415  6.689
##  [11] 14.656 10.570  6.627 15.090 11.015 10.441 15.020 18.286 19.395 11.649
##  [21] 24.118 25.371 21.136 21.569 21.285 22.993 23.213 24.792 23.917 25.043
##  [31] 19.489 26.015 24.937 23.291 25.416 25.529 22.640 17.662 15.555  9.030
##  [41]  8.894  8.836  6.844  6.945  7.751  4.235  1.931  3.741  8.694  8.608
##  [51]  7.915  5.566  6.425  3.708  9.371  9.361 15.902 15.707 17.005 15.614
##  [61] 15.319 10.068 11.933  6.832 15.041 10.487  7.243 15.277 10.858 18.712
##  [71] 11.264 24.263 25.307 20.473 20.739 20.992 21.640 22.225 24.601 24.627
##  [81] 25.847 22.245 24.643 24.741 21.690 17.635 16.175 23.978  4.540  4.184
##  [91]  7.797  7.711  7.420  6.226  5.841  6.334  5.450  3.174  5.105  6.731
## [101]  4.468  6.062  7.432  7.349  6.898  7.006  5.257  4.745 10.553 10.602
## [111] 16.060 15.839 15.949 15.192  8.192 14.349 11.025  8.726 14.552 11.753
## [121] 11.326 16.413 18.005 18.978  9.836 25.068 25.740 21.119 21.280 21.792
## [131] 21.443 22.368 24.090 24.934 26.304
```

