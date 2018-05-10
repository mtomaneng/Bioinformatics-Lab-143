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



















