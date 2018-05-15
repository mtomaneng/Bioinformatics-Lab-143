---
title: "Class13"
author: "Matthew Tomaneng"
date: "5/15/2018"
output: 
  html_document: 
    keep_md: yes
---



## Genome Informatics (Part 1)
#Section 1
Q1: What are those 4 candidate SNPs? 

Q3: What is the location of rs8067378 and what are the different alleles for rs8067378? 
A:Chromosome 17: 39895095 (forward strand) A/G variants (43% G) 

Q4: What are the downstream genes for rs8067378? Any genes named ZPBP2,
GSDMB, and ORMDL3? 
A: GSDMB and ORMDL3 are downstream

Q5: What proportion of the Mexican Ancestry in Los Angeles sample population (MXL)
are homozygous for the asthma associated SNP (G|G)?
[HINT: You can download a CVS file for this population from ENSEMBLE and use the R
functions read.csv(), and table() to answer this question] 

```r
mxl <- read.table("Mexican_genome.csv", header = TRUE)
table(mxl$X..Genotype)
```

```
## 
## ,"A|A","ALL, ,"A|G","ALL, ,"G|A","ALL, ,"G|G","ALL, 
##           22           21           12            9
```

```r
#table gives us a summary of the different outputs for a column
nrow(mxl)
```

```
## [1] 64
```

```r
#number of total genomes present
9/64
```

```
## [1] 0.140625
```

```r
# 14% is the answer
```

Q6. Back on the ENSEMBLE page, search for the particular sample HG00109. This is a
male from the GBR population group. What is the genotype for this sample?
A: G|G 


#Section 2: Initial RNA-Seq analysis 
Look more closely at sample HG00109 with G|G genotype for this SNP (from GBR population)
Q7: How many sequences are there in the first file?
What is the file size and format of the data? Make sure
the format is fastqsanger here!
A: 3863 sequences

Q8:  Does the first sequence have good quality? 

```r
#install.packages("seqinr")
#install.packages("gtools")
library(seqinr)
library(gtools)
phred <- asc(s2c("DDDDCDEDCDDDDBBDDDCC@")) - 33
phred
```

```
##  D  D  D  D  C  D  E  D  C  D  D  D  D  B  B  D  D  D  C  C  @ 
## 35 35 35 35 34 35 36 35 34 35 35 35 35 33 33 35 35 35 34 34 31
```
A: Yes decent quality


```r
#http://129.114.104.32/u/spongebob/h/class-13-matt
##This is a link to the history of using the supercomputer "Galaxy"
```


##Section 4 Population Analysis

```r
#Take file and convert into a file R can use
pop <- read.csv("population_asthma.txt")
summary(pop)
```

```
##                  sample.geno.exp
##  1 HG00367 A/G 28.96038  :  1   
##  10 HG00115 A/G 33.85374 :  1   
##  100 NA19189 A/G 30.63079:  1   
##  101 HG00155 A/G 19.1042 :  1   
##  102 HG00111 A/A 40.82922:  1   
##  103 NA12827 A/G 25.70962:  1   
##  (Other)                 :456
```


```r
table(pop$geno)
```

```
## < table of extent 0 >
```

```r
nrow(pop)
```

```
## [1] 462
```

```r
#LEt's plot this
```


```r
summary( pop$exp[ (pop$geno =="G/G")])
```

```
## Length  Class   Mode 
##      0   NULL   NULL
```

```r
 summary( pop$exp[ (pop$geno =="A/A")])
```

```
## Length  Class   Mode 
##      0   NULL   NULL
```

```r
 summary( pop$exp[ (pop$geno =="A/G")])
```

```
## Length  Class   Mode 
##      0   NULL   NULL
```



```r
#boxplot(pop ~ geno, data = pop, notch = TRUE)
```



```r
#install.packages("ggplot2")
library(ggplot2)
#ggplot(pop, aes(exp, fill=geno)) + geom_boxplot()
```










