---
title: "Class15"
author: "Matthew Tomaneng"
date: "5/22/2018"
output: 
  html_document: 
    keep_md: yes
---


##Section 1
Use DESeq2 again


```r
#library(dplyr)
library(DESeq2)
```

```
## Loading required package: S4Vectors
```

```
## Loading required package: stats4
```

```
## Loading required package: BiocGenerics
```

```
## Loading required package: parallel
```

```
## 
## Attaching package: 'BiocGenerics'
```

```
## The following objects are masked from 'package:parallel':
## 
##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
##     clusterExport, clusterMap, parApply, parCapply, parLapply,
##     parLapplyLB, parRapply, parSapply, parSapplyLB
```

```
## The following objects are masked from 'package:stats':
## 
##     IQR, mad, sd, var, xtabs
```

```
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, cbind, colMeans,
##     colnames, colSums, do.call, duplicated, eval, evalq, Filter,
##     Find, get, grep, grepl, intersect, is.unsorted, lapply,
##     lengths, Map, mapply, match, mget, order, paste, pmax,
##     pmax.int, pmin, pmin.int, Position, rank, rbind, Reduce,
##     rowMeans, rownames, rowSums, sapply, setdiff, sort, table,
##     tapply, union, unique, unsplit, which, which.max, which.min
```

```
## 
## Attaching package: 'S4Vectors'
```

```
## The following object is masked from 'package:base':
## 
##     expand.grid
```

```
## Loading required package: IRanges
```

```
## Loading required package: GenomicRanges
```

```
## Loading required package: GenomeInfoDb
```

```
## Loading required package: SummarizedExperiment
```

```
## Loading required package: Biobase
```

```
## Welcome to Bioconductor
## 
##     Vignettes contain introductory material; view with
##     'browseVignettes()'. To cite Bioconductor, see
##     'citation("Biobase")', and for packages 'citation("pkgname")'.
```

```
## Loading required package: DelayedArray
```

```
## Loading required package: matrixStats
```

```
## 
## Attaching package: 'matrixStats'
```

```
## The following objects are masked from 'package:Biobase':
## 
##     anyMissing, rowMedians
```

```
## 
## Attaching package: 'DelayedArray'
```

```
## The following objects are masked from 'package:matrixStats':
## 
##     colMaxs, colMins, colRanges, rowMaxs, rowMins, rowRanges
```

```
## The following object is masked from 'package:base':
## 
##     apply
```

load data files

```r
metaFile <- "data/GSE37704_metadata.csv"
countFile <- "data/GSE37704_featurecounts.csv"

# Import metadata and take a peak
colData = read.csv(metaFile, row.names=1)
head(colData)
```

```
##               condition
## SRR493366 control_sirna
## SRR493367 control_sirna
## SRR493368 control_sirna
## SRR493369      hoxa1_kd
## SRR493370      hoxa1_kd
## SRR493371      hoxa1_kd
```


```r
# Import countdata
countData = read.csv(countFile, row.names=1)
head(countData)
```

```
##                 length SRR493366 SRR493367 SRR493368 SRR493369 SRR493370
## ENSG00000186092    918         0         0         0         0         0
## ENSG00000279928    718         0         0         0         0         0
## ENSG00000279457   1982        23        28        29        29        28
## ENSG00000278566    939         0         0         0         0         0
## ENSG00000273547    939         0         0         0         0         0
## ENSG00000187634   3214       124       123       205       207       212
##                 SRR493371
## ENSG00000186092         0
## ENSG00000279928         0
## ENSG00000279457        46
## ENSG00000278566         0
## ENSG00000273547         0
## ENSG00000187634       258
```

We want to get rid of the odd "Length" column

```r
# Note we need to remove the odd first $length col
countData <- as.matrix(countData[,-1])
head(countData)
```

```
##                 SRR493366 SRR493367 SRR493368 SRR493369 SRR493370
## ENSG00000186092         0         0         0         0         0
## ENSG00000279928         0         0         0         0         0
## ENSG00000279457        23        28        29        29        28
## ENSG00000278566         0         0         0         0         0
## ENSG00000273547         0         0         0         0         0
## ENSG00000187634       124       123       205       207       212
##                 SRR493371
## ENSG00000186092         0
## ENSG00000279928         0
## ENSG00000279457        46
## ENSG00000278566         0
## ENSG00000273547         0
## ENSG00000187634       258
```

Filter out 0's in data

```r
# Filter count data where you have 0 read count across all samples.
countData = countData[rowSums(countData)>1, ]
head(countData)
```

```
##                 SRR493366 SRR493367 SRR493368 SRR493369 SRR493370
## ENSG00000279457        23        28        29        29        28
## ENSG00000187634       124       123       205       207       212
## ENSG00000188976      1637      1831      2383      1226      1326
## ENSG00000187961       120       153       180       236       255
## ENSG00000187583        24        48        65        44        48
## ENSG00000187642         4         9        16        14        16
##                 SRR493371
## ENSG00000279457        46
## ENSG00000187634       258
## ENSG00000188976      1504
## ENSG00000187961       357
## ENSG00000187583        64
## ENSG00000187642        16
```

Now our data is usable for DESeq

```r
dds = DESeqDataSetFromMatrix(countData=countData,
                             colData=colData,
                             design=~condition)
dds = DESeq(dds)
```

```
## estimating size factors
```

```
## estimating dispersions
```

```
## gene-wise dispersion estimates
```

```
## mean-dispersion relationship
```

```
## final dispersion estimates
```

```
## fitting model and testing
```

Lets checkk out dds

```r
dds
```

```
## class: DESeqDataSet 
## dim: 15280 6 
## metadata(1): version
## assays(3): counts mu cooks
## rownames(15280): ENSG00000279457 ENSG00000187634 ...
##   ENSG00000276345 ENSG00000271254
## rowData names(21): baseMean baseVar ... deviance maxCooks
## colnames(6): SRR493366 SRR493367 ... SRR493370 SRR493371
## colData names(2): condition sizeFactor
```

#We are able to work with specific genes now
We will look at “hoxa1_kd” and “control_sirna”

```r
resultsNames(dds)
```

```
## [1] "Intercept"                          
## [2] "condition_hoxa1_kd_vs_control_sirna"
```

```r
res = results(dds, contrast=c("condition", "hoxa1_kd", "control_sirna"))
```

Reorder p-values and call a summary to get a sense of which is up and down regulated

```r
res = res[order(res$pvalue),]
summary(res)
```

```
## 
## out of 15280 with nonzero total read count
## adjusted p-value < 0.1
## LFC > 0 (up)     : 4352, 28% 
## LFC < 0 (down)   : 4400, 29% 
## outliers [1]     : 0, 0% 
## low counts [2]   : 590, 3.9% 
## (mean count < 1)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
```

Our data is labeled with Ensembl notation but we need Entrez Gene ID's

```r
#This is calling out the data tables with the translation names
library("AnnotationDbi")
library("org.Hs.eg.db")
```

```
## 
```

```r
columns(org.Hs.eg.db)
```

```
##  [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT" 
##  [5] "ENSEMBLTRANS" "ENTREZID"     "ENZYME"       "EVIDENCE"    
##  [9] "EVIDENCEALL"  "GENENAME"     "GO"           "GOALL"       
## [13] "IPI"          "MAP"          "OMIM"         "ONTOLOGY"    
## [17] "ONTOLOGYALL"  "PATH"         "PFAM"         "PMID"        
## [21] "PROSITE"      "REFSEQ"       "SYMBOL"       "UCSCKG"      
## [25] "UNIGENE"      "UNIPROT"
```

We can use this to translate our Ensebl to Entrez

```r
#Transform symbol into the Ensembl ID
res$symbol = mapIds(org.Hs.eg.db,
                    keys=row.names(res), 
                    column="SYMBOL",
                    keytype="ENSEMBL",
                    multiVals="first")
```

```
## 'select()' returned 1:many mapping between keys and columns
```

```r
#Transform EntrezID to Ensembl
res$entrez = mapIds(org.Hs.eg.db,
                    keys=row.names(res), 
                    column="ENTREZID",
                    keytype="ENSEMBL",
                    multiVals="first")
```

```
## 'select()' returned 1:many mapping between keys and columns
```

```r
#Transform the genename to Ensembl
res$name =   mapIds(org.Hs.eg.db,
                    keys=row.names(res), 
                    column="GENENAME",
                    keytype="ENSEMBL",
                    multiVals="first")
```

```
## 'select()' returned 1:many mapping between keys and columns
```

```r
head(res, 10)
```

```
## log2 fold change (MLE): condition hoxa1_kd vs control_sirna 
## Wald test p-value: condition hoxa1 kd vs control sirna 
## DataFrame with 10 rows and 9 columns
##                  baseMean log2FoldChange      lfcSE      stat    pvalue
##                 <numeric>      <numeric>  <numeric> <numeric> <numeric>
## ENSG00000117519  4483.627      -2.422719 0.06001850 -40.36620         0
## ENSG00000183508  2053.881       3.201955 0.07241968  44.21388         0
## ENSG00000159176  5692.463      -2.313737 0.05757255 -40.18820         0
## ENSG00000150938  7442.986      -2.059631 0.05386627 -38.23601         0
## ENSG00000116016  4423.947      -1.888019 0.04318301 -43.72134         0
## ENSG00000136068  3796.127      -1.649792 0.04394825 -37.53942         0
## ENSG00000164251  2348.770       3.344508 0.06907610  48.41773         0
## ENSG00000124766  2576.653       2.392288 0.06171493  38.76352         0
## ENSG00000124762 28106.119       1.832258 0.03892405  47.07264         0
## ENSG00000106366 43719.126      -1.844046 0.04194432 -43.96415         0
##                      padj      symbol      entrez
##                 <numeric> <character> <character>
## ENSG00000117519         0        CNN3        1266
## ENSG00000183508         0      FAM46C       54855
## ENSG00000159176         0       CSRP1        1465
## ENSG00000150938         0       CRIM1       51232
## ENSG00000116016         0       EPAS1        2034
## ENSG00000136068         0        FLNB        2317
## ENSG00000164251         0       F2RL1        2150
## ENSG00000124766         0        SOX4        6659
## ENSG00000124762         0      CDKN1A        1026
## ENSG00000106366         0    SERPINE1        5054
##                                                        name
##                                                 <character>
## ENSG00000117519                                  calponin 3
## ENSG00000183508 family with sequence similarity 46 member C
## ENSG00000159176         cysteine and glycine rich protein 1
## ENSG00000150938 cysteine rich transmembrane BMP regulator 1
## ENSG00000116016            endothelial PAS domain protein 1
## ENSG00000136068                                   filamin B
## ENSG00000164251                 F2R like trypsin receptor 1
## ENSG00000124766                                   SRY-box 4
## ENSG00000124762        cyclin dependent kinase inhibitor 1A
## ENSG00000106366                    serpin family E member 1
```

##Section 2: Pathway Analysis
Download packages we need

```r
#source("http://bioconductor.org/biocLite.R")
#biocLite( c("pathview", "gage", "gageData") )
```


```r
library(pathview)
```

```
## ##############################################################################
## Pathview is an open source software package distributed under GNU General
## Public License version 3 (GPLv3). Details of GPLv3 is available at
## http://www.gnu.org/licenses/gpl-3.0.html. Particullary, users are required to
## formally cite the original Pathview paper (not just mention it) in publications
## or products. For details, do citation("pathview") within R.
## 
## The pathview downloads and uses KEGG data. Non-academic uses may require a KEGG
## license agreement (details at http://www.kegg.jp/kegg/legal.html).
## ##############################################################################
```


```r
library(gage)
library(gageData)

data(kegg.sets.hs)
data(sigmet.idx.hs)

kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]
head(kegg.sets.hs, 3)
```

```
## $`hsa00232 Caffeine metabolism`
## [1] "10"   "1544" "1548" "1549" "1553" "7498" "9"   
## 
## $`hsa00983 Drug metabolism - other enzymes`
##  [1] "10"     "1066"   "10720"  "10941"  "151531" "1548"   "1549"  
##  [8] "1551"   "1553"   "1576"   "1577"   "1806"   "1807"   "1890"  
## [15] "221223" "2990"   "3251"   "3614"   "3615"   "3704"   "51733" 
## [22] "54490"  "54575"  "54576"  "54577"  "54578"  "54579"  "54600" 
## [29] "54657"  "54658"  "54659"  "54963"  "574537" "64816"  "7083"  
## [36] "7084"   "7172"   "7363"   "7364"   "7365"   "7366"   "7367"  
## [43] "7371"   "7372"   "7378"   "7498"   "79799"  "83549"  "8824"  
## [50] "8833"   "9"      "978"   
## 
## $`hsa00230 Purine metabolism`
##   [1] "100"    "10201"  "10606"  "10621"  "10622"  "10623"  "107"   
##   [8] "10714"  "108"    "10846"  "109"    "111"    "11128"  "11164" 
##  [15] "112"    "113"    "114"    "115"    "122481" "122622" "124583"
##  [22] "132"    "158"    "159"    "1633"   "171568" "1716"   "196883"
##  [29] "203"    "204"    "205"    "221823" "2272"   "22978"  "23649" 
##  [36] "246721" "25885"  "2618"   "26289"  "270"    "271"    "27115" 
##  [43] "272"    "2766"   "2977"   "2982"   "2983"   "2984"   "2986"  
##  [50] "2987"   "29922"  "3000"   "30833"  "30834"  "318"    "3251"  
##  [57] "353"    "3614"   "3615"   "3704"   "377841" "471"    "4830"  
##  [64] "4831"   "4832"   "4833"   "4860"   "4881"   "4882"   "4907"  
##  [71] "50484"  "50940"  "51082"  "51251"  "51292"  "5136"   "5137"  
##  [78] "5138"   "5139"   "5140"   "5141"   "5142"   "5143"   "5144"  
##  [85] "5145"   "5146"   "5147"   "5148"   "5149"   "5150"   "5151"  
##  [92] "5152"   "5153"   "5158"   "5167"   "5169"   "51728"  "5198"  
##  [99] "5236"   "5313"   "5315"   "53343"  "54107"  "5422"   "5424"  
## [106] "5425"   "5426"   "5427"   "5430"   "5431"   "5432"   "5433"  
## [113] "5434"   "5435"   "5436"   "5437"   "5438"   "5439"   "5440"  
## [120] "5441"   "5471"   "548644" "55276"  "5557"   "5558"   "55703" 
## [127] "55811"  "55821"  "5631"   "5634"   "56655"  "56953"  "56985" 
## [134] "57804"  "58497"  "6240"   "6241"   "64425"  "646625" "654364"
## [141] "661"    "7498"   "8382"   "84172"  "84265"  "84284"  "84618" 
## [148] "8622"   "8654"   "87178"  "8833"   "9060"   "9061"   "93034" 
## [155] "953"    "9533"   "954"    "955"    "956"    "957"    "9583"  
## [162] "9615"
```

gage() function requires a named vector of fold changes, where the names of the values are the Entrez gene IDs

```r
foldchanges = res$log2FoldChange
names(foldchanges) = res$entrez
head(foldchanges)
```

```
##      1266     54855      1465     51232      2034      2317 
## -2.422719  3.201955 -2.313737 -2.059631 -1.888019 -1.649792
```

Lets use Gage function. This will give us overlap in our data and other experimental data

```r
# Get the results
keggres = gage(foldchanges, gsets=kegg.sets.hs, same.dir=TRUE)
```

Lets look at the result object. It is a list with three elements (“greater”, “less” and “stats”).

```r
attributes(keggres)
```

```
## $names
## [1] "greater" "less"    "stats"
```

So it is a list object (you can check it with str(keggres)) and we can use the dollar syntax to access a named element, e.g.

```r
#Upregulated
head(keggres$greater)
```

```
##                                         p.geomean stat.mean       p.val
## hsa04640 Hematopoietic cell lineage   0.002709366  2.857393 0.002709366
## hsa04630 Jak-STAT signaling pathway   0.005655916  2.557207 0.005655916
## hsa04142 Lysosome                     0.008948808  2.384783 0.008948808
## hsa00140 Steroid hormone biosynthesis 0.009619717  2.432105 0.009619717
## hsa04740 Olfactory transduction       0.014450242  2.239717 0.014450242
## hsa04916 Melanogenesis                0.022339115  2.023074 0.022339115
##                                           q.val set.size        exp1
## hsa04640 Hematopoietic cell lineage   0.3847887       49 0.002709366
## hsa04630 Jak-STAT signaling pathway   0.3847887      103 0.005655916
## hsa04142 Lysosome                     0.3847887      117 0.008948808
## hsa00140 Steroid hormone biosynthesis 0.3847887       26 0.009619717
## hsa04740 Olfactory transduction       0.4624078       39 0.014450242
## hsa04916 Melanogenesis                0.5297970       85 0.022339115
```


```r
#Downregulated
head(keggres$less)
```

```
##                                      p.geomean stat.mean        p.val
## hsa04110 Cell cycle               1.004024e-05 -4.353447 1.004024e-05
## hsa03030 DNA replication          8.909718e-05 -3.968605 8.909718e-05
## hsa03013 RNA transport            1.471026e-03 -3.007785 1.471026e-03
## hsa04114 Oocyte meiosis           1.987557e-03 -2.915377 1.987557e-03
## hsa03440 Homologous recombination 2.942017e-03 -2.868137 2.942017e-03
## hsa00240 Pyrimidine metabolism    5.800212e-03 -2.549616 5.800212e-03
##                                         q.val set.size         exp1
## hsa04110 Cell cycle               0.001606438      120 1.004024e-05
## hsa03030 DNA replication          0.007127774       36 8.909718e-05
## hsa03013 RNA transport            0.078454709      143 1.471026e-03
## hsa04114 Oocyte meiosis           0.079502292       98 1.987557e-03
## hsa03440 Homologous recombination 0.094144560       28 2.942017e-03
## hsa00240 Pyrimidine metabolism    0.138500584       95 5.800212e-03
```

Combine the two data tables view together

```r
lapply(keggres, head)
```

```
## $greater
##                                         p.geomean stat.mean       p.val
## hsa04640 Hematopoietic cell lineage   0.002709366  2.857393 0.002709366
## hsa04630 Jak-STAT signaling pathway   0.005655916  2.557207 0.005655916
## hsa04142 Lysosome                     0.008948808  2.384783 0.008948808
## hsa00140 Steroid hormone biosynthesis 0.009619717  2.432105 0.009619717
## hsa04740 Olfactory transduction       0.014450242  2.239717 0.014450242
## hsa04916 Melanogenesis                0.022339115  2.023074 0.022339115
##                                           q.val set.size        exp1
## hsa04640 Hematopoietic cell lineage   0.3847887       49 0.002709366
## hsa04630 Jak-STAT signaling pathway   0.3847887      103 0.005655916
## hsa04142 Lysosome                     0.3847887      117 0.008948808
## hsa00140 Steroid hormone biosynthesis 0.3847887       26 0.009619717
## hsa04740 Olfactory transduction       0.4624078       39 0.014450242
## hsa04916 Melanogenesis                0.5297970       85 0.022339115
## 
## $less
##                                      p.geomean stat.mean        p.val
## hsa04110 Cell cycle               1.004024e-05 -4.353447 1.004024e-05
## hsa03030 DNA replication          8.909718e-05 -3.968605 8.909718e-05
## hsa03013 RNA transport            1.471026e-03 -3.007785 1.471026e-03
## hsa04114 Oocyte meiosis           1.987557e-03 -2.915377 1.987557e-03
## hsa03440 Homologous recombination 2.942017e-03 -2.868137 2.942017e-03
## hsa00240 Pyrimidine metabolism    5.800212e-03 -2.549616 5.800212e-03
##                                         q.val set.size         exp1
## hsa04110 Cell cycle               0.001606438      120 1.004024e-05
## hsa03030 DNA replication          0.007127774       36 8.909718e-05
## hsa03013 RNA transport            0.078454709      143 1.471026e-03
## hsa04114 Oocyte meiosis           0.079502292       98 1.987557e-03
## hsa03440 Homologous recombination 0.094144560       28 2.942017e-03
## hsa00240 Pyrimidine metabolism    0.138500584       95 5.800212e-03
## 
## $stats
##                                       stat.mean     exp1
## hsa04640 Hematopoietic cell lineage    2.857393 2.857393
## hsa04630 Jak-STAT signaling pathway    2.557207 2.557207
## hsa04142 Lysosome                      2.384783 2.384783
## hsa00140 Steroid hormone biosynthesis  2.432105 2.432105
## hsa04740 Olfactory transduction        2.239717 2.239717
## hsa04916 Melanogenesis                 2.023074 2.023074
```


```r
## Sanity check displaying all pathways data
pathways = data.frame(id=rownames(keggres$greater), keggres$greater)
head(pathways)
```

```
##                                                                          id
## hsa04640 Hematopoietic cell lineage     hsa04640 Hematopoietic cell lineage
## hsa04630 Jak-STAT signaling pathway     hsa04630 Jak-STAT signaling pathway
## hsa04142 Lysosome                                         hsa04142 Lysosome
## hsa00140 Steroid hormone biosynthesis hsa00140 Steroid hormone biosynthesis
## hsa04740 Olfactory transduction             hsa04740 Olfactory transduction
## hsa04916 Melanogenesis                               hsa04916 Melanogenesis
##                                         p.geomean stat.mean       p.val
## hsa04640 Hematopoietic cell lineage   0.002709366  2.857393 0.002709366
## hsa04630 Jak-STAT signaling pathway   0.005655916  2.557207 0.005655916
## hsa04142 Lysosome                     0.008948808  2.384783 0.008948808
## hsa00140 Steroid hormone biosynthesis 0.009619717  2.432105 0.009619717
## hsa04740 Olfactory transduction       0.014450242  2.239717 0.014450242
## hsa04916 Melanogenesis                0.022339115  2.023074 0.022339115
##                                           q.val set.size        exp1
## hsa04640 Hematopoietic cell lineage   0.3847887       49 0.002709366
## hsa04630 Jak-STAT signaling pathway   0.3847887      103 0.005655916
## hsa04142 Lysosome                     0.3847887      117 0.008948808
## hsa00140 Steroid hormone biosynthesis 0.3847887       26 0.009619717
## hsa04740 Olfactory transduction       0.4624078       39 0.014450242
## hsa04916 Melanogenesis                0.5297970       85 0.022339115
```

Let's try to use the pathview function

```r
pathview(gene.data=foldchanges, pathway.id="hsa04110")
```

```
## 'select()' returned 1:1 mapping between keys and columns
```

```
## Info: Working in directory /Users/matthewtomaneng/Documents/Bioinformatics BIMM 143/bimm143/Class15
```

```
## Info: Writing image file hsa04110.pathview.png
```

View the pathway

```r
# A different PDF based output of the same data
pathview(gene.data=foldchanges, pathway.id="hsa04110", kegg.native=FALSE)
```

```
## 'select()' returned 1:1 mapping between keys and columns
```

```
## Info: Working in directory /Users/matthewtomaneng/Documents/Bioinformatics BIMM 143/bimm143/Class15
```

```
## Info: Writing image file hsa04110.pathview.pdf
```


![](hsa04110.pathview.png)


Extracting just a small snippet of the important data

```r
## Focus on top 5 upregulated pathways here for demo purposes only
keggrespathways <- rownames(keggres$greater)[1:5]

# Extract the IDs part of each string
keggresids = substr(keggrespathways, start=1, stop=8)
keggresids
```

```
## [1] "hsa04640" "hsa04630" "hsa04142" "hsa00140" "hsa04740"
```

Finally, lets pass these IDs in keggresids to the pathview() function to draw plots for all the top 5 pathways.



















