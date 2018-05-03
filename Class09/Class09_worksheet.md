---
title: "Class09"
author: "Matthew Tomaneng"
date: "5/1/2018"
output:
  html_document:
    keep_md: yes 
---




```r
wisc.df <- read.csv("~/Downloads/WisconsinCancer.csv")
```


#converting features of wisc.df

```r
wisc.data <- as.matrix(wisc.df[,3:32])
row.names(wisc.data) <- wisc.df$id
head(wisc.data)
```

```
##          radius_mean texture_mean perimeter_mean area_mean smoothness_mean
## 842302         17.99        10.38         122.80    1001.0         0.11840
## 842517         20.57        17.77         132.90    1326.0         0.08474
## 84300903       19.69        21.25         130.00    1203.0         0.10960
## 84348301       11.42        20.38          77.58     386.1         0.14250
## 84358402       20.29        14.34         135.10    1297.0         0.10030
## 843786         12.45        15.70          82.57     477.1         0.12780
##          compactness_mean concavity_mean concave.points_mean symmetry_mean
## 842302            0.27760         0.3001             0.14710        0.2419
## 842517            0.07864         0.0869             0.07017        0.1812
## 84300903          0.15990         0.1974             0.12790        0.2069
## 84348301          0.28390         0.2414             0.10520        0.2597
## 84358402          0.13280         0.1980             0.10430        0.1809
## 843786            0.17000         0.1578             0.08089        0.2087
##          fractal_dimension_mean radius_se texture_se perimeter_se area_se
## 842302                  0.07871    1.0950     0.9053        8.589  153.40
## 842517                  0.05667    0.5435     0.7339        3.398   74.08
## 84300903                0.05999    0.7456     0.7869        4.585   94.03
## 84348301                0.09744    0.4956     1.1560        3.445   27.23
## 84358402                0.05883    0.7572     0.7813        5.438   94.44
## 843786                  0.07613    0.3345     0.8902        2.217   27.19
##          smoothness_se compactness_se concavity_se concave.points_se
## 842302        0.006399        0.04904      0.05373           0.01587
## 842517        0.005225        0.01308      0.01860           0.01340
## 84300903      0.006150        0.04006      0.03832           0.02058
## 84348301      0.009110        0.07458      0.05661           0.01867
## 84358402      0.011490        0.02461      0.05688           0.01885
## 843786        0.007510        0.03345      0.03672           0.01137
##          symmetry_se fractal_dimension_se radius_worst texture_worst
## 842302       0.03003             0.006193        25.38         17.33
## 842517       0.01389             0.003532        24.99         23.41
## 84300903     0.02250             0.004571        23.57         25.53
## 84348301     0.05963             0.009208        14.91         26.50
## 84358402     0.01756             0.005115        22.54         16.67
## 843786       0.02165             0.005082        15.47         23.75
##          perimeter_worst area_worst smoothness_worst compactness_worst
## 842302            184.60     2019.0           0.1622            0.6656
## 842517            158.80     1956.0           0.1238            0.1866
## 84300903          152.50     1709.0           0.1444            0.4245
## 84348301           98.87      567.7           0.2098            0.8663
## 84358402          152.20     1575.0           0.1374            0.2050
## 843786            103.40      741.6           0.1791            0.5249
##          concavity_worst concave.points_worst symmetry_worst
## 842302            0.7119               0.2654         0.4601
## 842517            0.2416               0.1860         0.2750
## 84300903          0.4504               0.2430         0.3613
## 84348301          0.6869               0.2575         0.6638
## 84358402          0.4000               0.1625         0.2364
## 843786            0.5355               0.1741         0.3985
##          fractal_dimension_worst
## 842302                   0.11890
## 842517                   0.08902
## 84300903                 0.08758
## 84348301                 0.17300
## 84358402                 0.07678
## 843786                   0.12440
```
# How many M's do we have

```r
table(wisc.df$diagnosis)
```

```
## 
##   B   M 
## 357 212
```
#Making a new vector called diagnosis

```r
diagnosis <- as.numeric(wisc.df$diagnosis == "M")
table(diagnosis)
```

```
## diagnosis
##   0   1 
## 357 212
```
#Should match previous table



## Exploratory data analysis
Q1. How many observations are in this dataset?

```r
nrow(wisc.df)
```

```
## [1] 569
```

Q2. How many variables/features in the data are suffixed with _mean?

```r
grep("mean", colnames(wisc.data))
```

```
##  [1]  1  2  3  4  5  6  7  8  9 10
```

```r
#produces a drop down of the differnt variables with _mean in the name
#12 is how many "_mean" columns we have
```

Q3. How many of the observations have a malignant diagnosis?

```r
sum(diagnosis)
```

```
## [1] 212
```

```r
#212 Malignant diagnosis 
```

##Section 2
#Check if rescaling is needed

```r
# Check column means and standard deviations
colMeans(wisc.data)
```

```
##             radius_mean            texture_mean          perimeter_mean 
##            1.412729e+01            1.928965e+01            9.196903e+01 
##               area_mean         smoothness_mean        compactness_mean 
##            6.548891e+02            9.636028e-02            1.043410e-01 
##          concavity_mean     concave.points_mean           symmetry_mean 
##            8.879932e-02            4.891915e-02            1.811619e-01 
##  fractal_dimension_mean               radius_se              texture_se 
##            6.279761e-02            4.051721e-01            1.216853e+00 
##            perimeter_se                 area_se           smoothness_se 
##            2.866059e+00            4.033708e+01            7.040979e-03 
##          compactness_se            concavity_se       concave.points_se 
##            2.547814e-02            3.189372e-02            1.179614e-02 
##             symmetry_se    fractal_dimension_se            radius_worst 
##            2.054230e-02            3.794904e-03            1.626919e+01 
##           texture_worst         perimeter_worst              area_worst 
##            2.567722e+01            1.072612e+02            8.805831e+02 
##        smoothness_worst       compactness_worst         concavity_worst 
##            1.323686e-01            2.542650e-01            2.721885e-01 
##    concave.points_worst          symmetry_worst fractal_dimension_worst 
##            1.146062e-01            2.900756e-01            8.394582e-02
```

```r
apply(wisc.data,2,sd)
```

```
##             radius_mean            texture_mean          perimeter_mean 
##            3.524049e+00            4.301036e+00            2.429898e+01 
##               area_mean         smoothness_mean        compactness_mean 
##            3.519141e+02            1.406413e-02            5.281276e-02 
##          concavity_mean     concave.points_mean           symmetry_mean 
##            7.971981e-02            3.880284e-02            2.741428e-02 
##  fractal_dimension_mean               radius_se              texture_se 
##            7.060363e-03            2.773127e-01            5.516484e-01 
##            perimeter_se                 area_se           smoothness_se 
##            2.021855e+00            4.549101e+01            3.002518e-03 
##          compactness_se            concavity_se       concave.points_se 
##            1.790818e-02            3.018606e-02            6.170285e-03 
##             symmetry_se    fractal_dimension_se            radius_worst 
##            8.266372e-03            2.646071e-03            4.833242e+00 
##           texture_worst         perimeter_worst              area_worst 
##            6.146258e+00            3.360254e+01            5.693570e+02 
##        smoothness_worst       compactness_worst         concavity_worst 
##            2.283243e-02            1.573365e-01            2.086243e-01 
##    concave.points_worst          symmetry_worst fractal_dimension_worst 
##            6.573234e-02            6.186747e-02            1.806127e-02
```

#Using PCA on wisc.data

```r
wisc.pr <- prcomp(wisc.data, scale = TRUE)
summary(wisc.pr)
```

```
## Importance of components:
##                           PC1    PC2     PC3     PC4     PC5     PC6
## Standard deviation     3.6444 2.3857 1.67867 1.40735 1.28403 1.09880
## Proportion of Variance 0.4427 0.1897 0.09393 0.06602 0.05496 0.04025
## Cumulative Proportion  0.4427 0.6324 0.72636 0.79239 0.84734 0.88759
##                            PC7     PC8    PC9    PC10   PC11    PC12
## Standard deviation     0.82172 0.69037 0.6457 0.59219 0.5421 0.51104
## Proportion of Variance 0.02251 0.01589 0.0139 0.01169 0.0098 0.00871
## Cumulative Proportion  0.91010 0.92598 0.9399 0.95157 0.9614 0.97007
##                           PC13    PC14    PC15    PC16    PC17    PC18
## Standard deviation     0.49128 0.39624 0.30681 0.28260 0.24372 0.22939
## Proportion of Variance 0.00805 0.00523 0.00314 0.00266 0.00198 0.00175
## Cumulative Proportion  0.97812 0.98335 0.98649 0.98915 0.99113 0.99288
##                           PC19    PC20   PC21    PC22    PC23   PC24
## Standard deviation     0.22244 0.17652 0.1731 0.16565 0.15602 0.1344
## Proportion of Variance 0.00165 0.00104 0.0010 0.00091 0.00081 0.0006
## Cumulative Proportion  0.99453 0.99557 0.9966 0.99749 0.99830 0.9989
##                           PC25    PC26    PC27    PC28    PC29    PC30
## Standard deviation     0.12442 0.09043 0.08307 0.03987 0.02736 0.01153
## Proportion of Variance 0.00052 0.00027 0.00023 0.00005 0.00002 0.00000
## Cumulative Proportion  0.99942 0.99969 0.99992 0.99997 1.00000 1.00000
```

#Q4. From your results, what proportion of the original variance is captured by the first principal components (PC1)?
There is a .4427 original variance

#Q5. How many principal components (PCs) are required to describe at least 70% of the original variance in the data?
All the way to PC3 will have = 72.6%

#How many principal components (PCs) are required to describe at least 90% of the original variance in the data?
You would need all the way up to PC7 to get 91% of the original variance

##Interpreting the Results

```r
#Try to plot our info
biplot(wisc.pr, scale = TRUE)
```

![](Class09_worksheet_files/figure-html/unnamed-chunk-10-1.png)<!-- -->
#Q7 What stands out to you about this plot? Is it easy or difficult to understand? Why?
This plot is too convoluted with data to actually read any individual point.

```r
# Scatter plot observations by components 1 and 2
attributes
```

```
## function (obj)  .Primitive("attributes")
```

```r
#We want "X"
plot(wisc.pr$x[, c(1,2)], col = (diagnosis +1),
     xlab = "PC1", ylab = "PC2")
```

![](Class09_worksheet_files/figure-html/unnamed-chunk-11-1.png)<!-- -->

#Q8. Repeat the same for principal components 1 and 3. What do you notice about these plots?


```r
# Repeat for components 1 and 3
plot(wisc.pr$x[, c(1, 3)], col = (diagnosis + 1), 
     xlab = "PC1", ylab = "PC3")
```

![](Class09_worksheet_files/figure-html/unnamed-chunk-12-1.png)<!-- -->
All the data is now on a higher range vertially. There is a greater spread


```r
pr.var <- wisc.pr$sdev^2
pve <- pr.var/sum(pr.var)
# Alternative scree plot of the same data, note data driven y-axis
barplot(pve, ylab = "Precent of Variance Explained", 
#Names.arg labels x axis values (PC#)
     names.arg=paste0("PC",1:length(pve)), las=2, axes = FALSE)
axis(2, at=pve, labels=round(pve,2)*100 )
```

![](Class09_worksheet_files/figure-html/unnamed-chunk-13-1.png)<!-- -->

##Variance Scree Plot

```r
#Calculate varianceof each component
pr.var <- wisc.pr$sdev^2
pve <- pr.var/sum(pr.var)

# Plot cumulative proportion of variance explained
plot( pve , xlab = "Principal Component", 
     ylab = "Cumulative Proportion of Variance Explained", 
     ylim = c(0, 1), type = "o")
```

![](Class09_worksheet_files/figure-html/unnamed-chunk-14-1.png)<!-- -->

##Clustering

```r
wisc.clust <-hclust(dist(scale(wisc.data)))

wisc.hclust.clusters <- cutree(wisc.clust, k=4)

#Plot the histogram
plot(wisc.clust)
abline(h = 20, col = "red")
```

![](Class09_worksheet_files/figure-html/unnamed-chunk-15-1.png)<!-- -->

How may M or 1 (cancer) and 0 (non-cancer) arein each cluster?


```r
table(wisc.hclust.clusters, diagnosis)
```

```
##                     diagnosis
## wisc.hclust.clusters   0   1
##                    1  12 165
##                    2   2   5
##                    3 343  40
##                    4   0   2
```

#Q12. Can you find a better cluster vs diagnoses match with by cutting into a different number of clusters between 2 and 10?

##K means

```r
wisc.km <- kmeans(scale(wisc.data), centers= 2, nstart = 20)
table(wisc.km$cluster, diagnosis)
```

```
##    diagnosis
##       0   1
##   1  14 175
##   2 343  37
```



Is this better than histogram?


```r
table(wisc.hclust.clusters , wisc.km$cluster)
```

```
##                     
## wisc.hclust.clusters   1   2
##                    1 160  17
##                    2   7   0
##                    3  20 363
##                    4   2   0
```

##Try again except with more PC values


```r
## Use the distance along the first 7 PCs for clustering i.e. wisc.pr$x[, 1:7]
wisc.pr.hclust <- hclust(dist(wisc.pr$x[,1:7]), method= "complete")
plot(wisc.pr.hclust)
```

![](Class09_worksheet_files/figure-html/unnamed-chunk-19-1.png)<!-- -->
Cut the tree down

```r
wisc.pr.hclust.clusters <- cutree(wisc.pr.hclust, k= 4)

table(wisc.pr.hclust.clusters, diagnosis)
```

```
##                        diagnosis
## wisc.pr.hclust.clusters   0   1
##                       1   5 113
##                       2 350  97
##                       3   2   0
##                       4   0   2
```
Compare all 3 methods

```r
table(wisc.km$cluster, diagnosis)
```

```
##    diagnosis
##       0   1
##   1  14 175
##   2 343  37
```

```r
table(wisc.hclust.clusters, diagnosis)
```

```
##                     diagnosis
## wisc.hclust.clusters   0   1
##                    1  12 165
##                    2   2   5
##                    3 343  40
##                    4   0   2
```

##Specificity vs. Selectivity
Sensitivity
- km seems to have the best sensitivity with a score of .7783 on sucessfully identifying sick patients
Selectivity
- The PC7 trial had the best specificity by rejecting most of the benign cases
