---
title: "Bioniformatics Class 7"
author: "Matthew Tomaneng"
date: "4/24/2018"
output: 
  html_document: 
    keep_md: yes
---


##Functions Again
Revisiting funtions from class 6

```r
source("http://tinyurl.com/rescale-R")
```

Lets see if we can use this function


```r
rescale(1:10)
```

```
##  [1] 0.0000000 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556 0.6666667
##  [8] 0.7777778 0.8888889 1.0000000
```

Looks good lets break it!


```r
#rescale(c(1:10, "string"))
```

Lets try the new **rescale2** function


```r
#rescale2(c(1:10, "string"))
```

##Lets make a function to count number of NAs
Write a new function to check for NAs in two inputs


```r
#Start by searching for an existing function that does what you want
#example data
x <- c( 1, 2, NA, 3, NA)
y <- c(NA, 3, NA, 3, 4)
```



```r
which(is.na(x))
```

```
## [1] 3 5
```
How many total NAs

```r
sum(is.na(x) & is.na(y))
```

```
## [1] 1
```

Let's test if this works for both vectors at once

```r
sum(is.na(c(x, y)))
```

```
## [1] 4
```
#It works!
Now Let's Combine the Functions to make our function

```r
num_na <- function(x, y)
  #Number of NAs in the vectors
  sum(is.na(c(x, y)))
```
Let's test it this works


```r
num_na(x, y)
```

```
## [1] 4
```

Test if it works for other situations

```r
x <- c(NA, NA, NA)
y1 <- c( 1, NA, NA)
y2 <- c( 1, NA, NA, NA)
num_na(x, y2)
```

```
## [1] 6
```

I messed up, doing it from both_na3

```r
x <- c( 1, 2, NA, 3, NA)
y <- c(NA, 3, NA, 3, 4)
both_na3(x, y)
```

```
## Found 1 NA's at position(s):3
```

```
## $number
## [1] 1
## 
## $which
## [1] 3
```

##Function for gene intersection

```r
df1
```

```
##     IDs exp
## 1 gene1   2
## 2 gene2   1
## 3 gene3   1
```

```r
df2
```

```
##     IDs exp
## 1 gene2  -2
## 2 gene4  NA
## 3 gene3   1
## 4 gene5   2
```

```r
x <- df1$IDs
y <- df2$IDs
```


```r
x
```

```
## [1] "gene1" "gene2" "gene3"
```

```r
y
```

```
## [1] "gene2" "gene4" "gene3" "gene5"
```

Found a function called Intersect

```r
intersect(x, y)
```

```
## [1] "gene2" "gene3"
```
Not good enough, don't have position data

```r
#Value Matching funciton
x%in%y
```

```
## [1] FALSE  TRUE  TRUE
```


```r
y%in%x
```

```
## [1]  TRUE FALSE  TRUE FALSE
```
Getting a little more information when it comes to position of match

```r
#We can do better by getting the actual names of the genes
x[x%in%y]
```

```
## [1] "gene2" "gene3"
```
Since X is a vector, putting the values into sqaure brackets will return the value


```r
#Now lets combine this info into a single table
#Put info into two columns 
cbind(x[x%in%y], y[y%in%x])
```

```
##      [,1]    [,2]   
## [1,] "gene2" "gene2"
## [2,] "gene3" "gene3"
```

So now lets make this into a function


```r
#function(x, y) { 
#   cbind( x[ x %in% y ], y[ y %in% x ] )
#} Already made in the enviornment
#Not very clear as well
```

Lets try the merge function

```r
merge(df1, df2, by = "IDs")
```

```
##     IDs exp.x exp.y
## 1 gene2     1    -2
## 2 gene3     1     1
```
##MERGE IS A SUPER HELPFUL FUNCTION













