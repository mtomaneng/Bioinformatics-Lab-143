# Class 6 Bioinformatics

# Functions

add <- function(x, y =1) {
  #sum x and y inputs
  x + y
}

# Rescaling a range of values
rescale <- function(x) {
  rng <-range(x, na.rm = TRUE)
  (x - rng[1]) / (rng[2] - rng[1])
}

#Update Rescale to do work more often
rescale2 <- function(x, na.rm=TRUE, plot=FALSE)  {
  if(na.rm) {
    rng <-range(x, na.rm=na.rm)
  } else {
    rng <-range(x)
  }
  print("Hello")
  answer <- (x - rng[1]) / (rng[2] - rng[1])
  print("is it me you are looking for?")
  if(plot) {
    plot(answer, typ="b", lwd=4)
  }
  print("I can see it in ...")
}

#3rd version with return function inside
rescale3 <- function(x, na.rm=TRUE, plot=FALSE) {
  if(na.rm) {
    rng <-range(x, na.rm=TRUE)
  } else {
    rng <-range(x)
  }
  print("Hello")
  answer <- (x - rng[1]) / (rng[2] - rng[1])
  return(answer)
  print("is it me you are looking for?")
  if(plot) {
    plot(answer, typ="b", lwd=4)
  }
  print("I can see ")
}
#Test it out
rescale(c(1,2,NA, 3,10))

