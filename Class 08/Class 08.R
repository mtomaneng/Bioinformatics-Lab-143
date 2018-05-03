#Unstructured learning: K clustering

# Generate some example data for clustering
tmp <- c(rnorm(30,-3), rnorm(30,3))
x <- cbind(x=tmp, y=rev(tmp))
plot(x)
# Use the kmeans() function setting k to 2 and nstart=20
# Inspect/print the results
# Q. How many points are in each cluster?
# Q. What ‘component’ of your result object details
# - cluster size?
# - cluster assignment/membership?
# - cluster center?
# Plot x colored by the kmeans cluster assignment and
# add cluster centers as blue points
# Q. Repeat for k=3, which one has the better total SS? 
