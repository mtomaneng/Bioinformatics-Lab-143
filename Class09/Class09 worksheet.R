wisc.df <- read.csv("~/Downloads/WisconsinCancer.csv")

#converting features of wisc.df
wisc.data <- as.matrix(wisc.df)

#Set the row names of wisc.data
row.names(wisc.data) <- wisc.df$id
head(wisc.data)

#Making a new vector called diagnosis
diagnosis <- as.numeric(wisc.data$diagnosis == "M", )
