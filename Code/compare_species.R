d <- read.csv("compare.csv", header = FALSE)
long <- d[,1]
short <-d[,2]

long <- sub("_", " ", long)

short <- short[! short %in% c("")]

x <- setdiff(long, short)


