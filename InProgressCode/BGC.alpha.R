setwd("/Users/alexjvr/2016RADAnalysis/1.2_Phylo/BGC/EASTrun4")

##alpha1
a1 <- read.csv("a.4c.out", header = F)

colnames(a1) <- c("y1", "y2", "y3", "y4")  ##add column names
head(a1)

a1$x <- 1:nrow(a1) ##add column x
head(a1)


library("ggplot2")
library("reshape")

a1.melt <- melt(a1, id.vars=5)
head(a1.melt)

ggplot(a1.melt, aes(x, value, colour=variable)) +
  geom_point() +
  scale_color_manual(values = c("black", "black", "black", "black"))+
  xlab("mcmc") +
  ylab("alpha1") +
  ggtitle("EASTsubset alpha BGC run1") +
  theme(legend.position="none")

##beta1
b1 <- read.csv("b.4c.out", header = F)

colnames(b1) <- c("y1", "y2", "y3", "y4")  ##add column names
head(b1)

b1$x <- 1:nrow(b1) ##add column x
head(b1)


library("ggplot2")
library("reshape")

b1.melt <- melt(b1, id.vars=5)
head(b1.melt)

ggplot(b1.melt, aes(x, value, colour=variable)) +
  geom_point() +
  scale_color_manual(values = c("black", "black", "black", "black"))+
  xlab("mcmc") +
  ylab("beta1") +
  ggtitle("EASTsubset beta BGC run1") +
  theme(legend.position="none")

##LnL1  
##This uses a function tcsv that transforms the data and reads in the header. I.e. can keep the numeric format, while specifying a text header
###Add "LnL," to the start of the file using nano in bash

read.tcsv = function(file, header=TRUE, sep=",", ...) {
  
  n = max(count.fields(file, sep=sep), na.rm=TRUE)
  x = readLines(file)
  
  .splitvar = function(x, sep, n) {
    var = unlist(strsplit(x, split=sep))
    length(var) = n
    return(var)
  }
  
  x = do.call(cbind, lapply(x, .splitvar, sep=sep, n=n))
  x = apply(x, 1, paste, collapse=sep) 
  out = read.csv(text=x, sep=sep, header=header, ...)
  return(out)
  
}

LnL1 <- read.tcsv("LnL.4c.out")
head(LnL1)

LnL1$x <- 1:nrow(LnL1) ##add column x
head(LnL1)

ggplot(LnL1, aes(x, LnL1$LnL, colour="black")) +
  geom_point() +
  scale_color_manual(values = c("black"))+
  xlab("mcmc") +
  ylab("LnL1") +
  ggtitle("EASTsubset Log likelihood BGC run1") +
  theme(legend.position="none")

