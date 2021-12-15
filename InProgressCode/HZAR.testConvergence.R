####Plots of BGC
##test for convergence

setwd("/Users/alexjvr/2016RADAnalysis/1.2_Phylo/BGC/EASTall.newinput_20161006")

##HI
hi1 <- read.csv("hiEASTall.1.out", header = F)

colnames(hi1) <- c("y1", "y2", "y3", "y4")  ##add column names
head(hi1)

hi1$x <- 1:nrow(hi1) ##add column x
head(hi1)


library("ggplot2")
library("reshape")

hi1.melt <- melt(hi1, id.vars=5)
head(hi1.melt)

ggplot(hi1.melt, aes(x, value, colour=variable)) +
  geom_point() +
  scale_color_manual(values = c("black", "black", "black", "black"))+
  xlab("mcmc") +
  ylab("Hybrid Index") +
  ggtitle("EASTsubset HI -I 0; -u 0.05 run6.3b") +
  theme(legend.position="none")


##alpha1
a1 <- read.csv("aEASTall.1.out", header = F)

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
  ggtitle("EASTsubset alpha -I 0 -u 0.2 run6.1a") +
  theme(legend.position="none")

##beta1
b1 <- read.csv("bEASTall.1.out", header = F)

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
  ggtitle("EASTsubset beta -I 0 -u 0.2 run6.1a") +
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

LnL1 <- read.tcsv("LnL.EASTall.1.out")
head(LnL1)

LnL1$x <- 1:nrow(LnL1) ##add column x
head(LnL1)

ggplot(LnL1, aes(x, LnL1$LnL, colour="black")) +
  geom_point() +
  scale_color_manual(values = c("black"))+
  xlab("mcmc") +
  ylab("LnL1") +
  ggtitle("EASTsubset Log likelihood run8.6") +
  theme(legend.position="none")

