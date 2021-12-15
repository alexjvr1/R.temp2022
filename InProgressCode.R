###BGC plots

setwd("/Users/alexjvr/2016RADAnalysis/1.2_Phylo/BGC/EASTrun4/")

##alpha1
a1 <- read.csv("a.4b.out", header = F)

colnames(a1) <- c("mean", "median", "CI_LB", "CI_UB")  ##add column names
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
  ggtitle("EASTsubset alpha BGC run4.2") +
  theme(legend.position="none")

##beta1
b1 <- read.csv("b.4b.out", header = F)
head(b1)

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
  ggtitle("EASTsubset beta BGC run4.2") +
  theme(legend.position="none")

##LnL1  
##This uses a function tcsv that transforms the data and reads in the header. I.e. can keep the numeric format, while specifying a text header

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

##nano, add "LnL," to start of the file

LnL1 <- read.tcsv("LnL.4b.out")
head(LnL1)

LnL1$x <- 1:nrow(LnL1) ##add column x
head(LnL1)

ggplot(LnL1, aes(x, LnL1$LnL, colour="black")) +
  geom_point() +
  scale_color_manual(values = c("black"))+
  xlab("mcmc") +
  ylab("LnL1") +
  ggtitle("EASTsubset Log likelihood BGC run4.2") +
  theme(legend.position="none")


######
###Assess results

##Plot alpha
##excess loci

setwd("/Users/alexjvr/2016RADAnalysis/1.2_Phylo/BGC/EASTrun4/")

##alpha1
a1 <- read.csv("a.4b.out", header = F)

colnames(a1) <- c("mean", "median", "CI_LB", "CI_UB")  ##add column names
head(a1)

a1$x <- 1:nrow(a1) ##add column x
head(a1)

CIlimits <- aes (ymax=a1$CI_UB, ymin=a1$CI_LB)

library("ggplot2")
library("reshape")

ggplot(a1, aes(x=a1$x, y=a1$mean)) +
  geom_point() +
  #geom_pointrange(CIlimits, width=0.2) +
  scale_color_manual(values = c("black", "black", "black", "black"))+
  xlab("x") +
  ylab("alpha + CI") +
  ggtitle("EASTsubset alpha BGC run4 + CI") +
  theme(legend.position="none")


##Identify outliers
##2.5% tails

setwd("/Users/alexjvr/2016RADAnalysis/1.2_Phylo/BGC/EASTrun4/")

a1 <- read.csv("a.4b.out", header = F)

colnames(a1) <- c("mean", "median", "CI_LB", "CI_UB")  ##add column names
head(a1)

a1$x <- 1:nrow(a1) ##add column x
head(a1)


library("ggplot2")
library("reshape")
#Calculate the variance

alpha.postDist <- a1$mean

tau.alpha <- (sqrt(var(alpha.postDist)))
tau.alpha  ##SD used for the calculation of outliers

alpha.pos <- tau.alpha*1.96
alpha.pos

alpha.neg <- 0-(tau.alpha*1.96)
alpha.neg



CIlimits <- aes (ymax=a1$CI_UB, ymin=a1$CI_LB)

library("ggplot2")
library("reshape")

ggplot(a1, aes(x=a1$x, y=a1$mean, color=a1$mean)) +
  geom_point(aes(color = cut(a1$mean, c(-Inf, -26.87914, 26.87924, Inf)))) +
  scale_color_manual(name = "mean", values = c("(-Inf,-26.87914]" = "green",
                                               "(-26.87914,26.87914]" = "black",
                                               "(26.87914, Inf]" = "purple")) +
  xlab("x") +
  ylab("alpha + CI") +
  ggtitle("EASTsubset alpha BGC run4 + CI") +
  theme(legend.position="none")

#scale_color_manual(colours = c("green", "black", "purple"), values=c(Inf,26.87914,-26.87914,-Inf)) +
#geom_pointrange(CIlimits, width=0.2) +




