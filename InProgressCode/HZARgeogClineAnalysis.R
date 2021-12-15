##HZAR geographic cline analysis

#What does the input file look like? 

library(hzar)

setwd("~/2016RADAnalysis/1.2_Phylo/HZAR/")
freq <- read.table("plink.frq.strat.2.txt", header = T, stringsAsFactors = F)
head(freq)
str(freq)

#freq <- as.data.frame(freq, stringsAsFactors = F)
#str(freq)


freq$CHR = NULL##get rid of some of the columns
freq$A1 = NULL
freq$A2 = NULL
freq$MAC = NULL
##freq$MAF = NULL ##for the test

##Step 2
#freq$n.samples <- with(freq, NCHROBS/2) #Calculate the number of individuals typed for each locus per pop
#head(freq)
freq$NCHROBS = NULL #delete this column

#colnames(freq)[colnames(freq)=="MAF"] <- "B" ##rename the MAF to B. (This specifies a particular column)


#freq <- freq[,c(1,2,4,3)]#change the order of the columns
#head(freq)


library("ggplot2")
library("reshape2")

freq2 <- melt(freq, id.vars = c("CLST", "SNP"), variable_name = c("MAF"))
str(freq2)
head(freq2)


freq3 <- dcast(freq2, formula= CLST ~ SNP)
head(freq3)
write.csv(freq3, file="freq3.csv")


##Now do the sampe with n.samples
##Step 2: nsamples

freq <- read.table("plink.frq.strat.2.txt", header = T, stringsAsFactors = F)
head(freq)

freq$n.samples <- with(freq, NCHROBS/2) #Calculate the number of individuals typed for each locus per pop
head(freq)
freq$NCHROBS = NULL #delete this column

freq$CHR = NULL##get rid of some of the columns
freq$A1 = NULL
freq$A2 = NULL
freq$MAC = NULL
freq$MAF = NULL

head(freq)

library("ggplot2")
library("reshape2")


freq2 <- melt(freq, id.vars = c("CLST", "SNP"), variable_name = c("n.samples"))
str(freq2)
head(freq2)


freq3 <- dcast(freq2, formula= CLST ~ SNP)
head(freq3)
write.csv(freq3, file="nsamples.csv")

##Now do the sampe with n.samples
##Step 3: Major Allele freq

freq <- read.table("plink.frq.strat.2.txt", header = T, stringsAsFactors = F)
head(freq)

freq$A <- with(freq, 1-MAF) ##add column with major allele freq.
head(freq) 

freq$NCHROBS = NULL #delete this column
freq$CHR = NULL##get rid of some of the columns
freq$A1 = NULL
freq$A2 = NULL
freq$MAC = NULL
freq$MAF = NULL

head(freq)

library("ggplot2")
library("reshape2")


freq2 <- melt(freq, id.vars = c("CLST", "SNP"), variable_name = c("A"))
str(freq2)
head(freq2)


freq3 <- dcast(freq2, formula= CLST ~ SNP)
head(freq3)
write.csv(freq3, file="MajorAllele.csv")





