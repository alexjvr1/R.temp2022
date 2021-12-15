##HZAR example to run

setwd("~/2016RADAnalysis/1.2_Phylo/HZAR/")

library(hzar)
data("manakinMolecular")
print(manakinMolecular)

write.table(manakinMolecular, # The data we just loaded
            file="mknExLoci.txt", # The file to overwrite
            col.names=TRUE, # The columns are named
            row.names=FALSE, # The rows are not named
            sep="\t", # The file will be tab-delimited
            quote=TRUE) # Use quotes as needed.

## As we no longer need the in-memory copy, drop the local reference
manakinMolecular <- NULL

## Save all plots in a series of png files
#png(width=900, height=900, res=200, family="Arial", filename="mExPlot
#%03d.png",pointsize=8)



## A typical chain length. This value is the default setting in the package.
chainLength=1e5;

##each model runs off a different seed
mainSeed=
  list(A=c(596,528,124,978,544,99), 
       B=c(528,124,978,544,99,596), 
       C=c(124,978,544,99,596,528))

#install.packages("doMC")
if(require(doMC)){
## If you have doMC, use foreach in parallel mode
## to speed up computation.
registerDoMC() } else {## Use foreach in sequential mode
registerDoSEQ();
}

## Molecular Analysis

## Load example Molecular data from the data table.
manakinMolecular <- read.table("mknExLoci.txt",header=TRUE)
  
## Print sample data
print(manakinMolecular)

## ## Picking an allele for a locus
## useAlleles <- "ada.A";
## ua.nSamples <- "ada.nSamples"

 ## Blank out space in memory to hold molecular analysis
 if(length(apropos("^mkn$",ignore.case=FALSE)) == 0 ||
    !is.list(mkn) ) mkn <- list()

snp.freq <- read.table("snp.freq.txt", header = T)
print(snp.freq)


loci2 <- read.csv("loci.mkn", header=F)$V1
loci3 <- as.character(loci2)
loci3


obsData <-
  sapply(loci3,
          function(snpIter) hzar.doMolecularData1DPops(
            snp.freq$distance,
            snp.freq[[paste(snpIter,"A",sep=".")]],
            snp.freq[[paste(snpIter,"nSamples",sep=".")]],
            rownames(snp.freq)),
          simplify=FALSE)
summary(obsData)





## We are doing just the one allele at one locus, but it is
 ## good to stay organized.
 mkn$AdaA <- list();
 ## Space to hold the observed data
   mkn$AdaA$obs <- list();
 ## Space to hold the models to fit
   mkn$AdaA$models <- list();
 ## Space to hold the compiled fit requests
   mkn$AdaA$fitRs <- list();
 ## Space to hold the output data chains
   mkn$AdaA$runs <- list();
 
 ## Space to hold the analysed data
   mkn$AdaA$analysis <- list();

head(mkn)
 
 mkn$AdaA$obs <-
   hzar.doMolecularData1DPops(manakinMolecular$distance,
                                 manakinMolecular$ada.A,
                                 manakinMolecular$ada.nSamples); 
 ##graph of the observed data
 hzar.plot.obsData(mkn$AdaA$obs);
 
 ##Make helper function
 
 mkn.loadAdaAmodel <- function(scaling,tails,
                               id=paste(scaling,tails,sep="."))
   mkn$AdaA$models[[id]] <<- hzar.makeCline1DFreq(mkn$AdaA$obs, scaling, tails)

 mkn.loadAdaAmodel("fixed","none","modelI");
 mkn.loadAdaAmodel("free" ,"none","modelII");
 mkn.loadAdaAmodel("free" ,"both","modelIII");
 mkn.loadAdaAmodel("modelIV")
 
 #Check the default settings
print(mkn$AdaA$models)
 
##modify all the models to focus on the region where the data were collected
##example data this is between 0 and 570km
##CH.EAST 0-85km
 
 mkn$AdaA$models <- sapply(mkn$AdaA$models,
                            hzar.model.addBoxReq,
                            -30 , 600,
                            simplify=FALSE)
 
 ## Check the updated settings
 print(mkn$AdaA$models)
 
 ## Compile each of the models to prepare for fitting
 mkn$AdaA$fitRs$init <- sapply(mkn$AdaA$models,
                                 hzar.first.fitRequest.old.ML,
                                 obsData=mkn$AdaA$obs,
                                 verbose=FALSE,
                                 simplify=FALSE)
 
##results  ##update the settings for the fitter if desired
 
chainLength <- 1e5
 
 mkn$AdaA$fitRs$init$modelI$mcmcParam$chainLength <-
   chainLength; #1e5
 mkn$AdaA$fitRs$init$modelI$mcmcParam$burnin <-
   chainLength %/% 10; #1e4
 mkn$AdaA$fitRs$init$modelI$mcmcParam$seed[[1]] <-
   mainSeed$A
 

mkn$AdaA$fitRs$init$modelII$mcmcParam$chainLength <-
    chainLength; #1e5
 mkn$AdaA$fitRs$init$modelII$mcmcParam$burnin <-
    chainLength %/% 10; #1e4
 mkn$AdaA$fitRs$init$modelII$mcmcParam$seed[[1]] <-
    mainSeed$B
 
   mkn$AdaA$fitRs$init$modelIII$mcmcParam$chainLength <-
    chainLength; #1e5
 mkn$AdaA$fitRs$init$modelIII$mcmcParam$burnin <-
    chainLength %/% 10; #1e4
 mkn$AdaA$fitRs$init$modelIII$mcmcParam$seed[[1]] <-
    mainSeed$C
   
   
##Check fit request settings
   
print(mkn$AdaA$fitRs$init)

##Run just one of the models for an initial chain
mkn$AdaA$runs$init <- list()
mkn$AdaA$runs$init$modelI <- hzar.doFit(mkn$AdaA$fitRs$init$modelI)

##plot the trace
plot(hzar.mcmc.bindLL(mkn$AdaA$runs$init$modelI))

##Run another model for an initial chain
mkn$AdaA$runs$init$modelII <- hzar.doFit(mkn$AdaA$fitRs$init$modelII)

#Plot the trace
plot(hzar.mcmc.bindLL(mkn$AdaA$runs$init$modelII))
 
##Run another model for an initial chain
mkn$AdaA$runs$init$modelIII <- hzar.doFit(mkn$AdaA$fitRs$init$modelIII)

#Plot the trace
plot(hzar.mcmc.bindLL(mkn$AdaA$runs$init$modelIII))


