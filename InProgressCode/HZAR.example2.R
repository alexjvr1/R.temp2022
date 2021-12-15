##HZAR example to run
##example data to use
#1. snp.freq.vcf = data file
#2. example.loci.csv = list of loci



setwd("~/2016RADAnalysis/1.2_Phylo/HZAR/")


library(hzar)
input.hzar <- read.table("snp.frq", header = T)

write.table(input.hzar, # The data we just loaded
            file="input.hzar.txt", # The file to overwrite
            col.names=TRUE, # The columns are named
            row.names=FALSE, # The rows are not named
            sep="\t", # The file will be tab-delimited
            quote=TRUE) # Use quotes as needed.

summary(input.hzar)
## As we no longer need the in-memory copy, drop the local reference
input.hzar <- NULL

## Save all plots in a series of png files
png(width=900, height=900, res=200, family="Arial", filename="mExPlot
    %03d.png",pointsize=8)



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
input.hzar <- read.table("input.hzar.txt",header=TRUE)

## Print data
print(input.hzar)

## ## Picking an allele for a locus
## useAlleles <- "ada.A";
## ua.nSamples <- "ada.nSamples"

## Blank out space in memory to hold molecular analysis
if(length(apropos("^hzar.test$",ignore.case=FALSE)) == 0 ||
   !is.list(input.hzar) ) hzar.test <- list()

## We are doing just the one allele at one locus, but it is
## good to stay organized.
## Space to hold the observed data
#CZ$CHN <- list();
#CZ$CHN$obs <- list();
## Space to hold the models to fit
#CZ$CHN$models <- list();
## Space to hold the compiled fit requests
#CZ$CHN$fitRs <- list();
## Space to hold the output data chains
#CZ$CHN$runs <- list();

## Space to hold the analysed data
#CZ$CHN$analysis <- list();

#head(CZ)

#CZ$CHN$obs <-
 # hzar.doMolecularData1DPops(EAST$distance,
  #                           EAST$mtDNA.CHN1,
   #                          EAST$mtDNA.nAlleles); 


###List of all the alleles to be analysed
##see http://elizabethderryberry.tulane.edu/derryberrylab/Software/Entries/2014/10/1_Boilerplate_loader_code_for_molecular_data.html

loci <- read.table("example.loci")$V1  ##read in a list of the loci names 
loci <- as.character(loci)
loci

obsData.example <-
  sapply(loci,
         function(snpIter) hzar.doMolecularData1DPops(
          input.hzar$distance,
           input.hzar[[paste(snpIter,"A",sep=".")]],
           input.hzar[[paste(snpIter,"nSamples",sep=".")]],
           rownames(input.hzar)),
         simplify=FALSE)
summary(obsData.example)




##BOILERPLATE FOR COMPILING MODELS FOR MOLECULAR DATA
##http://elizabethderryberry.tulane.edu/derryberrylab/Software/Entries/2014/10/1_Boilerplate_for_compiling_models_for_molecular_data.html
##This creates an object with 15 models specified to be fit to the molecular data


## Use parallel to speed code up: 
library(parallel)

## Determine range of observed data 
##   (assumes site distances in 'snp.freq$distance'
obs.range <- extendrange(input.hzar$distance)


## This code assumes that:
##  'loci' are the loci you wish to load (eg X66, X95, X15)
##  'obsData' is a named list of hzar.obsData objects, 
##     such as the result of the boilerplate loader code

## Assume cache of result in "cache.fq.models.dat.gz"
if(file.exists("example.cache.fq.models.dat.gz")){
  load("example.cache.fq.models.dat.gz") }else{
    fq.models <- local({
      mkFQModel <- function(locusID,scaling,tails){
        res <-
          hzar.model.addBoxReq(hzar.makeCline1DFreq(
            data=obsData.example[[locusID]],
            scaling=scaling,tails=tails),
            low=obs.range[[1]],high=obs.range[[2]])
        hzar.first.fitRequest.old.ML(res,obsData.example[[locusID]],verbose=FALSE)
      }
      allModels <- expand.grid(s=c("none","fixed","free"),
                               t=c("none","left","right","mirror","both"))
      mkFQModelA <- function(s,t){
        res <- mclapply(loci,mkFQModel,scaling=s,tails=t,mc.cores=3)
        
        ##res <- lapply(res, hzar.first.fitRequest.old.ML)
        names(res) <- paste(loci,s,t,sep=".")
        res
      }
      res <- list()
      for(iter in 1:nrow(allModels))
        res <- c(res,mkFQModelA(allModels$s[[iter]],allModels$t[[iter]]))
      print(    res.names <- names(res))
      res <- hzar.multiFitRequest(res)
      res.names -> names(res)
      res
    })
    ## Save cache of result in "cache.fq.models.dat.gz"
    save(fq.models,file="example.cache.fq.models.dat.gz") 
  }


summary(fq.models)

##RUNNING MODELS & MODEL SELECTION
###Boilerplate code to run models, do model selection, and then analyze the selected models (using the subplex helper snippet). Includes caching.

## This code assumes:
##  'obsData.example' is a named list of hzar.obsData objects, 
##     such as the result of the boilerplate loader code
##  'fq.models' is named list of hzar.fitRequest objects,
##     such as the result of the boilerplate molecular model compile code
library(doMC)
registerDoMC(3)

##started at 14:30-14:57
if(file.exists("example.cache.fq.init.run.dat.gz")){
  load("example.cache.fq.init.run.dat.gz")
}else{
  hzar.doFit.multi(fq.models,doPar=TRUE)->fq.init.run;
  save(fq.init.run,file="example.cache.fq.init.run.dat.gz")
}

###Functions that need defining

library(hzar)
#install.packages("subplex")
library(subplex)

temp.doMLE <- function(fitR){
  mP <- fitR$modelParam
  res <- subplex(mP$init,
                 function(theta) -fitR$llFunc(theta),
                 control=list(parscale=sapply(names(mP$init),
                                              function(x) (mP$upper[[x]]-mP$lower[[x]])/50)))
  if(attr(fitR, "fit.success"))
    fitR$mcmcRaw <- rbind(fitR$mcmcRaw ,res$par)
  else{
    fitR$mcmcRaw <- res$par
    fitR$modelParam$init <- as.list(res$par)
  }
  fitR
}


temp.getMLE <- function(fitR){
  mP <- fitR$modelParam
  res <- subplex(mP$init,
                 function(theta) -fitR$llFunc(theta),
                 control=list(parscale=sapply(names(mP$init),
                                              function(x) (mP$upper[[x]]-mP$lower[[x]])/50)))
  hzar.gen.cline(res$par, fitR)
}

###
##FIT MODELS  started 15:03 - 15:04
fq.sel.dG <- list()

summary(obsData.example)

for(iter in 1:length(obsData.example)) {
  fq.tmp <- local({
    load("example.cache.fq.init.run.dat.gz")
    foreach(run=fq.init.run) %:%
      when(hzar.sameObsData(run,obsData.example[[iter]])) %dopar% hzar.fit2DataGroup(temp.doMLE(run))
  })
  fq.tmp <- hzar.copyModelLabels(fq.models,
                                 hzar.make.obsDataGroup( fq.tmp ));
  
  
  print(fq.aic <- hzar.AICc.hzar.obsDataGroup(fq.tmp,show.count=TRUE));
  fq.sel.dG <- c( fq.sel.dG,fq.tmp$data.groups[ which.min(fq.aic$AICc)]);
  print(names(fq.sel.dG))
  rm(fq.tmp)
  print(gc())
}

summary()
head(fq.init.run)

##started 15:04 - 15:11

if(file.exists("example.cache.fq.chains.dat.gz")){
  load("example.cache.fq.chains.dat.gz")
}else{
  names(fq.models) -> names(fq.init.run)
  fq.sel.run <- fq.init.run[names(fq.sel.dG)]
  fq.chains.init <- foreach(fq.run=fq.sel.run) %dopar% {
    hzar.next.fitRequest(fq.run) }
  fq.chains.init <- lapply(fq.chains.init,function(x) {
    x$modelParam$tune <- sapply(x$modelParam$tune,
                                function(y) 0.8*y, simplify=FALSE);
    x })
  fq.chains <-
    hzar.doChain.multi(doPar=TRUE,
                       hzar.multiFitRequest( fq.chains.init,
                                             each=3,
                                             baseSeed=NULL,
                                             baseChannel=60,adjChannel=10))
  save(fq.chains,file="example.cache.fq.chains.dat.gz")
}




