##Run HZAR
##import final input file

setwd("~/2016RADAnalysis/1.2_Phylo/HZAR/")

library(hzar)

hzar1 <- read.table("InputEAST.txt", header = T)
head(hzar1)

#create .csv with only allele frequencies included
Allelefreq <- read.csv("Allelefreq.csv", header=T)
summary(Allelefreq)


for (i in Allelefreq){
  EAST.Allele01 <- hzar.doMolecularData1DPops(hzar1$Dist,
                                              Allelefreq$i,
                                              hzar1$nsamples)}
hzar.plot.obsData(EAST.Allele01);

HZAR <- function(file, col){
  EAST.A1model <- hzar.makeCline1DFreq(EAST.Allele01, scaling="fixed",tails="none");
  EAST.A1Amodel <-
    hzar.model.addBoxReq(EAST.A1model,-30,600);
  EAST.A1modelFitR <-
    hzar.first.fitRequest.old.ML(model=EAST.A1model ,
                                 EAST.Allele01,
                                 verbose=FALSE);
  EAST.A1modelFitR$mcmcParam$chainLength <- 2e3;
  EAST.A1modelFitR$mcmcParam$burnin <- 5e2;
  EAST.A1modelFit <- hzar.doFit(EAST.A1modelFitR)
  plot(hzar.mcmc.bindLL(EAST.A1modelFit))
  EAST.A1modelData <-
    hzar.dataGroup.add(EAST.A1modelFit);
  ## Not run:
  EAST.A1modelData <-
    hzar.dataGroup.add(
      EAST.A1modelData,
      hzar.chain.doSeq(hzar.next.fitRequest(EAST.A1modelFit)));
  hzar.plot.cline(EAST.A1modelData);
  hzar.plot.fzCline(EAST.A1modelData);
  
  ## End(Not run)
  print(hzar.getLLCutParam(EAST.A1modelData,c("center","width")));
  EAST.A1modelNull <- hzar.dataGroup.null(EAST.Allele01);
  EAST.A1AdGs <- list(clineModel = EAST.A1modelData,
                      nullModel = EAST.A1modelNull);
  EAST.A1oDG <- hzar.make.obsDataGroup(EAST.A1AdGs);
  EAST.A1oDG <- hzar.copyModelLabels(EAST.A1AdGs,EAST.A1oDG);
  hzar.plot.cline(EAST.A1oDG);
  print(hzar.AICc.hzar.obsDataGroup(EAST.A1oDG));
}

for (i in 1:ncol(hzar1)){
  HZAR(file,col=i)
}


