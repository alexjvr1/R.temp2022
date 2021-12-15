####Check conversion script for BGC
#############

setwd("/Users/alexjvr/2016RADAnalysis/1.2_Phylo/input.files/BGCinput/")

genepop <- data.table::fread("genepop.BGC.txt",
                             header = FALSE, sep = "\t",
                             stringsAsFactors = FALSE)


header <- genepop[1,]
if(length(gregexpr(',', header, fixed=F)[[1]])>1){
  lociheader <- strsplit(header,",")
  lociheader <- gsub(" ","",unlist(lociheader))
  #remove the first column of loci names
  genepop <- as.vector(genepop)
  genepop <- genepop[-1,]
  genepop <- c(lociheader,genepop)
  genepop <- as.data.table(genepop,stringsAsFactors = FALSE)
}
header ###this is the first row. i.e. if the locus names are in one row and separated my ","

## Stacks version information
stacks.version <- genepop[1,] #this could be blank or any other source. First row is ignored by genepop


genepop <- genepop[-1,]
colnames(genepop) <- "data"  ##adds colname "data" to the first column
head(genepop)


Pops  <-  which(genepop$data == "Pop" | genepop$data =="pop" | genepop$data == "POP")
Pops ##check that the line numbers make sense. i.e. double the nr of indivs for each pop
npops  <-  1:length(Pops)
npops ##should be 3 if named Admix, P1, P2
  
## separate the data into the column headers and the rest
ColumnData <- genepop$data[1:(Pops[1]-1)]
ColumnData <- gsub("\r","",ColumnData)#remove any hidden carriage returns
summary(ColumnData)   ##length should be equal to nr of loci
snpData <- genepop[Pops[1]:NROW(genepop),]
head(snpData)

#Get a datafile with just the snp data no pops
tempPops <- which(snpData$data=="Pop"| snpData$data =="pop" | snpData$data == "POP") ## Changed because we allowed
tempPops
snpData <- snpData[-tempPops,]
head(snpData)

#separate the snpdata
temp <- as.data.frame(do.call(rbind, strsplit(snpData$data," ")))
summary(temp)

#data format check
if(unique(temp[,2])!="," | !length(which(temp[,3]==""))>1){
  stop("Genepop sampleID delimiter not in proper format. Ensure sampleIDs are separated from loci by ' ,  ' (space comma space space). Function stopped.",call. = FALSE)
}
temp2 <- temp[,4:length(temp)] #split characters by spaces

#Contingency to see if R read in the top line as the "stacks version"
if (length(temp2)!=length(ColumnData)){colnames(temp2) <- c(stacks.version,ColumnData)}
if (length(temp2)==length(ColumnData)){colnames(temp2) <- ColumnData}
if (length(temp2)!=length(ColumnData)){stacks.version="No STACKS version specified"}

#stacks version character
stacks.version <- as.character(stacks.version)

## Get the population names (prior to the _ in the Sample ID)
NamePops <- temp[,1] # Sample names of each
NameExtract <- substr(NamePops,1,regexpr("_",NamePops)-1)  ##each sample name should start with "pop_" for this to work


#convert the snp data into character format to get rid of factor levels
temp2[] <- lapply(temp2, as.character)

#allele coding length
alleleEx <- max(sapply(temp2[,1],FUN=function(x){nchar(as.character(x[!is.na(x)]))})) #presumed allele length
alelleEx ##should be 6 if you used pgdspider for the conversion to genepop


#check to make sure the allele length is a even number
if(!alleleEx %% 2 ==0){stop(paste("The length of each allele is assumed to be equal (e.g. loci - 001001 with 001 for each allele), but a max loci length of", alleleEx, "was detected. Please check data."))}

#get the allele values summary header
firstAllele <-  as.data.frame(sapply(temp2,function(x)as.numeric(as.character(substring(x,1,(alleleEx/2))))))
secondAllele <-  as.data.frame(sapply(temp2,function(x)as.numeric(as.character(substring(x,(alleleEx/2)+1,alleleEx)))))

# switch from combined allele in one row to two allele each in their own row according to a given locus
holdframe <- rbind(firstAllele,secondAllele)#create a dummy data frame
holdframe[!1:nrow(holdframe) %% 2 == 0,] = firstAllele #odd rows
holdframe[1:nrow(holdframe) %% 2 == 0,] = secondAllele #even rows

holdframe[holdframe==0]= -9 # replace missing values with -9

groupvec <- NameExtract
for (i in 1:length(unique(NameExtract))) # replace pop names with numbers
{
  groupvec[which(groupvec==unique(NameExtract)[i])] = i
}

holdframe=cbind(rep(NamePops,each=2),rep(groupvec,each=2),rep(NameExtract,each=2),holdframe)
colnames(holdframe)[1:3]=c("ID","PopID","Pop")

popdef <- data.frame(pops=c("P1", "P2", "Admix"), group=c("P1", "P2", "Admixed"))  ##create popdef. First bracket is pop names. Second bracket corresponds to the BGC group names. Order doesn't seem to matter, and this is needed only at pop level (not indiv)


##if(is.character(popdef)){popdef <- utils::read.csv("BGC_popdef.csv",header=F)} #if popdef is a path then read it in 

#Extract the parental data and admixed data
P1_raw <- holdframe[which(holdframe$Pop %in% popdef[which(popdef[,2]=="P1"),1]),]#Parental 1
P2_raw <- holdframe[which(holdframe$Pop %in% popdef[which(popdef[,2]=="P2"),1]),]#Parental 2
P3_raw <- holdframe[which(holdframe$Pop %in% popdef[which(popdef[,2]=="Admixed"),1]),]#Admixed

summary(P1_raw) ##check popID that only the right pop is represented. SampleID can also be checked
# names of the snps
snpnames <- colnames(temp2)  ##names of loci from column1

#map to be used for missing alleles (if any present across loci)
Allele_Map <- data.frame(SNP=snpnames,
                         Allele1=rep(999,length(snpnames)),
                         Allele2=rep(999,length(snpnames))) # 999 is a dummy placeholder

for(i in 1:length(snpnames)){
  #unique alleles for a given snp (locus)
  alleleVals <- as.data.frame(table(as.character(c(P1_raw[,snpnames[i]],P2_raw[,snpnames[i]],P3_raw[,snpnames[i]]))))
  
  # if there is missing data (-9) delete it as a possibe allele
  if(length(which(alleleVals[,1]==(-9)))>0){
    alleleVals <- alleleVals[-which(alleleVals[,1]==(-9)),]
  }
  
  Allele_Map[i,"Allele1"]=as.character(alleleVals[1,1])
  Allele_Map[i,"Allele2"]=as.character(alleleVals[2,1])
}

#NULL vectors
P1_BGC <- NULL
P2_BGC <- NULL

for(i in snpnames){
  # grab vector of alleles and delete replace missing values (-9) with NA
  P1_alleles <- P1_raw[,i];P1_alleles[which(P1_alleles==-9)]=NA
  P2_alleles <- P2_raw[,i];P2_alleles[which(P2_alleles==-9)]=NA
  
  #If the population only has one allele for a given locus then a zero and the allele have be be added
  if(length(table(P1_alleles))==1|sum(is.na(P1_alleles))==length(P1_alleles)){
    if(length(table(P1_alleles))==1){
      hold <- as.data.frame(table(P1_alleles))
      hold[,1] <- as.character(hold[,1])
      hold <- rbind(hold,c(setdiff(as.numeric(Allele_Map[which(Allele_Map$SNP==i),c("Allele1","Allele2")]),hold[1,1]),0)) #add in the extra value
      hold <- hold[order(hold[,1]),] #sort the right order from a conventional table output
      P1_alleles <- hold[,2]
      rm(hold)} else {P1_alleles <- c(0,0)}
  } else {P1_alleles <- as.character(as.data.frame(table(P1_alleles))[,2])}
  
  if(length(table(P2_alleles))==1 | sum(is.na(P2_alleles))==length(P2_alleles)){
    if(length(table(P2_alleles))==1){
      hold <- as.data.frame(table(P2_alleles))
      hold[,1] <- as.character(hold[,1])
      hold <- rbind(hold,c(setdiff(as.numeric(Allele_Map[which(Allele_Map$SNP==i),c("Allele1","Allele2")]),hold[1,1]),0)) #add in the extra value
      hold <- hold[order(hold[,1]),] #sort the right order from a conventional table output
      P2_alleles <- hold[,2]
      rm(hold)} else{P2_alleles <- c(0,0)}
  } else {P2_alleles <- as.character(as.data.frame(table(P2_alleles))[,2])}
  
  
  #for a given locus get the format for BGC
  P1_temp <- c(paste("locus_",i,sep=""),paste(P1_alleles[1],P1_alleles[2],sep=" "))
  P2_temp <- c(paste("locus_",i,sep=""),paste(P2_alleles[1],P2_alleles[2],sep=" "))
  
  #Combine output sequentially for each locus
  P1_BGC <- c(P1_BGC,P1_temp)
  P2_BGC <- c(P2_BGC,P2_temp)
}

path="/Users/alexjvr/2016RADAnalysis/1.2_Phylo/input.files/BGCinput/"  ##specify path
fname="EASTall.new"  ##specify start of name

##Save output for BGC formatted for the parental populations ------------
if(substring(path,nchar(path))!="/"){path=paste0(path,"/")}

utils::write.table(x = P1_BGC,file=paste0(path,fname,"_Parental1_BGC.txt",sep=""),
                   sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)

utils::write.table(x = P2_BGC,file=paste0(path,fname,"_Parental2_BGC.txt",sep=""),
                   sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)


#Convert the admixed data to BGC format --------------

###specify majorminor function
majorminor<- function(Vec, allele_length = 6){
  if(sum(is.na(Vec))!=length(Vec)){
    Vec <- as.character(Vec)
    firstAllele <-  as.data.frame(sapply(Vec,function(x)as.character(substring(x,1,(allele_length/2)))),stringsAsFactors = F)
    secondAllele <-  as.data.frame(sapply(Vec,function(x)as.character(substring(x,(allele_length/2)+1,allele_length))),stringsAsFactors = F)
    x=c(firstAllele[,1],secondAllele[,1])
    AlleleMajor <- as.character(names(which(table(x)==max(table(x)))))
    AlleleMinor <- as.character(names(which(table(x)==min(table(x)))))
    if(length(AlleleMajor)>1){AlleleMajor <- AlleleMajor[1];AlleleMinor=AlleleMinor[2]}
    
    Vec[is.na(Vec)]="-9 -9" #missing data
    Vec=gsub(paste0(AlleleMajor,AlleleMajor),"2 0",Vec) #homozygous major
    Vec=gsub(paste0(AlleleMajor,AlleleMinor),"1 1",Vec) #heterozygous
    Vec=gsub(paste0(AlleleMinor,AlleleMajor),"1 1",Vec) #heterozygous
    Vec=gsub(paste0(AlleleMinor,AlleleMinor),"0 2",Vec) #homozygous minor
  } else {Vec=rep("0 0",length(Vec))}
  return(Vec)
}


#subset data for admixed populations
missingfix<- function(x){ #create functions for apply loop
  hold=x
  hold[grep("000",hold)]=NA
  return(hold)}

#Remove Alleles with missing data and replace with NA
temp3 <- apply(temp2,2,missingfix)

#convert to zygosity format (2 0 - homozygous major, 0 2 - homozygous minor, 1 1 - heterozygous, -9 -9 - missing data )
temp4 <- apply(temp3,2,FUN = function(x) {majorminor(x,allele_length=alleleEx)})

MixedStruct <- temp4[which(NameExtract %in% popdef[which(popdef[,2]=="Admixed"),1]),]
MixedPops <- NameExtract[which(NameExtract %in% popdef[which(popdef[,2]=="Admixed"),1])]

#the number of individuals for all populations but the last (Pop tagged to the end)
PopLengths <- table(MixedPops)[-length(table(MixedPops))]

if(length(table(MixedPops))==2){PopPosition = PopLengths+1}

if(length(table(MixedPops))>2){
  PopPosition <- c(PopLengths[1]+1,rep(NA,(length(PopLengths)-1)))
  for (i in 2:length(PopLengths)){
    PopPosition[i] <- PopLengths[i]+PopPosition[i-1]
  }
}

#Insert the population labels
if(length(table(MixedPops))!=1){
  temp5 <- apply(MixedStruct,2,function(x){insert_vals(x,breaks=PopPosition,
                                                       newVal=paste0("pop_",unique(MixedPops)[2:length(unique(MixedPops))]))})} else {
                                                         temp5 <- MixedStruct}

temp5=as.data.frame(temp5,stringsAsFactors = FALSE)

#Add the "locus_" and first "pop_" labels
temp6=as.matrix(rbind(paste0("locus_",colnames(temp5)),
                      rep(paste0("pop_",unique(MixedPops)[1]),length(temp5)),
                      temp5))

#redim as a single vector
MixedData=as.vector(temp6)

##Save output for BGC formatted for the parental and mixed populations ------------
utils::write.table(x = MixedData,file=paste(path,fname,"_Admixed_BGC.txt",sep=""),
                   sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)

} #end of function
