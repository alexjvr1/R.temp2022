# raw ddRAD data to SNP calls

Data are already demultiplexed, trimmed for adapter, and where a sample was sequenced on multiple lanes the sequences were concatenated (xxcat.fq.trim.gz)

All analyses conducted here: /newhome/bzzjrb/R.temp/

Draft Ref genome from [Sanger](ftp://ngs.sanger.ac.uk/scratch/project/grit/VGP/aRanTem1)


## Map to the genome

1. Index the genome
```
/newhome/bzzjrb/R.temp/RefGenome

#Index the reference genome if needed. Check if the *fasta.fai* file exists in the SpeciesName/RefGenome/ folder in your local directory. If not, run the indexing code. 

#index reference genome
module load apps/bwa-0.7.15
bwa index RefGenome/*fa

```

```
#Create files with input names
## CH
ls CH1027/*gz >> CH.names
sed -i s:CH1027/::g CH.names


## SE
#We're pointing to two input folders so I'll leave the path in the sample names folder
ls SE193/*gz >> SE.names

sed -i 's:SE193/::g' *names
sed -i 's:01d_musAll_merged/::g' *names

#make output directories. 
mkdir 02a_CH_mapped
mkdir 02a_SE_mapped


#Check that you're pointing to the correct reference genome

#Check that the file separator makes sense: 
##sample_name=`echo ${NAME1} | awk -F "_R" '{print $1}'`
#Change the -F "xxx" according to the file names. 
#e.g the above works for files named as follows: 
#HS-01-2016-26_L007_cutadapt_filtered_R2.fastq.gz
#we want only the first part of this name to carry through. 
```

## 
