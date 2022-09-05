#!/bin/bash

################ Parameters #####################
#where your index file will be stored
index_dir= "./SalmonIndexing_Metadata" 
#Reference genome
refGenome="./references/hg38.fa.gz"
#Reference transcriptome
refTranscriptome="./references/refMrna.fa.gz"
#Salmon software path
Salmon= "salmon"
#Threads used to run
threads= "4"
######### Get decoys from reference genome ########### 
# For obtaining the decoys.txt file we followed the process described at https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/
# Reference genome and transcriptome used for the analaysis is detailed at Mazuelas et al. 2022 

echo "***********************************************************"
echo "*** `date` Getting Salmon index  ****"
echo "***********************************************************"

mkdir -p $index_dir
cd $index_dir


grep "^>" <(gunzip -c $refGenome) | cut -d " " -f 1 > decoys.txt
sed -i.bak -e 's/>//g' decoys.txt

# concatenating transcriptome and genome reference file for index
zcat $refTransciptome $refGenome > gentrome.fa.gz

#Get Salmon index
$Samlon index -t ./gentrome.fa.gz -d ./decoys.txt -p $threads -i salmon_index -k 31

echo "*** `date` Salmon index obtained  ****"

