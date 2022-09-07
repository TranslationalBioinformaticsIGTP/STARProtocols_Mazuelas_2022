#!/bin/bash

################ Parameters #####################
#where your index file will be stored
index_dir="./SalmonIndex" 
#Reference genome
refGenome="./references/hg38.fa.gz"
#Reference transcriptome
refTranscriptome="./references/refMrna.fa.gz"
#Salmon software path
Salmon="salmon"
#Threads used to run
threads="1"
######### Get decoys from reference genome ########### 
# For obtaining the decoys.txt file we followed the process described at https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/
# Reference genome and transcriptome used for the analaysis is detailed at Mazuelas et al. 2022 

echo "***********************************************************"
echo "*** `date` Getting Salmon index  ****"
echo "***********************************************************"

mkdir -p $index_dir


grep "^>" <(gunzip -c $refGenome) | cut -d " " -f 1 > $index_dir/decoys.txt
sed -i.bak -e 's/>//g' $index_dir/decoys.txt

# concatenating transcriptome and genome reference file for index
zcat $refTranscriptome $refGenome > $index_dir/gentrome.fa.gz

#Get Salmon index
$Salmon index -t $index_dir/gentrome.fa.gz -d $index_dir/decoys.txt -p $threads -i $index_dir/salmon_index 

echo "*** `date` Salmon index obtained  ****"

