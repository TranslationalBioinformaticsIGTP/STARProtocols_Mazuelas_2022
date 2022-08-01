#!/bin/bash

################ Parameters #####################
#where your index file will be stored
index_dir= "./SalmonIndexing_Metadata" 
#Reference genome 
refGenome="./Genomes/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
#Reference transcriptome
refTranscriptome="./Transcriptomes/hg38_refMrna.fa.gz"
#Salmon software
Salmon= "/software/debian-10/bio/salmon-1.8.0/bin/salmon"

######### Get decoys from reference genome ########### 
# For obtaining the decoys.txt file we followed the process described at https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/
# Reference genome and transcriptome used for the analaysis is detailed at Mazuelas et al. 2022 

echo "***********************************************************"
echo "*** `date` Getting Salmon index  ****"
echo "***********************************************************"

cd $index_dir


grep "^>" <(gunzip -c $refGenome) | cut -d " " -f 1 > decoys.txt
sed -i.bak -e 's/>//g' decoys.txt

# concatenating transcriptome and genome reference file for index
zcat $refTransciptome $refGenome > gentrome.fa.gz

#Get Salmon index
$Samlon index -t ./gentrome.fa.gz -d ./decoys.txt -p 4 -i salmon_index -k 31

echo "*** `date` Salmon index obtained  ****"

