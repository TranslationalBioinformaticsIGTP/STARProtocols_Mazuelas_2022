######################################## README ######################################

Pipeline for performing QC analysis of neurofibromaspheres differentiation.
This analysis is performed in a Linux system

 This directory contains the necesary files to perform a QC of a neurospheres differentiation
 RNA-seq data necesary to perfome this QC is available at EGA under the accession number EGAS00001005907.
 Figure 3E published at Mazuelas et al. 2022 (https://doi.org/10.1016/j.celrep.2022.110385) is a representation of the purpouse of this QC.

 Before starting the QC analysis, it is necesary to follow the next steps:
 #1- Download R version 4.2.0 https://cran.r-project.org/ and R-Studio from https://www.rstudio.com/products/rstudio/download/#download 
 #2- Download Bioconductor v3.15 or supperior https://www.bioconductor.org/install/
 #3- Download Salmon v.1.9.0 from https://salmon.readthedocs.io/en/latest/index.html (salmon is only available for Linux systems)
 #4- Download the reference genome hg38.fa.gz from UCSC https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/?C=M;O=A
 #5- Download the reference transcriptome refMrna.fa.gz from UCSC https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/?C=M;O=A
 #6- Fill up the QCheatmap_pipeline_parameters.yaml with the necesary information. There are some parameters you should fill up befor starting the pipeline.
 #7- Create a Salmon index using STARProtocols_SalmonIndex.sh contained in this directory for processing your data
 #8- Pseudoalign your fastq files using  STARProtocols_Salmon_Alignment.R
 #9- Proceed to create a QC heatmap plot to see if your neurofibromaspheres are well differentiated. 
 *** Please, refer to Mazuelas et al. 2022 Figure 3 E(https://doi.org/10.1016/j.celrep.2022.110385) for a reference of a good neurofibromasphere differenciation.