# README 

Pipeline for performing a Quality Control (QC) analysis of neurofibromaspheres differentiation.

**WARNING!** This analysis is performed in a Linux system

This directory contains the necessary files to perform a QC of a new generation of neurospheres differentiation. There are diffent files for processing your RNA-seq data and compare the expression of your
neurofibromaspheres with the expression of the differentiantion publised at Mazuelas et al. 2022 (https://doi.org/10.1016/j.celrep.2022.110385). Figure 3E of this paper is a representation of how the differentiation of your spheres should look like.

In the case you have enough space in your computer and computation capacity, you can also the raw data of Mazuelas et al. 2022 available at EGA under the acession number EGAS00001005907, but it is not strickly necesary for performing this QC analysis.



For performing a good QC analysis, it is necessary to follow carefully the next steps:

 
##1- Download R version 4.2.0 https://cran.r-project.org/ and R-Studio from https://www.rstudio.com/products/rstudio/download/#download 

##2- Download Bioconductor v3.15 or superior https://www.bioconductor.org/install/

##3- Download Salmon v.1.9.0 from https://salmon.readthedocs.io/en/latest/index.html **(salmon is only available for Linux systems)**

##4- Download the reference genome hg38.fa.gz from UCSC https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/?C=M;O=A

##5- Download the reference transcriptome refMrna.fa.gz from UCSC https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/?C=M;O=A

##6- Fill up the QCheatmap_pipeline_parameters.yaml with the necesary information. There are some parameters you should fill up befor starting the pipeline.

##7- Create a Salmon index using STARProtocols_SalmonIndex.sh contained in this directory for processing your data

##8- Pseudoalign your fastq files using  Rscript STARProtocols_Salmon_Alignment.R in your command line.

##9- Proceed to create a QC heatmap plot to see if your neurofibromaspheres are well differentiated runnig Rscript STARProtocols_QC.R in your command line. 

*Please, refer to Mazuelas et al. 2022 Figure 3 E(https://doi.org/10.1016/j.celrep.2022.110385) for a reference of a good neurofibromasphere differenciation.*

