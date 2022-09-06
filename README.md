# README

Pipeline for performing a Quality Control (QC) analysis of neurofibromaspheres differentiation.

**IMPORTANT!** To run this analysis you need a Linux system

This directory contains the necessary files to perform QC of a new generation of neurofibromaspheres. There are different files for processing your RNA-seq data and compare the expression of your neurofibromaspheres with the expression of the differentiation published at Mazuelas et al. 2022 (<https://doi.org/10.1016/j.celrep.2022.110385>). Figure 3E in this paper is a representation of how the differentiation of your spheres should look like.

In the case you have enough space in your computer and computational capacity, you can also download the raw data of Mazuelas et al. 2022 available at EGA under the accession number EGAS00001005907, but it is not strictly necessary for performing this QC analysis.

To perform the QC analysis, please follow these steps:

#### Install Prerequisites

1.  Download and install R version 4.2.0 <https://cran.r-project.org/>.

2.  **OPTIONAL** Download and install Rstudio from <https://www.rstudio.com/products/rstudio/download/#download>. This will provide a friendly environment to run R code and generate figures.

3.  Install Bioconductor v3.15 or later <https://www.bioconductor.org/install/>

4.  Download and install Salmon v.1.9.0. Follow the instruction at <https://github.com/COMBINE-lab/salmon/releases/tag/v1.9.0> **(salmon is only available for Linux systems)**

5. Download the necesary code release v0.0.1 from XXX
 
#### Download Reference data and prepare the Salmon index

6.  Download the reference genome `hg38.fa.gz` from UCSC:

```
mkdir -p references && wget 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz' -O references/hg38.fa.gz
```

7.  Download the reference transcriptome `refMrna.fa.gz` from UCSC:
```
mkdir -p references && wget 'ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/refMrna.fa.gz' -O references/refMrna.fa.gz
```
8.  Create a Salmon index using `STARProtocols_SalmonIndex.sh` contained in this directory for processing your data
```
. ./STARProtocols_SalmonIndex.sh
```

#### Run QC pipeline

9.  **WARNING!** Add at the `QCheatmap_pipeline_parameters.yaml` the directory location where your fastq files are stored.

10.  Pseudoalign your fastq files using `STARProtocols_Salmon_Alignment.R` in your command line.
```
Rscript STARProtocols_Salmon_Alignment.R
```
11.  Please fill up `Sample.Info.Guide.csv` with your own sample information (fastq file name and sample name).

12.  Proceed to create a QC heatmap plot to see if your neurofibromaspheres are well differentiated running `STARProtocols_QC.R` in your command line.
```
Rscript STARProtocols_QC.R
```


*Please, refer to Figure 3E in Mazuelas et al. 2022 (<https://doi.org/10.1016/j.celrep.2022.110385>) for an example of a good neurofibromasphere differentiation.*
