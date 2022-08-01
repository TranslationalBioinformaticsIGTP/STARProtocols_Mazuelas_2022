######################################## README ######################################

Pipeline for performing QC analysis of neurofibromaspheres differentiation

 This directory contains the necesary files to perform a QC of a neurospheres differentiation
 RNA-seq data necesary to perfome this QC is available at EGA under the accession number EGAS00001005907.
 Figure 3E published at Mazuelas et al. 2022 (https://doi.org/10.1016/j.celrep.2022.110385) is a representation of the purpouse of this QC.

 Before starting the QC analysis, it is necesary to follow the next steps:
 1- Create a Salmon index using STARProtocols_SalmonIndex.sh contained in this directory
 2- Pseudoalign your fastq files using  STARProtocols_Salmon_Alignment.R
 3- Proceed to create a QC heatmap plot to see if your neurofibromaspheres are well differentiated. 
 *** Please, refer to Mazuelas et al. 2022 (https://doi.org/10.1016/j.celrep.2022.110385) for a reference of a good neurofibromasphere differenciation.