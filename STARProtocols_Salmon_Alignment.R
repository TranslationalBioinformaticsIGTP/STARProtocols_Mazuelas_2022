#######################################################################
#                        Salmon Alignment                             #
#######################################################################

###################################################################################################
############            ATENTION! BEFORE PseudoAlign the samples           ########################
############  YOU MUST generate your salmon index using STARProtocols_SalmonIndex.sh  #############
###################################################################################################

####### Packages Needed
if(!require("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2")
library(DESeq2)
if(!require("org.Hs.eg.db", quietly = TRUE)) BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
if(!require("yaml", quietly = TRUE)) utils::install.packages("yaml")
library(yaml)

####### Parameters

#All samples parameters
params <- yaml.load_file("./Parameters/QCheatmap_pipeline_parameters.yaml")

salmon.dir <- params$salmon.dir
results.dir <- "results"
if(!file.exists(results.dir)) dir.create(results.dir)
if(!file.exists(salmon.dir)) dir.create(salmon.dir)


#Sample data fastq files
fastq.dir <- params$fastq.dir  # Fastq data dir
fastq.files <- list.files(fastq.dir)

# Getting sample names from data files
file.names <- gsub("_[0-9]\\.[^.]+\\.[^.]+","",fastq.files)  #regex to delete i.e. "_1.fastq.gz"
file.names <- unique(file.names)

# Salmon alignement and quantification parameters
salmon.soft.dir <- params$salmon.soft.dir #path where salmon is installed
transcript.index <- params$transcript.index # Index previusly generated
output.suffix <- params$output.suffix
output.quants <- params$output.quants # Directory to save Salmon quant results
threads <- params$threads # Number of threads

###################   Functions   ################### 
source(file ="./utils.R")

### salmonAlignment
salmonAlignment <- function(sample.name, salmon.soft.dir,
                            file1.suffix, file2.suffix,
                            fastq.dir,
                            transcript.index,
                            output.suffix,
                            output.quants, 
                            libtype = "IU",
                            threads = 4,
                            verbose = TRUE){
  now.msg("Salmon starting...")
  file.name <- file.path(output.quants, paste0(sample.name, output.suffix))
  
  if (!file.exists(file.name)){
    message("missing file", file.name)
    full.command <- paste0(salmon.soft.dir, " quant -i",
                           transcript.index," -l ",
                           libtype, " -1 ", 
                           file.path(fastq.dir,
                           paste0(sample.name,
                           file1.suffix)), " -2 ", 
                           file.path(fastq.dir,
                                     paste0(sample.name,
                                     file2.suffix)), " --validateMappings -p ",
                           as.character(threads),
                           " -o ", 
                           file.path(output.quants,
                           paste0(sample.name, output.suffix)))
    system(full.command, wait = TRUE)
  }
}


###################                     Salmon alignment                      #####################
###################################################################################################
############            ATENTION! BEFORE PseudoAlign the samples           ########################
############  YOU MUST generate your salmon index using STARProtocols_SalmonIndex.sh  #############
###################################################################################################

# ## Executing Salmon alignement by Selective alignement and quantification
i=1
for(i in seq_len(length(file.names))){
  sn <- file.names[i]
  sn.files <- fastq.files[grepl(sn, fastq.files)]
  file.suffix <- gsub(sn,"", sn.files)
  salmonAlignment(sample.name = sn, salmon.soft.dir = salmon.soft.dir,
                  file1.suffix = file.suffix[1],
                  file2.suffix = file.suffix[2],
                  fastq.dir = fastq.dir,
                  transcript.index = transcript.index,
                  output.suffix = output.suffix,
                  output.quants = salmon.dir,
                  threads = threads
  )
}
now.msg("   Salmon done")
