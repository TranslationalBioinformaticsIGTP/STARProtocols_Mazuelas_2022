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

analysis.dir <- params$analysis.dir
sample.data <- read.table(file = params$sample.data.file , header = T, sep = "\t", stringsAsFactors = FALSE, comment.char = "")
file.names <- sample.data$File.Name #Names of the Fastq file samples: this must be contained in the data.frame sample.data information.

# Salmon alignement and quantification parameters
salmonDir <- params$salmon.software #path where salmon is installed
file1.suffix <- params$file1.suffix
file2.suffix <- params$file2.suffix
fastqdir <- params$fastqdir  # Fastq data dir
transcript.index <- params$transcript.index # Index previusly generated
output.suffix <- params$output.suffix
output.quants <- params$output.quants # Directory to save Salmon quant results
threads <- params$threads # Number of threads

###################   Functions   ################### 
source(file ="./utils.R")

### salmonAlignment
salmonAlignment <- function(sample.name, salmonDir,
                            file1.suffix, file2.suffix,
                            fastqdir,
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
    full.command <- paste0(salmonDir, " quant -i",
                           transcript.index," -l ",
                           libtype, " -1 ", 
                           fastqdir,
                           sample.name,
                           file1.suffix, " -2 ", 
                           fastqdir, sample.name, 
                           file2.suffix, " --validateMappings -p ",
                           as.character(threads),
                           " -o ", 
                           output.quants,
                           sample.name, output.suffix )
    system(full.command, wait = TRUE)
  }
}


###################                     Salmon alignment                      #####################
###################################################################################################
############            ATENTION! BEFORE PseudoAlign the samples           ########################
############  YOU MUST generate your salmon index using STARProtocols_SalmonIndex.sh  #############
###################################################################################################

# ## Executing Salmon alignement by Selective alignement and quantification
for(i in seq_len(length(file.names))){
  sn <- file.names[i]
  salmonAlignment(sample.name = sn, salmonDir,
                  file1.suffix = file1.suffix,
                  file2.suffix = file2.suffix,
                  fastqdir = fastqdir,
                  transcript.index = transcript.index,
                  output.suffix = output.suffix,
                  output.quants = output.quants,
                  threads = threads
  )
}
now.msg("   Salmon done")
