############################ Functions RNA-Seq #########################
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

####importQuantData
importQuantsData <- function(quant.files, orgdb, orgdb.keytype, org.columns, tximpot.type= "salmon", verbose = TRUE, ...){
  now.msg("Tximport starting...", verbose = verbose)
  #translating transctipts to gene symbol
  k <- keys(orgdb, keytype = org.keytype)
  tx2gene <- AnnotationDbi::select(x = orgdb,
                                   keys = k,
                                   columns = org.columns,
                                   keytype = org.keytype)
  
  
  #importing qunats data
  txi.salmon <- tximport(files = salmonquants.fl,
                         type = "salmon",
                         tx2gene = tx2gene, ...)
  now.msg("  Tximport done", verbose = verbose)
  return(txi.salmon)
}

###selectDataFromTximport  

selectDataFromTximport <- function(tximport, sample_names){
  tximport$abundance <- data.matrix(data.frame(tximport$abundance)[colnames(tximport$abundance) %in% sample_names])
  colnames(tximport$abundance)<- sample_names
  tximport$counts <-data.matrix(data.frame(tximport$counts)[colnames(tximport$counts) %in% sample_names])
  colnames(tximport$counts)<- sample_names
  tximport$length <- data.matrix(data.frame(tximport$length)[colnames(tximport$length) %in% sample_names])
  colnames(tximport$length)<- sample_names
  return(tximport)
}



#' getFilteredDDS
#'
#' @description
#'
#' getFilteredDDS function gets a filtered DESeqDataSet object where the design formula is added depending on
#' the test performed.
#'
#' @details
#'
#' getFilteredDDS function gets a filtered DESeqDataSet object. This funtion contains DESeqDataSetFromTximport function
#' from DESeq2 package, used to store the input values from tximport list. Depending on the test the user wants to do,
#' the design formula of DESeqDataSetFromTximport will vary. If the user like studying more than one sample group, LRT test
#' will be performed. If the interaction between different sample groups is analysed, then this function will join in
#' a unique sample group variable as one condition to study. If the user like to study one condition (sample_group parameter), it will
#' be obtained a simple formula design using any test.
#'
#' @usage
#'
#' getFilteredDDS(tximport, samples_group, samples_df, filter_min_reads = 5,
#' filter_min_samples  = 1, test = "Wald", interaction = FALSE, verbose = TRUE)
#'
#' @param tximport (list) tximport object produced by tximport function
#' @param samples_group (character) name of the factor from samples_df to contrast their levels. If sample_group is
#' more than one, LRT  test must be performed.
#' @param samples_df (data.frame) data frame containing all information about the samples.
#' @param filter_min_reads (numeric) minimum of count reads to filter (default: filter_min_reads = 5)
#' @param filter_min_samples (numeric)minimun of columns to apply filter_min_reads(default: filter_min_samples = 1).
#' @param test (character) test that is going to be performed in order to get DESeq design formula.
#' it could be whether "Wald" or "LRT" tests (default: test = "Wald").
#' @param interaction (logical) only used if LRT test is performed and sample_group are two.
#' Whether the design formula study the interaction of sample_group. (default: ineraction = FALSE)
#' @param verbose (logical) whether the user like some messages over the process.
#'
#'
#' @return
#'
#' getFilteredDDS's main output is a filtered DESeqDataSet object ready to get DEG.
#'
#' @export getFilteredDDS
#'
#'
getFilteredDDS <- function(tximport, samples_group, samples_df, filter_min_reads = 5,
                           filter_min_samples  = 1, test = "Wald", interaction = FALSE,
                           verbose = TRUE){
  
  now.msg(" Filtering DESeqDataSet...", verbose = verbose)
  
  if(test == "Wald"){
    
    if(length(samples_group) == 1){
      
      # Generation of a factor to contrast
      cont<-samples_df[,samples_group]
      df<- cbind(samples_df, design = as.character(cont), stringsAsFactors=FALSE)
      
      #DESFromTXI
      dds <- DESeqDataSetFromTximport(txi = tximport, colData = df, design = ~ design)
      
      
    }else{
      stop("Wald test is only for one conditon and sample_group's parameter has more than one")
    }
  }
  
  if(test == "LRT"){
    
    msg(" generating LRT test design for contrast  ...", verbose = verbose)
    
    # Generation of a factor to contrast
    
    if (samples_group > 1 && interaction == FALSE){
     
       # Generation of a factor to contrast
      
      cont<-samples_df[,samples_group]
      df<- cbind(samples_df, design = cont, stringsAsFactors=FALSE)
      
      #DESFromTXI
      
      for( i in 1:(length(samples_group)-1)){
        dds <- DESeqDataSetFromTximport(txi = tximport, colData = df,
                                        design =formula(paste0("~","design.",samples_group[i],
                                                               "+", "design.", samples_group[i+1])))
      }
      
    }else if ( samples_group == 1 && interaction == FALSE){
      # Generation of a factor to contrast
      
      cont<-samples_df[,samples_group]
      df<- cbind(samples_df, design = as.character(cont), stringsAsFactors=FALSE)
        
      #DESFromTXI
      dds <- DESeqDataSetFromTximport(txi = tximport, colData = df, design = ~ design)
        
      
      
    }else if (samples_group == 1 && interaction == TRUE){
      
      msg(" Preparing LRT design for interaction contrast...", verbose = verbose)
      
      for( i in 1:(length(samples_group)-1)){
      df$design = factor(paste0(df[,paste0("design.",samples_group[i])],"_",
                                df[,paste0("design.",samples_group[i+1])]))
      }
      
      dds <- DESeqDataSetFromTximport(txi = tximport, colData = df, 
                                      design = ~design)
    }
  }
  
  #estimateSizeFactor
  dds <- estimateSizeFactors(dds)
  
  # Filtering data
  msg("Fitering data...", verbose = verbose)
  keep <- rowSums(counts(dds) >= filter_min_reads) >= filter_min_samples
  filtered_dds <- dds[keep,]
  
  now.msg(" DESeqDataSet Filtered", verbose = verbose)
  
  return(filtered_dds)
  
}





#we create a function to calculate z-score
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

#Functions to wrap labels of barplot 
## Core wrapping function
wrap.it <- function(x, len)
{ 
  sapply(x, function(y) paste(strwrap(y, len), 
                              collapse = "\n"), 
         USE.NAMES = FALSE)
}


## Call this function with a list or vector
wrap.labels <- function(x, len)
{
  if (is.list(x))
  {
    lapply(x, wrap.it, len)
  } else {
    wrap.it(x, len)
  }
}
