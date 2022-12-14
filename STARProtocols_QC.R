##########################################################
#           STAR Protocols Mazuelas et al. 2022          #
##########################################################
####### Packages Needed
message("Loading libraries...")
if(!suppressPackageStartupMessages(require("BiocManager", quietly = TRUE))) install.packages("BiocManager")


if(!suppressPackageStartupMessages(require("DESeq2", quietly = TRUE))) BiocManager::install("DESeq2")
suppressPackageStartupMessages(library(DESeq2))
if(!suppressPackageStartupMessages(require("tximport", quietly = TRUE))) BiocManager::install("tximport")
suppressPackageStartupMessages(library(tximport))
if(!suppressPackageStartupMessages(require("org.Hs.eg.db", quietly = TRUE))) BiocManager::install("org.Hs.eg.db")
suppressPackageStartupMessages(library(org.Hs.eg.db))
if(!suppressPackageStartupMessages(require("yaml", quietly = TRUE))) utils::install.packages("yaml")
suppressPackageStartupMessages(library(yaml))
if(!suppressPackageStartupMessages(require("ggplot2", quietly = TRUE))) utils::install.packages("ggplot2")
suppressPackageStartupMessages(library(ggplot2))
if(!suppressPackageStartupMessages(require("ggbeeswarm", quietly = TRUE))) utils::install.packages("ggbeeswarm")
suppressPackageStartupMessages(library(ggbeeswarm))
if(!suppressPackageStartupMessages(require("pheatmap", quietly = TRUE))) utils::install.packages("pheatmap")
suppressPackageStartupMessages(library(pheatmap))
if(!suppressPackageStartupMessages(require("RColorBrewer", quietly = TRUE))) utils::install.packages("RColorBrewer")
suppressPackageStartupMessages(library(RColorBrewer))
if(!suppressPackageStartupMessages(require("gplots", quietly = TRUE))) utils::install.packages("gplots")
suppressPackageStartupMessages(library(gplots))

#######################  Loading Functions #########################
source(file = "./rna_seq_Functions.R")

###### Parameters ######

# QC parameters:Modify this parameter file before starting with the QC analysis.
params <- yaml.load_file("./Parameters/QCheatmap_pipeline_parameters.yaml")

#Sample data fastq files
fastq.dir <- params$fastq.dir  # Fastq data dir
fastq.files <- list.files(fastq.dir)

# Getting sample names from Fastq data files
file.names <- unique(gsub("_[0-9]\\.[^.]+\\.[^.]+","",fastq.files))  #regex to delete i.e. "_1.fastq.gz"
sample.names <- unique(file.names) 

#sample data information
sample.data <- data.frame(Sample.Name = sample.names, sample.group = rep_len(x = c(1,2), length.out = length(sample.names)))#Adding sample.group column

# Loading the Mazuelas et al. 2022 sample count expression of roadmap markers (QCroadmap counts).
counts.roadmap <- as.matrix(read.table(params$counts.roadmap, header = T, sep = "\t"))


# Tximport paramenters
output.quants <- params$salmon.dir
orgdb <- org.Hs.eg.db
org.columns <- params$org.columns
org.keytype <- params$org.keytype

# DESeq2 parameters
filt.min.reads <- params$filt.min.reads
filt.min.samples <-params$filt.min.samples
pvalue <- params$pvalue

#Color heatmap
color.plate <-  bluered(80)

# Genes Roadmap pipeline
gene.markers <- c() 
for(i in seq_len(length(params$stages))){
  gr <- params$stages[i]
  mks<- read.table(file = file.path(params$roadmap.markers.dir, paste0("up_", gr,".txt")), header = FALSE,stringsAsFactors = FALSE)[,1]
  # mks <- markers[[gr]]
  names(mks) <- rep(gr, length(mks))
  gene.markers <- c(gene.markers,mks)
  
}

###################################################################################################
############            ATENTION! BEFORE STARTING THE QC PROTOCOL          ########################
############  YOU MUST PROCESS YOUR DATA USING STAR PROTOCOLS SALMON ALIGNMENT FILE  ##############
###################################################################################################


##################              Tximport              #####################
salmonquants.fl <-file.path(output.quants, paste0(file.names,"_quant"), "quant.sf")
names(salmonquants.fl)<- file.names
txi.salmon <- importQuantsData(quant.files = salmonquants.fl, orgdb = orgdb,
                               orgdb.keytype = orgdb.keytype, org.columns = org.columns,
                               tximpot.type = "salmon")
#Change Sample name
colnames(txi.salmon$abundance) <- sample.names
colnames(txi.salmon$counts) <- sample.names
colnames(txi.salmon$length) <- sample.names


#####################     Procesing the data for the heatmap QC sphere formation    ##############################
# we recommend filtering your data using the same parametes as we used for our data (filt.min.reads = 5; filt.min.samples = 1)
filtered.dds <- suppressWarnings(getFilteredDDS(tximport = txi.salmon,
                                                samples_group = "sample.group", 
                                                samples_df = sample.data, 
                                                filter_min_reads = filt.min.reads, 
                                                filter_min_samples = filt.min.samples))

counts.deseq <- counts(filtered.dds)[rownames(filtered.dds) %in% rownames(counts.roadmap),]


# Join your processed data to QC counts.roadmap to compare your differentiation expression results
counts.deseq <- cbind(counts.roadmap[rownames(counts.roadmap)%in% rownames(counts.deseq),], counts.deseq)
 
#Perform a rlog transformation of the joined counts.
dds.rlog <- rlog(counts.deseq)

####### Data counts to plot #####
gene.markers <- gene.markers[gene.markers %in% rownames(dds.rlog)]
count.data <- dds.rlog[gene.markers,] # sorting rows as FiPS WT 2D markers

# Data normalization to plot in the heatmap
data_subset_norm <- t(apply(count.data, 1, cal_z_score))

#without iPSC
data_subset_norm <- data_subset_norm[,!grepl(x=colnames(data_subset_norm),pattern = "PSC")]

# Adjusting the zscore color associated with the different heatmaps
paletteLength <-length(color.plate)
myBreaks <- c(seq(min(data_subset_norm), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(data_subset_norm)/paletteLength, max(data_subset_norm), length.out=floor(paletteLength/2)))


#heatmap of QC Spheres stage specific markers
if(!file.exists(params$heatmap.dir)) dir.create(params$heatmap.dir)

png(filename = file.path(params$heatmap.dir,"QC_NeurospheresFormation_RoadmapExpression.png"),width=1000, height= 800)
pheatmap(data_subset_norm,
         cluster_cols = F,
         cluster_rows =F,
         fontsize = 20,
         color = bluered(80),
         fontsize_col = 25,
         show_rownames = FALSE,
         legend = TRUE,
         legend_breaks = c(round(min(data_subset_norm)),round(max(data_subset_norm))),
         breaks =myBreaks ,
         legend_labels = c("Down","Up"),
         annotation_names_row = FALSE,
         annotation_names_col = FALSE,
         lwd =3,
         treeheight_row = 40,
         treeheight_col = 40,
         angle_col = 90, 
         margins=c(50,50,50,50)

)
tt <- dev.off()

now.msg("QC Heatmap of neurofibromaspheres obtained")
