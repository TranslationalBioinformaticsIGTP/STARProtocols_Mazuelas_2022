##########################################################
#           STAR Protocols Mazuelas et al. 2022          #
##########################################################
####### Packages Needed
library(DESeq2)
library(tximport)
library(org.Hs.eg.db)
library(yaml)
library(ggplot2)
library(ggbeeswarm)
library(pheatmap)
library(viridis)
library(apeglm)
library(ggrepel)
library(RColorBrewer)
library(gplots)

#######################  Loading Functions #########################
source(file = "./rna_seq_Functions.R")

###### Parameters ######

# QC parameters:Modify this parameter file before starting with the QC analysis.
params <- yaml.load_file("./Parameters/QCheatmap_pipeline_parameters.yaml")

# Loading your sample data information 
sample.data <- read.table(file = params$sample.data.file , header = T, sep = ",", stringsAsFactors = FALSE, comment.char = "")

# Loading the Mazuelas et al. 2022 sample count expression of roadmap markers (QCroadmap counts).
counts.roadmap <- as.matrix(read.table(params$counts.roadmap, header = T, sep = "\t"))

# Getting the file names of the analysis: this must be contained in the data.frame sample.data information.
file.names <- sample.data$File.Name #names of the fastq files
sample.names <- sample.data$Sample.Name #names of the samples

# Tximport paramenters
output.quants <- params$output.quants
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
for(i in seq_len(length(stages))){
  gr <- stages[i]
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
filtered.dds <- getFilteredDDS(tximport = txi.salmon,
                               samples_group = samples.group, 
                               samples_df = sample.data, 
                               filter_min_reads = filt.min.reads, 
                               filter_min_samples = filt.min.samples)

counts.deseq <- counts(filtered.dds)[rownames(filtered.dds) %in% gene.markers,]


# Join your processed data to QC counts.roadmap to compare your differentiation expression results
counts.deseq <- cbind(counts.roadmap,counts.deseq)

#Perform a rlog transformation of the joined counts.
dds.rlog <- rlog(counts.deseq)

####### Data counts to plot #####
gene.markers <- gene.markers[gene.markers %in% rownames(dds.rlog)]
count.data <- dds.rlog[gene.markers,] # sorting rows as FiPS WT 2D markers
colnames(count.data)

# Data normalization to plot in the heatmap
data_subset_norm <- t(apply(count.data, 1, cal_z_score))
colnames(data_subset_norm)

#without iPSC
data_subset_norm <- data_subset_norm[,!grepl(x=colnames(data_subset_norm),pattern = "PSC")]

# Adjusting the zscore color associated with the different heatmaps
paletteLength <-length(color.plate)
myBreaks <- c(seq(min(data_subset_norm), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(data_subset_norm)/paletteLength, max(data_subset_norm), length.out=floor(paletteLength/2)))


#heatmap of QC Spheres stage specific markers
png(filename = file.path(params$heatmap.dir,"QC_NeurospheresFormation_RoadmapExpression.png"),width=1000, heigth= 800)
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
dev.off()
