################### CUT&TAG DATA ANALYSIS MODULE #################

# install.packages("rlist")
# install.packages("DiffBind")
# install.packages("ChIPQC")
# if (!requireNamespace("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager",dependencies=TRUE)
# }

# BiocManager::install()
# BiocManager::install("Matrix", force=TRUE,dependencies=TRUE)
# BiocManager::install("ChIPQC", force=TRUE, dependencies=TRUE)

# ChIPQC need following code to work
# https://www.biostars.org/p/357154

library("BiocParallel")
register(DoparParam())
registered()
bpparam("SerialParam")

library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library(ChIPQC)
library(DiffBind)
library(rlist)

# set and go to outputs folder
setwd("C:\\Users\\Jack Fan\\Documents\\R\\CUT_Tag\\outputs")


# read in samples
test <- dba(sampleSheet='C:\\Users\\Jack Fan\\Documents\\R\\CUT_Tag\\inputs\\subsamples.txt')
test

pdf("raw_heatmapPlot.pdf")
# png("raw_heatmapPlot.png")
plot(test)
dev.off()

test_qc = ChIPQC(test,annotation="hg38")

# calculate a binding matrix with scores based on read counts for every sample (affinity scores), 
# rather than confidence scores for only those peaks called in a specific sample (occupancy scores).
# take some time, please save the result to save time
test_counts <- dba.count(test)
test_counts
list.save(test_counts,'test.rds')

# test_counts <- readRDS('test.rds')

# FRiP, which stands for Fraction of Reads in Peaks, indicates which samples show more enrichment overall
info <- dba.show(test_counts)
libsizes <- data.frame(LibReads=info$Reads, FRiP=info$FRiP, PeakReads=round(info$Reads * info$FRiP))
rownames(libsizes) <- info$ID
libsizes

pdf("PCAplot_counted.pdf")
dba.plotPCA(test_counts, DBA_CONDITION, label=DBA_REPLICATE)
dev.off()

# plot a new correlation heatmap based on the count scores [make more sense!]
#pdf("heatmap_counted.pdf")
#png("heatmap_counted.png")
plot(test_counts)
# Normalizing the data
test_counts_norm <- dba.normalize(test_counts)
norm <- dba.normalize(test_counts_norm, bRetrieve=TRUE)
normlibs <- data.frame(FullLibSize=norm$lib.sizes, NormFacs=norm$norm.factors, NormLibSize=round(norm$lib.sizes/norm$norm.factors))
rownames(normlibs) <- info$ID
normlibs

# Establishing contrast under a condition
test_counts_norm <- dba.contrast(test_counts_norm, categories=DBA_CONDITION,minMembers = 2)

# Show all contrast with the condition
dba.show(test_counts_norm, bContrasts=TRUE)

# Performing the differential enrichment analysis
test_counts_norm <- dba.analyze(test_counts_norm, method=DBA_ALL_METHODS)

# Show result summary
dba.show(test_counts_norm, bContrasts=T)	

require(rlist)
list.save(test_counts_norm,'test_counts_norm.rds') # fixed error from diffbind_test.DB to test_counts_norm

# test_counts_norm <- readRDS('test_counts_norm.rds')

## Plot results with DiffBind methods
# MA plot
pdf("MAplot.pdf")
dba.plotMA(test_counts_norm)
dev.off()

# volcano plot
pdf("volcanoPlot.pdf")
dba.plotVolcano(test_counts_norm)
dev.off()

# box plot
pdf("boxPlot.pdf")
dba.plotBox(test_counts_norm)
dev.off()

pvals <- dba.plotBox(test_counts_norm)
pvals

# normalized heatmap, no special params
pdf("norm_heatmapPlot.pdf")
dba.plotHeatmap(test_counts_norm) # special params? contrast, correlations, scale, colScheme
dev.off()

# profile plots, need profileplyr to work ==> throws patchworkGrob version error #####
# profile plot, default (separate)
  # library(devtools)
  # devtools::install_github("thomasp85/patchwork") # redundant patchwork packages, make sure to install from github
  # install.packages("profileplyr", dependencies=TRUE) # make sure to setRepositories()

library(profileplyr)

pdf("profilePlot.pdf")
dba.plotProfile(test_counts_norm)
dev.off()

# profile plot, merged 

pdf("profilePlotMerged.pdf")
dba.plotProfile(test_counts_norm, merge=c(DBA_CONTROL))
dev.off()

# DBA significance analysis #####

## Frameshift/control contrast:
res_deseq <- dba.report(test_counts_norm, method=DBA_DESEQ2, contrast = 1, th=1)
head(res_deseq)

library("dplyr")
# Write to file
out <- as.data.frame(res_deseq)
write.table(out, file="fs_cntrl_deseq2_out.txt", sep="\t", quote=F, row.names=F)

# Create bed files for each keeping only significant peaks (p < 0.05)
fs_enrich <- out %>% filter(FDR < 0.05 & Fold > 0) #%>% select(seqnames, start, end)

# Write to file
write.table(fs_enrich, file="fs_enriched.txt", sep="\t", quote=F, row.names=F, col.names=F)

## Splicing/control contrast:
sp_deseq <- dba.report(test_counts_norm, method=DBA_DESEQ2, contrast = 2, th=1) # contrast numbers found with line 76
head(sp_deseq)

# File writing
sp_out <- as.data.frame(spl_deseq)
write.table(sp_out, file="sp_cntrl_deseq2_out.txt", sep="\t", quote=F, row.names=F)

# Filtering for significance into bed files
sp_enrich <- sp_out %>% filter(FDR < 0.05 & Fold > 0) 

#  Final file write
write.table(sp_enrich, file="sp_enriched.txt", sep="\t", quote=F, row.names=F, col.names=F)


#################### Annotation ###################
# install.packages("ChIPseeker")
# install.packages("TxDb.Hsapiens.UCSC.hg38.knownGene")
# install.packages("clusterProfiler")
# install.packages("org.Hs.eg.db")
# install.packages("DOSE")
# 
# # library("TxDb.Hsapiens.UCSC.hg38.knownGene")
# # library("ChIPseeker")
# # library("clusterProfiler")
# # library("org.Hs.eg.db")
# # library("DOSE")
# # library(ChIPpeakAnno)
# 
# gr1 <- toGRanges(fs_enrich, format="BED", header=FALSE)
# 
# ## loading Annotation packages
# annoData <- toGRanges(TxDb.Hsapiens.UCSC.hg38.knownGene, feature="gene")
# 
# annotatedPeak <- annotatePeakInBatch(gr1, AnnotationData = annoData)
# 
# 
# peaks_entrez_id = annotatedPeak$feature
# # save to disk
# write.csv(peaks_entrez_id,file="GREAT_inputs.csv")
# # use GREAT websites to annotate
# 
# 
# library(org.Hs.eg.db)
# over <- getEnrichedGO(peaks_entrez_id, orgAnn="org.Hs.eg.db",
#                       feature_id_type = "entrez_id",
#                       maxP=0.01, minGOterm=10,
#                       multiAdjMethod="BH",
#                       condense=FALSE)
# head(over[["bp"]][, -3])
# head(over[["cc"]][, -3])
# head(over[["mf"]][, -3])




