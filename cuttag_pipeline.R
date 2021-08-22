################### CUT&TAG DATA ANALYSIS MODULE FUNCTIONS #################

cut_tag_analysis <- function(projectDir, subSampleFile, plots=FALSE, plotTypes=NULL) {
  environ_prep(projectDir)
  info_reading(subsampleFile)
  
  tc_norm <- readRDS('test_counts_norm.rds')
  
  analysis(analysisFile = tc_norm, outputFileName = "cntrl_fs_deseq2_out", contrast = 1, bFlip = TRUE)
  analysis(analysisFile = tc_norm, outputFileName = "cntrl_sp_deseq2_out", contrast = 2, bFlip = TRUE)
  
  if (plots == TRUE) {
    dataplots(tc_norm,plotTypes,DBA_CONTROL)
  }
}

environ_prep <- function(projectDir) {
  library(BiocParallel)
  register(DoparParam())
  registered()
  bpparam("SerialParam")
  
  library(rlist)
  library(DiffBind)
  library(ChIPQC)
  library("TxDb.Hsapiens.UCSC.hg38.knownGene")
  library(dplyr)
  
  ifelse(!dir.exists(file.path(projectDir,"outputs")), dir.create(file.path(projectDir,"outputs")), FALSE)
}

info_reading <- function(subsampleFile) {
  if (file.exists('test.rds')) {
    test_counts <- readRDS('test.RDS')
  } else {
    test <- dba(sampleSheet=subsampleFile)
    
    test_qc = ChIPQC(test,annotation="hg38")
    
    test_counts <- dba.count(test)
    list.save(test_counts, 'test.rds')
  }
  
  if (file.exists('test_counts_norm.rds')) {
    test_counts_norm <- readRDS('test_counts_norm.rds')
  } else {
    test_counts_norm <- dba.normalize(test_counts)
    norm <- dba.normalize(test_counts_norm, bRetrieve=TRUE)
    test_counts_norm <- dba.contrast(test_counts_norm, categories=DBA_CONDITION,minMembers=2)
    
    test_counts_norm <- dba.analyze(test_counts_norm, method=DBA_ALL_METHODS)
    list.save(test_counts_norm, 'test_counts_norm.rds')
  }
}

analysis <- function(analysisFile, outputFileName, contrast, bFlip=FALSE) {
  res_deseq <- dba.report(analysisFile, method=DBA_DESEQ2, contrast = contrast, th = 1, bFlip = bFlip)
  out <- as.data.frame(res_deseq)
  write.table(out, file=paste(outputFileName,".txt",sep=""), sep="\t",quote=F, row.names=F)
  
  enrich <- out %>% filter(FDR < 0.05 & Fold > 0)
  write.table(enrich, file=paste(outputFileName,"_enriched.txt",sep=""),sep="\t",quote=F,row.names=F)
}

dataplots <- function(data, plots = "all", pvals = FALSE, profmerges = NULL) {
  if ("all" %in% plots) {
    plotMA(data)
    plotVolcano(data)
    plotBox(data)
    plotHeatmap(data)
    plotProfS(data)
    plotProfM(data)
  } else if ("MA" %in% plots) {
    plotMA(data)
  } else if ("volcano" %in% plots) {
    plotVolcano(data)
  } else if ("box" %in% plots) {
    plotBox(data, pvals = pvals)
  } else if ("heatmap" %in% plots) {
    plotHeatmap(data)
  } else if ("profplotS" %in% plots) {
    plotProfS(data)
  } else if ("profplotM" %in% plots) {
    plotProfM(data, profmerges)
  }
}

# plot functions #####
plotMA <- function(data) {
  pdf("MAplot.pdf")
  dba.plotMA(data)
  dev.off()
}

plotVolcano <- function(data) {
  pdf("volcanoPlot.pdf")
  dba.plotVolcano(data)
  dev.off()
}

plotBox <- function(data, pvals=FALSE) {
  pdf("boxPlot.pdf")
  dba.plotBox(data)
  dev.off()
  
  if (pvals == TRUE) {
    pvals <- dba.plotBox(test_counts_norm)
    pvals
  }
}

plotHeatmap <- function(data) {
  pdf("norm_heatmapPlot.pdf")
  dba.plotHeatmap(test_counts_norm) # special params? contrast, correlations, scale, colScheme
  dev.off()
}

plotProfS <- function(data) {
  pdf("profilePlot.pdf")
  dba.plotProfile(data)
  dev.off()
}

plotProfM <- function(data, merges) {
  pdf("profilePlotMerged.pdf")
  dba.plotProfile(test_counts_norm, merge = merges)
  dev.off()
}

plotPCA <- function(data) {
  pdf("PCAplot_counted.pdf")
  dba.plotPCA(data,DBA_CONDITION,label=DBA_REPLICATE)
  dev.off()
}
#####

cut_tag_analysis('C:\\Users\\Jack Fan\\Documents\\R\\CUT_Tag','C:\\Users\\Jack Fan\\Documents\\R\\CUT_Tag\\inputs\\subsamples.txt',TRUE,"all")
