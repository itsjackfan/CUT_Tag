## CUT_Tag Analysis and Data Processing for setd1a Gene
#### **Description**
Pipeline used to analyze binding sitemaps and frequencies regarding the setd1a gene. The pipeline comes with various functions to format, analyze, and plot data using various different packages (centrally, DiffBind and ChIPQC). 

#### **Usage**
There is 1 central function with various sub-functions that are used to carry out the different steps of the pipeline. Additionally, various plotting functions are also included, and generate various different types of DiffBind graphs in accordance with their name. The functions are as follows:

#### **Functions**
`cut_tag_analysis`  The central, do-all function that makes use of the other 4 sub-functions to read a subsample labeling text file into 7 different plot files, contrast analysis files, and enriched contrast files.      
  **-** `projectDir`: string, The project directory that contains all scripts, documentation, and input files that will be used with the CUT_Tag pipeline.    
  **-** `subSampleFile`: string, The address of a text file containing subsample labeling for all desired samples.    
  **-** `plots`: boolean, Whether or not plots should be included in the final output.       
  **-** `plotTypes`: string or vector of strings, all types of plots that should be included in the final output. Includes `MA`, `volcano`, `box`,`heatmap`,`profplotS`,`profplotM`, or `all` if all plots are desired. A separate function for PCA plots is included due to the difference in datasets.    
      
`environ_prep`  A function that automates the necessary package attachments and directory settings at the beginning of the analysis. Also creates an outputs folder and sets that as the working directory if one does not already exist.     
  **-** `projectDir`: string, The project directory that contains all scripts, documentation, and input files that will be used with the CUT_Tag pipeline.    
    
`info_reading`  A function that generates variables (or reads them from already existing files) based on the subsample labeling text file.
  **-** `subsampleFile`: string, The address of a text file containing subsample labeling for all desired samples.    
      
`analysis`  A function that takes objects generated from `info_reading` and performs DiffBind's significance analysis with options for different contrasts and orderings.    
  **-** `analysisFile`: string, The name of the file or object that the analysis is to be performed on. This is generally the `test_counts_norm` object from the `info_reading` function, though it may be another object that has had `dba.analyze` successfully done on it.    
  **-** `outputFileName`: string, The name of the file the outputs from the function will be saved to.    
  **-** `contrast`: integer, The contrast to be used when analysis is performed. The desired contrast can be performed by running ```{r} dba.show(data, bContrasts=T)``` on the data to be analyzed.    
  **-** `bFlip`: boolean, Whether or not a flip is to be conducted when analyzing the data.     
      
      
      
`dataplots` A function that takes objects generated from `info_reading` and performs the desired DiffBind plotting types, outputting to a PDF.    
  **-** `data`: string, The name of the file or object that is to be plotted. This is generally the `test_counts_norm` object from the `info_reading` function, though it may be another object that has not been normalized with `dba.normalize`.    
  **-** `plots`: string or vector of strings, All types of plots that should be included in the final output. Includes `MA`, `volcano`, `box`,`heatmap`,`profplotS`,`profplotM`, or `all` if all plots are desired. A separate function for PCA plots is included due to the difference in datasets. Defaults to `all`.    
  **-** `pvals`: boolean, Whether or not the pvals of the boxplot are to be saved to a separate object named `pvals`. Defaults to false.    
    
  There are also 7 separate functions that are provided for each of the plot types:    
  `plotMA`  A function used to generate an MA plot based on the given `data` object.       
  `plotVolcano` A function used to generate a volcano plot based on the given `data` object.     
  `plotBox` A function used to generate a box plot based on the given `data` object. Optional inclusion to save the p-values to a variable named `pvals`.      
  `plotHeatmap` A function used to generate a heatmap based on the given `data` object.      
  `plotProfS` A function used to generate a separate profile plot based on the given `data` object.    
  `plotProfM` A function used to generate a merged profile plot based on the given `data` object.     
  `plotPCA` A function used to generate a count-based PCA plot based on the given `data` object. NOTE: This is the only plot function that uses non-normalized, count-based data.     
