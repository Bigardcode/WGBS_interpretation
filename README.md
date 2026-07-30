# Interpretation of WGBS 

WGBS is a powerful tool for providing a comprehensive view of DNA methylation across the genome. Proper data processing, statistical analysis, and biological interpretation are crucial for drawing meaningful conclusions from the results. It can uncover insights into epigenetic regulation, disease mechanisms, and developmental biology.
The bioinformatics analysis and visualization of WGBS include basic data analysis such as methylation level and density analysis, motif analysis and advanced data analysis such as co-expression network analysis.

- A Methylation profiles detection 
Differentially methylated region (DMR) analysis
Differentially hyper-methylated (Hyper-DMR) and hypo-methylated (Hypo-DMR) 
Differentially methylated cytosines (DMCs) analysis
Differentially Methylated Position (DMP) analysis?

- B Biological Interpretation(CpG Annotation)
Transposable Element (TE) Analysis:
Promoter Analysis
Gene Body Analysis

- C Gene Set Enrichment
Biological pathway and gene ontology


---


    - Step_1.......................package installation ....................................................
     
     >>>>>>

     if (!requireNamespace("BiocManager", quietly = TRUE))install.packages("BiocManager")
     BiocManager::install("miceadds")
     save.image() 
     load("")


     - Step_2.........................Loading the_libraries ..........................

      >>>>>> 
     library(BiocManager)
     BiocManager::valid()

     # Main analysis package
     library(methylKit)
     library(genomation)
     library(GenomicRanges)
     library(tidyverse)
     library(AnnotationDbi)
     library(org.Mm.eg.db)
     library(GenomeInfoDb)
     library(GenomicRanges)

     # Annotation package
    library(readr)
    library(usethis)
    library(devtools)
    library(roxygen2)
    library(roxygen2)
    library(ggplot2)
    library(pacman)
    library(tidySummarizedExperiment)
    library(tools)
    library(limma)
    library(edgeR)
