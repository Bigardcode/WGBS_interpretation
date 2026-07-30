
![me](https://github.com/Bigardcode/WGBS_interpretation/assets/https://www.ivf-spain.com/storage/app/media/uploaded-files/epigenetica-en-1.jpg)

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
--- 
     if (!requireNamespace("BiocManager", quietly = TRUE))install.packages("BiocManager")
     BiocManager::install("miceadds")
     save.image() 
     load("")



- Step_2.........................Loading the_libraries ......................................................
---

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

    
- Step_3...........................set directory ............................................................
---

    setwd("D:/Susia_Student/raheleh/WGBS_part_11_Differential Methylation Analysis/mwethilkt")


- Step_4..............load datasets  ............................................
---   

    ##example 
    P6_1 <- read.delim("P6_1.bismark.cov.gz", header=FALSE, nrows=6)
    P6_1
    ##Define the list containing the Bismark coverage files.
    file.list <- list("P6_1.bismark.cov.gz", "P6_4.bismark.cov.gz", "P8_3.bismark.cov.gz",
    "P8_6.bismark.cov.gz")
    file.list
    # read the listed files into a methylRawList object ........................

    #options; 
    #mm10 = Mus musculus (house mouse)
    #treatment(a vector contatining 0 and 1 denoting which samples are control which samples are test)
    #mincov(minimum read coverage)
    #methylRawList object which stores methylation information per sample for each covered base
    #Coverage= 54+3 = 57
    myobj <- methRead(file.list,
    sample.id=list("Luminal_1","Luminal_2","Basal_1","Basal_2"),
    pipeline = "bismarkCoverage",
    assembly="mm10",
    treatment=c(1,1,0,0),
     mincov = 10)
     myobj


      #What type of data is stored here?
      head(myobj[[2]])
      class(myobj)

      
- Step_5..............Descriptive Statistics ....................................
---

    # Get percentile data by setting plot=FALSE
    getCoverageStats(myobj[[1]], plot=FALSE, both.strands=FALSE)
    #Get a histogram of the methylation percentage per sample
    #Here for sample 1
    #X-axis: Methylation levels (CpG methylation percentage). 
    #Y-axis: Frequency or count of CpG sites at each methylation level. 
    getMethylationStats(myobj[[1]], plot=TRUE, both.strands=FALSE)
    #x-axis representing the coverage depth
    #y-axis indicating the percentage or number of CpG sites within that depth. 
    #Get a histogram of the read coverage per sample
    getCoverageStats(myobj[[1]], plot=TRUE, both.strands=FALSE)

- Step_6..............Filter low coverage below 10...............................
---

    #discards bases that have coverage below 10 reads
    #discards the bases that have more than 99.9th percentile of coverage in each sample.

    myobj.filt <- filterByCoverage(myobj,
    lo.count=10,
    lo.perc=NULL,
    hi.count=NULL,
    hi.perc=99.9)
    getCoverageStats(myobj.filt [[1]], plot=TRUE, both.strands=FALSE)


- Step_7..............Normalization..............................................
---

     #Thus, if the coverage is quite similar across the samples, this step is not really essential, otherwise it might       be a good idea to normalize the data.

     myobj.filt.norm <- normalizeCoverage(myobj.filt, method = "median")
     getCoverageStats(myobj.filt.norm[[1]], plot=TRUE, both.strands=FALSE)


- Step_8..............Merging samples.................................................................
---
    #Merging samples....................
    library(dplyr)
    library(tidyverse)
    library(tidyr)
    library(rlang)
    library(RColorBrewer)
    library(knitr)
    ##we use :: notation to make sure unite() function from methylKit is called
    #unite: unite methylRawList to a single table
    meth =methylKit::unite(myobj, destrand=FALSE)
    meth
    class(meth)

- Step_9.............Differential Methylation. cytosine (DMC)...................................
---
    #test = c("F", "Chisq", "fast.fisher", "midPval")
 
    #adjust = c("SLIM", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr",
    #"none", "qvalue")

    #Table containing pairwise analysis of all significantly differentially methylated cytosines(DMC)
    myDiff=calculateDiffMeth(meth, test = ("fast.fisher"))
    head(myDiff)
    sorted_myDiff <- myDiff[order(myDiff$qvalue),]
    write.csv(soreted_myDiff, file="DMC_soreted_myDiff.csv", row.names=T)
    setwd("D:/Susia_Student/raheleh/WGBS_part_11_Differential Methylation Analysis/mwethilkt/result")
    write.csv(myDiff, file="myDiff.csv", row.names=T)

    #Simple volcano plot to get an overview of differential methylation
    plot(myDiff$meth.diff, -log10(myDiff$qvalue))
    abline(v=0)

    #The x-axis shows the DNA methylation difference 
    #and the y-axis shows the − log¹⁰p-value of each CpG site.

    #Overview of percentage hyper and hypo CpGs per chromosome.
    #Differential methylation by chromosome
    diffMethPerChr(myDiff)
    #after q-value calculation, we can select the differentially methylated 
    #regions/bases based on q-value and percent methylation difference cutoffs. 
    #Following bit selects the bases that have q-value<0.01 and percent methylation 
    #difference larger than 25%. If you specify type="hyper" or type="hypo" options,
    #you will get hyper-methylated or hypo-methylated regions/bases.

    #get hyper methylated bases................................................
    myDiff25p.hyper= getMethylDiff(myDiff,difference=25, qvalue=0.01, type="hyper")
    head(myDiff25p.hyper)
    myDiff25p.hyper <- myDiff25p.hyper[order(myDiff25p.hyper$qvalue),]
   
    #get hypo methylated bases
    myDiff25p.hypo=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hypo")
    head(myDiff25p.hypo)
    myDiff25p.hypo <- myDiff25p.hypo[order(myDiff25p.hypo$qvalue),]


     #get all differentially mentholated bases
     myDiff25p=getMethylDiff(myDiff,difference=25,qvalue=0.01)
     head(myDiff25p)
     myDiff25p <- myDiff25p[order(myDiff25p$qvalue),]
     myDiff25p
     write.csv(myDiff25p, file="DMC_myDiff25p.csv", row.names=T)
     diffMethPerChr(myDiff,plot=FALSE,qvalue.cutoff=0.01, meth.cutoff=25)
     setwd("D:/Susia_Student/raheleh/WGBS_part_11_Differential Methylation Analysis/mwethilkt/result"
     write.csv(myDiff25p.hyper, file="DMC_myDiff25p.hyper.csv", row.names=T)
     write.csv(myDiff25p.hypo, file="DMC_myDiff25p.hypo.csv", row.names=T)
     write.csv(myDiff25p, file="DMC_myDiff25p.csv", row.names=T)
