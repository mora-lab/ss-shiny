### packages to be used (using "BiocManager::install()")
installed = installed.packages()[,"Package"]

packages.required = c("GEOquery", "tidyverse", "GRAPE", "GSVAdata", 
                      "GSEABase","GSVA", "limma", "metap", "pheatmap",
                      "ggplot2", "shinythemes", "shiny", "R.utils", "DT",
                      "plotROC", "ggpubr", "reshape2","stringr", "plotly",
                      "KEGGdzPathwaysGEO", "KEGGandMetacoreDzPathwaysGEO")

for (package in packages.required){
  if (!(package %in% installed)) {
    #install.packages(package)
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    
    BiocManager::install(package, update = FALSE)
  }
}

suppressMessages(library(R.utils)); 
suppressMessages(library(stringr));
suppressMessages(library(GEOquery));
suppressMessages(library(tidyverse));
suppressMessages(library(reshape2));
suppressMessages(library(GRAPE));
suppressMessages(library(GSVAdata));
suppressMessages(library(GSEABase));
suppressMessages(library(GSVA));
suppressMessages(library(limma));
suppressMessages(library(metap));
suppressMessages(library(pheatmap));
suppressMessages(library(ggplot2));
suppressMessages(library(shinythemes));
suppressMessages(library(shiny));
suppressMessages(library(ggpubr));
suppressMessages(library(plotROC));
suppressMessages(library(plotly));
suppressMessages(library(DT));
suppressMessages(library(KEGGdzPathwaysGEO));
suppressMessages(library(KEGGandMetacoreDzPathwaysGEO));

### function to be used
source("R/basicfunctions.R")

### reference genesets
pathway = readRDS("data/pathway.RDS");

### samplefile data, targetGS = TarcaGS[[14]]
samplefile = readRDS("data/samplefile.RDS");
sample_pvalue = readRDS("data/sample_pvalue.RDS")
sample.result = readRDS("data/samplefile_SSP.RDS")

### Tarca's data
nTarcadata1 = readRDS("data/tarca1.RDS");
nTarcadata2 = readRDS("data/tarca2.RDS");
nTarcadata3 = readRDS("data/tarca3.RDS");
Tarcadata = c(nTarcadata1, nTarcadata2, nTarcadata3);
Tarcadata = Tarcadata[-24]

TarcaGS = readRDS("data/nTarcaGS.RDS")
Tarcadata_pvalue = readRDS("data/Tarcadata_pvalue.RDS")# t.RDS
names(Tarcadata_pvalue) = names(Tarcadata)
Tarca.result = readRDS("data/Tarcadata_SSP.RDS")
rownames(Tarca.result) = NULL

### ssbenchmark data
ssBenchdata = readRDS("data/ssBenchdata.RDS");
sstargetGS = readRDS("data/sstargetGS.RDS")



