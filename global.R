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

### load data in need
pathway = readRDS("data/pathway.RDS")


### sample data
samplefile = readRDS("data/samplefile.RDS")
sample_pvalue = readRDS("data/sample_pvalue .RDS")

### Tarca data
tarca1 = readRDS("data/Tarca1.RDS")
tarca2 = readRDS("data/Tarca2.RDS")
tarca3 = readRDS("data/Tarca3.RDS")
Tarcadata = c(tarca1, tarca2, tarca3)

Tarcadata_pvalue = readRDS("data/Tarcadata_pvalue.RDS")
TarcaGS = readRDS("data/TarcaGS.RDS")

