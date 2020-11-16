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
# source("R/get Tarcadata.R")

### load data in need
download.file("http://www.moralab.science/downloads/Chengshu_data/import_data.RData", destfile = "data/import_data.RData")
load("data/import_data.RData")



