### packages to be used (using "BiocManager::install()")
installed <- installed.packages()[,"Package"]

packages.required = c("GEOquery", "tidyverse", "GRAPE", "GSVAdata", 
                      "GSEABase","GSVA", "limma", "metap", "pheatmap",
                      "ggplot2", "shinythemes", "shiny", "R.utils", 
                      "plotROC", "ggpubr", "reshape2","stringr")

for (package in packages.required){
  if (!(package %in% installed)) {
    #install.packages(package)
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    
    BiocManager::install(package)
  }
}

library(R.utils);
# library(xlsx); 
library(stringr);
library(GEOquery);
library(tidyverse);
library(reshape2);
library(GRAPE);
library(GSVAdata);
library(GSEABase);
library(GSVA);
library(limma);
library(metap)
library(pheatmap);
library(ggplot2);
library(shinythemes)
library(shiny)
library(ggpubr)
library(plotROC)
library(DT)


### function to be used
source("R/basicfunctions.R")

msigKEGG <- readRDS("data/msigdbKEGG.RDS");
KEGGgscollection <- readRDS("data/KEGGgscollection.RDS");
pathway = list("msigKEGG" = msigKEGG,
               "KEGGgscollection" = KEGGgscollection)
samplefile <- readRDS("data/samplefile.RDS");
sampleGS <- readRDS("data/sampleGS.RDS")
# tarcafile <- readRDS("data/tarcafile.RDS");
# tarcaGS <- readRDS("data/tarcaGS.RDS")

ssBenchdata <- readRDS("data/ssBenchdata.RDS");
sstargetGS <- readRDS("data/sstargetGS.RDS")



