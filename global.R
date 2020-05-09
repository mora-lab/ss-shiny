### packages to be used (using "BiocManager::install()")
installed <- installed.packages()[,"Package"]

packages.required = c("GEOquery", "tidyverse", "GRAPE", "GSVAdata", 
                      "GSEABase","GSVA", "limma", "metap", "pheatmap",
                      "ggplot2", "shinythemes", "shiny", "R.utils", 
                      "plotROC", "ggpubr", "reshape2","stringr", "plotly")

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
library(plotly)
library(DT)


### function to be used
source("R/basicfunctions.R")

### reference genesets
pathway <- readRDS("data/pathway.RDS");

### samplefile data, targetGS = TarcaGS[[14]]
samplefile <- readRDS("data/samplefile.RDS");
sample_pvalue <- readRDS("data/sample_pvalue.RDS")
sample.result <- read.csv("results/samplefileSSP.csv", header = T);
rownames(sample.result) = NULL

### Tarca's data
nTarcadata1 <- readRDS("data/nTarcadata1.RDS");
nTarcadata2 <- readRDS("data/nTarcadata2.RDS");
Tarcadata <- c(nTarcadata1, nTarcadata2);
TarcaGS <- readRDS("data/nTarcaGS.RDS")
nTarcadata_pvalue <- readRDS("data/t.RDS")#nTarcadata_pvalue.RDS")
Tarca.result <- read.csv("results/TarcaSSP(259 pathways).csv", header = T);
rownames(Tarca.result) = NULL

### ssbenchmark data
ssBenchdata <- readRDS("data/ssBenchdata.RDS");
sstargetGS <- readRDS("data/sstargetGS.RDS")



