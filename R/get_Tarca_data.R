get_Tarcadata = function(expdata, GPL){  
  
  ### Download and filter microarray data
  exprSet = as.data.frame(exprs(expdata));
  pdata = as.data.frame(pData(expdata));
  pdata$Group = gsub("c", "1", pdata$Group)
  pdata$Group = gsub("d", "0", pdata$Group)
  
  ### Read or Download GPL file, and load it
  gpl = getGEO(GPL,destdir=".")
  gpl_data = Table(gpl);
  symbol = gpl_data$`Gene Symbol`;
  gpl_data$symbol = unlist(lapply(symbol, function(x) word(x, 1, sep="\\ /// ")));  # extract symbols
  gpl_data = gpl_data[, which(colnames(gpl_data) %in% c("ID","symbol"))];
  
  new_exprSet = exprSet
  new_exprSet$ID = as.character(rownames(new_exprSet));
  
  #### convert probe ids into 
  newdata = new_exprSet %>% 
    merge(.,gpl_data, by="ID") %>% #merge informatio of probe_id
    dplyr::select(-ID) %>% #remove extra column    
    dplyr::select(symbol, everything()) %>% #re-arrange
    dplyr::mutate(rowMean = rowMeans(.[grep("GSM", names(.))])) %>% #calculate means
    dplyr::arrange(desc(rowMean))  %>% #order the means
    distinct(symbol,.keep_all = T) %>% # keep the first symbol information
    dplyr::select(-rowMean) #  #remove the new colmun
  newdata = as.data.frame(newdata);
  newdata = newdata[which(newdata$symbol != ""),]
  rownames(newdata) = newdata$symbol ;
  newdata = dplyr::select(newdata, -"symbol") # remove extra column
  newdata = rbind(pdata$Group,newdata);
  rownames(newdata)[1] = "Normal" ; 
  newdata = na.omit(newdata);
  newdata = data.matrix(newdata);
  newdata
  }

mysets = data(package="KEGGdzPathwaysGEO")$results[,"Item"]
data(list=mysets)
mysets1 = data(package="KEGGandMetacoreDzPathwaysGEO")$results[,"Item"]
data(list=mysets1)

tarcadata96 = list("GSE9476" = GSE9476,
                   "GSE3585" = GSE3585,
                   "GSE6956AA" = GSE6956AA,
                   "GSE6956C" = GSE6956C,
                   "GSE781" =  GSE781,
                   "GSE1297" = GSE1297,
                   "GSE20291" = GSE20291,
                   "GSE20164" = GSE20164)
tarcadata570 = list("GSE14924_CD4"= GSE14924_CD4,
                    "GSE14924_CD8" = GSE14924_CD8,
                    "GSE24739_G0" = GSE24739_G0,
                    "GSE24739_G1" = GSE24739_G1,
                    "GSE8671" = GSE8671,
                    "GSE9348" = GSE9348,
                    "GSE4183" = GSE4183,
                    "GSE23878" = GSE23878,
                    "GSE4107" = GSE4107,
                    "GSE1145" = GSE1145,
                    "GSE7305" = GSE7305,
                    "GSE19728" = GSE19728,
                    "GSE21354" = GSE21354,
                    "GSE8762" = GSE8762,
                    "GSE14762" = GSE14762,
                    "GSE5281_EC" = GSE5281_EC,
                    "GSE5281_HIP" = GSE5281_HIP,
                    "GSE5281_VCX" = GSE5281_VCX,
                    "GSE16759" = GSE16759,
                    "GSE18842" = GSE18842,
                    "GSE19188" = GSE19188,
                    "GSE15471" = GSE15471,
                    "GSE16515" = GSE16515,
                    "GSE32676" = GSE32676,
                    "GSE20153"= GSE20153,
                    "GSE3467" = GSE3467,
                    "GSE3678" = GSE3678)				
ntarcadata96 = lapply(tarcadata96, function(x) get_Tarcadata(x, "GPL96"))
ntarcadata570 = lapply(tarcadata570, function(x) get_Tarcadata(x, "GPL570"))

data = c(ntarcadata96,ntarcadata570)
Tarcadata = list()
Tarcadata = list("GSE14924_CD4.AML" = data$GSE14924_CD4,
                 "GSE14924_CD8.AML" = data$GSE14924_CD8,
                 "GSE9476.AML" = data$GSE9476,
                 "GSE5281_EC.Alzheimer" = data$GSE5281_EC,
                 "GSE5281_HIP.Alzheimer" = data$GSE5281_HIP,
                 "GSE5281_VCX.Alzheimer" = data$GSE5281_VCX,
                 "GSE16759.Alzheimer" = data$GSE16759,
                 "GSE1297.Alzheimer" = data$GSE1297,
                 "GSE24739_G0.CML" = data$GSE24739_G0,
                 "GSE24739_G1.CML" = data$GSE24739_G1,
                 "GSE8671.Colorectal" = data$GSE8671,
                 "GSE9348.Colorectal" = data$GSE9348,
                 "GSE4183.Colorectal" = data$GSE4183,
                 "GSE23878.Colorectal" = data$GSE23878,
                 "GSE4107.Colorectal" = data$GSE4107,
                 "GSE1145.DC" = data$GSE1145,
                 "GSE3585.DC" = data$GSE3585,
                 "GSE7305.EC" = data$GSE7305,
                 "GSE19728.Glioma" = data$GSE19728,
                 "GSE21354.Glioma" = data$GSE21354,
                 "GSE8762.Huntington" = data$GSE8762,
                 "GSE18842.NSCLC" = data$GSE18842,
                 "GSE19188.NSCLC" = data$GSE19188,
                 "GSE16515.Pancreatic" = data$GSE16515,
                 "GSE32676.Pancreatic" = data$GSE32676,
                 "GSE20153.Parkinson" = data$GSE20153,
                 "GSE20291.Parkinson" = data$GSE20291,
                 "GSE20164.Parkinson" = data$GSE20164, 
                 "GSE6956AA.Prostate" = data$GSE6956AA,
                 "GSE6956C.Prostate" = data$GSE6956C,
                 "GSE14762.RCC" = data$GSE14762,
                 "GSE781.RCC" = data$GSE781,
                 "GSE3467.Thyroid" = data$GSE3467,
                 "GSE3678.Thyroid" = data$GSE3678)


### reference genesets
### pathway = readRDS("data/pathway.RDS");

### samplefile data, targetGS = TarcaGS[[14]]
### samplefile = readRDS("data/samplefile.RDS");
### sample_pvalue = readRDS("data/sample_pvalue.RDS")

### Tarca's data
### Tarcadata
### TarcaGS = readRDS("data/TarcaGS.RDS")
### Tarcadata_pvalue = readRDS("data/t.RDS")
### names(Tarcadata_pvalue) = names(Tarcadata)
