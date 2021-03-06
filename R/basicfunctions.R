read_data = function(file){
  if(grepl("\\.RDS$",file)[1]){ 
    data = readRDS(file)
  } else if(grepl("\\.gmt$",file)[1]){ 
    geneSetDB = readLines(file)                                ##read in the gmt file as a vector of lines
    geneSetDB = strsplit(geneSetDB,"\t")                       ##convert from vector of strings to a list
    names(geneSetDB) = sapply(geneSetDB,"[",1)                 ##move the names column as the names of the list
    geneSetDB = lapply(geneSetDB, "[",-1:-2)                   ##remove name and description columns
    geneSet = lapply(geneSetDB, function(x){x[which(x!="")]})##remove empty strings
    ### Sort the pathway, pathway with more genes comes first.
    geneSet.tmp = lengths(geneSet)
    geneSet.sort = sort(geneSet.tmp, decreasing = FALSE)
    data = list()
    for (i in 1:length(geneSet)){
      data[i] = geneSet[names(geneSet.sort)[i]]
      names(data)[i] = names(geneSet.sort)[i]
    }} 
  return(data)
}

run_GRAPE = function(expData, PathwayList){
    
    ### divide the data into two parts
    ConData = expData[-1,][ , which(expData[1,] == 1)] 
    ConData = as.matrix(ConData) 
    
    psmat = makeGRAPE_psMat(ConData, expData[-1,], PathwayList) 
    colnames(psmat) = colnames(expData) 
    psmat
  
}

run_limma = function(score.result, expData){
  
  design = c();
  contrast.matrix = c();
  fit = c();
  fit2 = c();
  results = c();
  
  ## matrix design
  design = model.matrix(~0+factor(expData[1,], levels=c("0","1")))
  colnames(design) = c("D0", "N1")
  rownames(design) = colnames(expData) 
  contrast.matrix = makeContrasts(N1-D0, levels=design)
  fit = lmFit(score.result, design)
  fit2 = contrasts.fit(fit, contrast.matrix)
  fit2 = eBayes(fit2)
  results = topTable(fit2, adjust="BH", n = Inf)
  results = as.matrix(results[order(rownames(results)),4])
  rownames(results) = rownames(score.result[order(rownames(score.result)),])
  colnames(results) = "P.value"
  results
  #write.table(results,paste0("limma_DEG.",contrast,".txt"),quote = F,sep = '\t',row.names = T)
}
comb_pval = function(P.res){
  
  com_pvals = c();
  sumlog_temp = matrix(0,nrow(P.res),1)
  sumz_temp = matrix(0,nrow(P.res),1)
  meanp_temp = matrix(0,nrow(P.res),1)
  ntemp = data.frame()
  
  for (i in 1:nrow(P.res)){
    sumlog_temp[i,1] = sumlog(as.numeric(P.res[i,]))$p
    sumz_temp[i,1] = sumz(as.numeric(P.res[i,]))$p
    meanp_temp[i,1] = meanp(as.numeric(P.res[i,]))$p
  }
  
  sumlog_temp = as.data.frame(sumlog_temp)
  sumz_temp = as.data.frame(sumz_temp)
  meanp_temp = as.data.frame(meanp_temp)
  colnames(sumlog_temp) = "sumlog"
  colnames(sumz_temp) = "sumz"
  colnames(meanp_temp) = "meanp"
  
  ntemp = cbind(P.res, sumlog_temp, sumz_temp, meanp_temp)
  rownames(ntemp) = rownames(P.res)
  return(ntemp)
}

run_methods = function(expData, pathway, GSEA.Methods){
  
  if (length(GSEA.Methods) == 1){stop("Please choose two GSEA methods!")}
  
  result = list()
  pval.res = list()
  output = list();
  
  data = expData[-1,];
  
  if ("PLAGE" %in% GSEA.Methods){
    showNotification(paste("Running PLAGE!"), duration = 20, type = "message")
    result$res.plage = GSVA::gsva(as.matrix(data), pathway$KEGGgscollection, method = "plage", 
                            mx.diff = FALSE, parallel.sz=2, abs.ranking = FALSE, verbose=FALSE)
    showNotification(paste("PLAGE is done successfully!"), duration = 20, type = "message")
    pval.res$pval.plage = run_limma(result$res.plage, expData)
  } 
  if ("ZSCORE" %in% GSEA.Methods){
    showNotification(paste("Running ZSCORE!"), duration = 20, type = "message")
    result$res.zscore = GSVA::gsva(as.matrix(data), pathway$KEGGgscollection, method = "zscore", 
                             mx.diff = FALSE, parallel.sz=2, abs.ranking = FALSE, verbose=FALSE)
    showNotification(paste("ZSCORE is done successfully!"), duration = 20, type = "message")
    pval.res$pval.zscore = run_limma(result$res.zscore, expData)
  } 
  if ("SSGSEA" %in% GSEA.Methods){
    showNotification(paste("Running SSGSEA!"), duration = 20, type = "message")
    result$res.ssgsea = GSVA::gsva(as.matrix(data), pathway$KEGGgscollection, method = "ssgsea", 
                             mx.diff = FALSE, parallel.sz=2, abs.ranking = FALSE, verbose=FALSE)
    showNotification(paste("SSGSEA is done successfully!"), duration = 20, type = "message")
    pval.res$pval.ssgsea = run_limma(result$res.ssgsea, expData)
  } 
  if ("GSVA" %in% GSEA.Methods){
    showNotification(paste("Running GSVA!"), duration = 20, type = "message")
    result$res.gsva = GSVA::gsva(as.matrix(data), pathway$KEGGgscollection, method = "gsva", 
                           mx.diff = FALSE, parallel.sz=2, abs.ranking = FALSE, verbose=FALSE)
    showNotification(paste("GSVA is done successfully!"), duration = 20, type = "message")
    pval.res$pval.gsva = run_limma(result$res.gsva, expData)
  } 
  if ("GRAPE" %in% GSEA.Methods){
    showNotification(paste("Running GRAPE!"), duration = 20, type = "message")
    result$res.GRAPE = run_GRAPE(expData, pathway$msigKEGG)
    showNotification(paste("GRAPE is done successfully!"), duration = 20, type = "message")
    pval.res$pval.GRAPE = run_limma(result$res.GRAPE, expData)
  }
  
  new.res.pval = c()
  new.res.pval = as.data.frame(pval.res)
  colnames(new.res.pval) = toupper(names(pval.res))
  colnames(new.res.pval) = gsub("PVAL.", "", colnames(new.res.pval))
  showNotification(paste("Computing combined pvalues!"), duration = 20, type = "message")
  combined.pvalue = comb_pval(new.res.pval)
  showNotification(paste("Combined pvalues has computed successfully!"), duration = 20, type = "message")
  
  output$score = result
  # output$pvalue.result = pval.res
  # output$summary.pvalue.result = new.res.pval
  output$pvalue.result = combined.pvalue
  return(output)
}

SSP_calculation = function(Pval.res, target.pathway){
  
  cond_below_0.05 = c() 
  cond_over_0.05 = c() 
  
  true_positives = list() 
  true_positives_data = c() 
  false_negatives = c() 
  false_negatives_data = c() 
  true_negatives_ids = list() 
  false_positives1_ids = list() 
  false_positives2_ids = list() 
  false_positives = list() 
  
  greater_than_0.05 = c() 
  less_than_0.05 = c() 
  true_negatives_ids = c() 
  false_positives = c() 
  sensitivity = c() 
  specificity = c() 
  precision = c() 
  
  sensitivity_result = c() 
  specificity_result = c() 
  precision_result = c() 
  
  Pval.res = as.matrix(Pval.res)
  for (i in 1:ncol(Pval.res)) {
    
    cond_below_0.05[[i]] = rownames(Pval.res[which(Pval.res[,i] <= 0.05),]) %in% toupper(names(target.pathway))
    cond_over_0.05[[i]] = rownames(Pval.res[which(Pval.res[,i] > 0.05),]) %in% toupper(names(target.pathway))
    
    true_positives_data  = Pval.res[which(cond_below_0.05[[i]]),]
    if(class(true_positives_data) == "matrix"){ true_positives = as.numeric(length(true_positives_data[,i])) }
    if(class(true_positives_data) == "numeric"){ true_positives = 1 }
    false_negatives_data  = Pval.res[which(cond_over_0.05[[i]]),]
    if(class(false_negatives_data) == "matrix"){ false_negatives = as.numeric(length(false_negatives_data[,i])) }
    if(class(false_negatives_data) == "numeric"){ false_negatives = 1 }
    
    ## Tool results' subsets on the basis of statistical significance.
    greater_than_0.05[[i]] = as.data.frame(Pval.res[which(Pval.res[,i] > 0.05), i]) 
    less_than_0.05[[i]] = as.data.frame(Pval.res[which(Pval.res[,i] <= 0.05),i]) 
    true_negatives_ids[[i]] = setdiff(rownames(greater_than_0.05[[i]]), names(target.pathway))  ## All ids that are there in the tool result with p > 0.05 and absent in the disease pool.
    
    false_positives2_ids = setdiff(rownames(less_than_0.05[[i]]), names(target.pathway)) 
    false_positives[[i]] = false_positives2_ids 
    
    sensitivity[[i]] = true_positives/(true_positives + false_negatives) 
    specificity[[i]] = length(true_negatives_ids[[i]])/(length(true_negatives_ids[[i]]) + length(false_positives[[i]])) 
    precision[[i]] = true_positives/(true_positives + length(false_positives[[i]]))
  }
  
  sensitivity.re = data.frame(t(sapply(sensitivity,c)))
  specificity.re = data.frame(t(sapply(specificity,c)))
  precision.re = data.frame(t(sapply(precision,c)))
  colnames(sensitivity.re) = colnames(specificity.re) = colnames(precision.re) = colnames(Pval.res)
  # rownames(sensitivity.re) = rownames(specificity.re) = rownames(precision.re) = "value"
  
  result = list("sensitivity.result" = sensitivity.re, 
                "specificity.result" = specificity.re,
                "precision.result" = precision.re)
  return(result)
}

newSSPresult = function(SSP_cal) {
  SSP = c();
  SSP.sen = t(SSP_cal$sensitivity.result)
  SSP.spe = t(SSP_cal$specificity.result)
  SSP.pre = t(SSP_cal$precision.result)
  SSP$Methods = rownames(SSP.sen)
  SSP$Sensitivity = SSP.sen
  SSP$Specificity = SSP.spe
  SSP$Precision = SSP.pre
  return(as.data.frame(SSP))
  
}

get_SSP = function(Pval.res, target.pathway){
  
  SSP.res = list()
  newSSP.res = list()
  newSSP.result = c()
  
  for( i in 1:length(Pval.res)){
    SSP.res[[i]] = SSP_calculation(Pval.res[[i]], target.pathway);  
    newSSP.res[[i]] = newSSPresult(SSP.res[[i]]);
    newSSP.res[[i]]$Datasets = names(Pval.res)[[i]];
    }
  
  newSSP.result = do.call(rbind, newSSP.res);
  # newSSP.result$Methods = gsub("PVAL.", "", toupper(newSSP.result$Methods))
  rownames(newSSP.result) = NULL
  return(newSSP.result)
}
