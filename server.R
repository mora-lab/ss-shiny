shinyServer(function(input, output, session){
  options(shiny.maxRequestSize=100*1024^2)
  source("global.R")
  note = NULL
  ####################################################
  ###################   Input data   #################
  ####################################################
  expdata <- reactive({
    data = c()
    if(input$exp == "sf"){
      data = samplefile
    } else if(input$exp == "td"){
      data = Tarcadata
    } else if(input$exp == "browse" & !is.null(input$userdata)){
      data = readRDS(input$userdata$datapath)
    }
    return(data)
  })
  
  targetgs <- reactive({
    tardata = c()
    if(input$targetpath == "sf1"){
      tardata = TarcaGS[[14]]
    }
    if(input$targetpath == "td1"){
      tardata = TarcaGS
    }
    if(input$targetpath == "browse1" & !is.null(input$userdata1)){
      tardata = read_data(input$userdata1$datapath)
    }
    return(tardata)
  })
  
  SSP.result = reactive({
    if (input$exp == "sf" && input$targetpath == "sf1") {
      data = sample.result
      if(length(input$GSEAMethods) == 5){
        result = data
      } else if(length(input$GSEAMethods) == 4){
        result = data[which(data$Methods %in% input$GSEAMethods),]
      } else if(length(input$GSEAMethods) == 3){
        result = data[which(data$Methods %in% input$GSEAMethods),]
      } 
      return(result)
    } else if (input$exp == "td" && input$targetpath == "td1") {
      data = Tarca.result
      if(length(input$GSEAMethods) == 5){
        result = data
      } else if(length(input$GSEAMethods) == 4){
        result = data[which(data$Methods %in% input$GSEAMethods),]
      } else if(length(input$GSEAMethods) == 3){
        result = data[which(data$Methods %in% input$GSEAMethods),]
      } 
      return(result)
    } else if (input$exp == "browse" && input$targetpath == "browse1" ){
      data = resOutput()$SSP.result
      if(length(input$GSEAMethods) == 5){
        result = data
      } else if(length(input$GSEAMethods) == 4){
        result = data[which(data$Methods %in% input$GSEAMethods),]
      } else if(length(input$GSEAMethods) == 3){
        result = data[which(data$Methods %in% input$GSEAMethods),]
      } 
      return(result)
    }
  })
  
  pval.res <- reactive({
    nsample_pvalue = list()
    pvalue = list()
	nTarcadata_pvalue = list()
	nbrowse_pvalue = list()
    if (input$exp == "sf" && input$targetpath == "sf1") {
      for (i in 1:length(sample_pvalue)){
        pvalue[[i]] = sample_pvalue[[i]][,which(colnames(sample_pvalue[[i]]) %in% input$GSEAMethods)]
        nsample_pvalue[[i]] = comb_pval(pvalue[[i]], method = input$comp)
      }
      nsample_pvalue = as.data.frame(nsample_pvalue)
      colnames(nsample_pvalue) = paste("combined.pval", names(sample_pvalue), sep="")
      return(nsample_pvalue)
    } else if(input$exp == "td" && input$targetpath == "td1") {
      for (i in 1:length(Tarcadata_pvalue)){
        pvalue[[i]] = Tarcadata_pvalue[[i]][,which(colnames(Tarcadata_pvalue[[i]]) %in% input$GSEAMethods)]
        nTarcadata_pvalue[[i]] = comb_pval(pvalue[[i]], method = input$comp)
      }
      nTarcadata_pvalue = as.data.frame(nTarcadata_pvalue)
      colnames(nTarcadata_pvalue) = paste("combined.pval", names(Tarcadata_pvalue),sep="")
      return(nTarcadata_pvalue)
    } else if(input$exp == "browse" && input$targetpath == "browse1" ){
		browse_pvalue = resOutput()$combined.pvalue
      for (i in 1:length(browse_pvalue)){
		pvalue[[i]] = browse_pvalue[[i]][,which(colnames(browse_pvalue[[i]]) %in% input$GSEAMethods)]
        nbrowse_pvalue[[i]] = comb_pval(pvalue[[i]], method = input$comp)
      }
      nbrowse_pvalue = as.data.frame(nbrowse_pvalue)
      colnames(nbrowse_pvalue) = paste("combined.pval", names(nbrowse_pvalue),sep="")
      return(nbrowse_pvalue)
    }    
  })
  ####################################################
  #############  Output of Data Preview  #############
  ####################################################
  output$expTable <-  DT::renderDataTable(DT::datatable({
    if(input$exp == "sf"){
      expdata()[[1]]
    } else if(input$exp == "td"){
      expdata()[[1]][[1]]
    } else if(input$exp == "browse" & !is.null(input$userdata)){
      expdata()[[1]]
    }
  }))

  output$Target <- renderPrint({
    if(input$targetpath == "sf1"){
      targetgs()
    } else if(input$targetpath == "td1"){
      targetgs()[[1]]
    } else if(input$targetpath == "browse1" & !is.null(input$userdata1)){
      targetgs()[[1]]
    }
    
  })
  ####################################################
  ################  Process Data   ###################
  ####################################################
  runner <- observeEvent(input$submit,{
    resOutput()
  })
  resOutput = eventReactive(input$submit, {
    if (input$exp == "browse" && input$targetpath == "browse1" ){
      note = showNotification(paste("Please wait for some minutes, the process is running!"), 
                              duration = 0, type = "message")
      t.res = list();
      p.res = list();
      score = list();
      pval.res = list();
      SSP.res = list();
      newSSP.res = list();
      ROCdata = c();
      nsamplefile = expdata();
      nsampleGS = targetgs();
      for (i in 1:length(nsamplefile)){
        temp1 = list();
        t.res[[i]] = run_methods(nsamplefile[[i]],
                                 pathway, 
                                 GSEA.Methods = input$GSEAMethods,
                                 pvalCombMethod =input$comp);
        p.res[[i]] = t.res[[i]]$summary.pvalue.result;
        score[[i]] = t.res[[i]]$score;
        pval.res[[i]] = t.res[[i]]$pvalue.result;
        combined.pvalue[[i]] = t.res[[i]]$combined.pvalue
        names(pval.res)[i] = names(score)[i] = names(p.res)[i] = names(t.res)[i] = names(nsamplefile)[i];
        # temp1 = word(names(nsamplefile)[[i]], 2, sep="\\.")
        # t = nsampleGS[[grep(temp2, names(t))]]
        showNotification(paste("Computing benchmark metrics data!"), duration = 20, type = "message");
        SSP.res[[i]] = SSP_calculation(p.res[[i]], nsampleGS);  
        newSSP.res[[i]] = newSSPresult(SSP.res[[i]]);
        # colnames(newSSP.res[[i]]) = c("sen")
        newSSP.res[[i]]$Datasets = names(p.res)[[i]];
        }
      newSSP.result = do.call(rbind, newSSP.res);
      newSSP.result$Methods = gsub("PVAL.", "", toupper(newSSP.result$Methods))
      rownames(newSSP.result) = NULL
      
      # plotROC
      # ROCdata = newSSP.result
      # ROCdata$FPR = 1-ROCdata$specificity
      # ROCdata = ROCdata[,c("Methods", "sensitivity", "FPR")]
      # ROCdata = gather(ROCdata, "Type", "Value", -Methods)
      
      result = list("score" = score,
                    "pvalue.result" = pval.res,
                    "combined.pvalue.results" = combined.pvalue,
                    "SSP.result" = newSSP.result)#"ROC.result" = ROCdata)
      
      removeNotification(note)
      showNotification(paste("Computation is done successfully, now you could download the results!"), duration = 20, type = "message")
      return(result)
   }
  })
  pt1 <- reactive({
    if (!input$sn) return(NULL)
    data = SSP.result()
    senplot = ggplot(data, aes(Methods,sensitivity, fill = Methods))+ 
      geom_boxplot(outlier.shape = NA,varwidth = TRUE)+
      geom_point(position = position_jitter(0.1)) +
      scale_y_sqrt()+
      xlab("Sensitivity")+
      ylab("P-values")
    senplot
  })
  pt2 <- reactive({
    if (!input$sp) return(NULL)
    data = SSP.result()
    speplot = ggplot(data, aes(Methods,specificity, fill = Methods))+
      geom_boxplot(outlier.shape = NA,varwidth = TRUE)+
      geom_point(position = position_jitter(0.1)) +
      scale_y_sqrt()+
      xlab("Specificity")+
      ylab("P-values")
    speplot
  })
  pt3 <- reactive({
    if (!input$pr) return(NULL)
    data = SSP.result()
    preplot = ggplot(data, aes(Methods,precision, fill = Methods))+
      geom_boxplot(outlier.shape = NA,varwidth = TRUE)+
      geom_point(position = position_jitter(0.1)) +
      scale_y_sqrt()+
      xlab("Precision")+
      ylab("P-values")
    preplot
    })
  # pt4 <- reactive({
  #     if (!input$roc) return(NULL)
  #   ggplot(resOutput()$ROC.result, aes(d = Type,m = Value, color = Methods))+ 
  #     geom_roc() + 
  #     style_roc(theme = theme_grey, xlab = "Specificity",ylab = "Sensitivity") + ggtitle("ROC") +
  #     ggsci::scale_color_lancet()
  # })
  ### Plots of Results
  output$text2 <-  renderText({
    t = paste("You have choose the methods:")
    print(c(t, input$GSEAMethods))
  })
  
  output$text3 <-  renderText({
    if(is.null(input$sn) && is.null(input$sp) && is.null(input$pr)) {
      t = paste("Please choose at least one benchmark metric:")
      print(t)
    }
  })
  
  output$comPval = renderPrint({
    if (length(input$GSEAMethods) == 2){
      cat("Please choose at least 3 methods!")
    } else {
      if (input$exp == "sf" && input$targetpath == "sf1") {
        head(pval.res(),10)
      } else if (input$exp == "td" && input$targetpath == "td1") {
        pval.res()[1:10,1:3]
      } else if (input$exp == "browse" && input$targetpath == "browse1" ){
        resOutput()$combined.pvalue.results
      }}
    
  })
  
  output$SSP = renderPrint({
    if (length(input$GSEAMethods) == 2){
      cat("Please choose at least 3 methods!")
    } else {
      if (input$exp == "sf" && input$targetpath == "sf1") {
        SSP.result()
      } else if (input$exp == "td" && input$targetpath == "td1") {
        head(SSP.result(),10)
      } else if (input$exp == "browse" && input$targetpath == "browse1" ){
        head(resOutput()$SSP.result,10)
      }}
    
  })
  
  output$cplot <- renderPlot({
    if (length(input$GSEAMethods) == 2){
      cat("Please choose at least 3 methods!")
    } else {
      ptlist <- list(pt1(),pt2(),pt3())#,pt4())
      # wtlist <- c(input$sn,input$sp,input$pr, input$roc)
      to_delete <- !sapply(ptlist,is.null)
      
      ptlist <- ptlist[to_delete] 
      # wtlist <- wtlist[to_delete]
      if (length(ptlist)==0) return(NULL)
      ggpubr::ggarrange(plotlist = ptlist)
    }
  })
  
  
  ### Download
  downdata <- reactive({
    data = c()
    if (input$exp == "sf" && input$targetpath == "sf1") {
      data = list("expdata" = expdata(),
                  "targetGS" = targetgs(),
                  "ReferenceGS" = pathway,
                  "combined.pvalue" = pval.res(),
                  "SSP.result" = SSP.result(),
                  "sensitivity.plot" = pt1(),
                  "specificity.plot" = pt2(),
                  "precision.plot" = pt3())
    } else if(input$exp == "td" && input$targetpath == "td1"){
      data = list("expdata" = expdata(),
                  "targetGS" = targetgs(),
                  "ReferenceGS" = pathway,
                  "combined.pvalue" = pval.res(),
                  "SSP.result" = SSP.result(),
                  "sensitivity.plot" = pt1(),
                  "specificity.plot" = pt2(),
                  "precision.plot" = pt3())
    } else if(input$exp == "browse" & !is.null(input$userdata)){
     data = list("expdata" = expdata(),
                  "targetGS" = targetgs(),
                  "ReferenceGS" = pathway,
                  "combined.pvalue" = resOutput()$combined.pvalue.results,
                  "SSP.result" = resOutput()$SSP.result(),
                  "sensitivity.plot" = pt1(),
                  "specificity.plot" = pt2(),
                  "precision.plot" = pt3())
    }
    return(data)
  })
  
  output$downloadData <- downloadHandler(
    filename = function() { 
      paste('result.RDS') 
    }, 
    
    content = function(con) { 
      saveRDS(downdata(), con)
    }
  )
  
})