shinyServer(function(input, output, session){
  options(shiny.maxRequestSize=100*1024^2)
  
  source("global.R")
  note = NULL
  observe({
    updateRadioButtons(session, "targetpath",
                       #label = paste("radioButtons label", x),enen
                       choices = input$exp,
                       selected = input$exp)
  })
  
  ####################################################
  ###################   Input data   #################
  ####################################################
  expdata <- reactive({
    data = c()
    if(input$exp == "Sample_file"){
      data = samplefile
    } else if(input$exp == "Tarca_datasets"){
      data = Tarcadata
    } else if(input$exp == "Browse" && !is.null(input$userdata)){
      data = readRDS(input$userdata$datapath)
    }
    return(data)
  })
  
  targetgs <- reactive({
    tardata = c()
    if(input$targetpath == "Sample_file"){
      tardata = TarcaGS$NSCLC
    } else if(input$targetpath == "Tarca_datasets"){
      tardata = TarcaGS
    } else if(input$targetpath == "Browse" && !is.null(input$userdata1)){
      tardata = read_data(input$userdata1$datapath)
    }
    return(tardata)
  })
  
  ####################################################
  #############  Output of Data Preview  #############
  ####################################################
  output$expTable <-  DT::renderDataTable(DT::datatable({
    if(input$exp == "Sample_file"){
      expdata()[[1]]
    } else if(input$exp == "Tarca_datasets"){
      expdata()[[1]]
    } else if(input$exp == "Browse" & !is.null(input$userdata)){
      expdata()[[1]]
    }
  }))
  
  output$Target <- renderPrint({
    if(input$targetpath == "Sample_file" && is.null(input$userdata1)){
      targetgs()
    } else if(input$targetpath == "Tarca_datasets" && is.null(input$userdata1)){
      targetgs()[[1]]
    } else if(input$targetpath == "Browse"&& !is.null(input$userdata1)){
      targetgs()[[1]]
    }
    
  })
  ####################################################
  ################  Process Data   ###################
  ####################################################
  manage.res = function(sample_pvalue){
    nsample_pvalue = c(); pvalue = c()
    nsample_pvalue = comb_pval(sample_pvalue)
    pvalue = nsample_pvalue[,which(colnames(nsample_pvalue) %in% c(input$GSEAMethods, input$comp))]
    return(pvalue)
  }
  runner <- observeEvent(input$submit,{
    resOutput()
  })
  resOutput = eventReactive(input$submit, {
    if (input$exp == "Browse"){
      note = showNotification(paste("Please wait a moment, all datasets are under process!"), 
                              duration = 0, type = "message")
      t.res = list(); p.res = list(); score = list(); SSP.res = list();
      newSSP.res = list(); new_data =c(); new_GS = c()
      # ROCdata = c();
      new_data = expdata();
      new_GS = targetgs();
      for (i in 1:length(new_data)){
        note1 = showNotification(paste("Dataset: ", names(new_data)[i], " is under processing!", sep = "" ), 
                                 duration = 0, type = "message")
        t.res[[i]] = run_methods(new_data[[i]],
                                 pathway, 
                                 GSEA.Methods = input$GSEAMethods);
        p.res[[i]] = t.res[[i]]$pvalue.result;
        score[[i]] = t.res[[i]]$score;
        names(score)[i] = names(p.res)[i] = names(t.res)[i] = names(new_data)[i];

        note2 = showNotification(paste("Computing benchmark metrics data!"), duration = 0, type = "message");
        temp = c();
        temp = word(names(new_data)[[i]], 2, sep="\\.")
        targetGS = new_GS[names(new_GS) == temp][[1]];
        
        p.res[[i]] = p.res[[i]][,which(colnames(p.res[[i]]) %in% c(input$GSEAMethods, input$comp))];
        
        SSP.res[[i]] = SSP_calculation(p.res[[i]], targetGS);  
        newSSP.res[[i]] = newSSPresult(SSP.res[[i]]);
        newSSP.res[[i]]$Datasets = names(p.res)[[i]];
        showNotification(paste("Dataset: ", names(new_data)[i], " has been processed successfully!", sep = "" ),
                         duration = 0, type = "message")
        removeNotification(note1)
        removeNotification(note2)
      }
      
      newSSP.result = do.call(rbind, newSSP.res);
      rownames(newSSP.result) = NULL
      
      # plotROC
      # ROCdata = newSSP.result
      # ROCdata$FPR = 1-ROCdata$specificity
      # ROCdata = ROCdata[,c("Methods", "sensitivity", "FPR")]
      # ROCdata = gather(ROCdata, "Type", "Value", -Methods)
      
      result = list("score" = score,
                    "pvalue.result" = p.res,
                    "SSP.result" = newSSP.result)#"ROC.result" = ROCdata)
      return(result)
    }
    removeNotification(note)
    showNotification(paste("Computation is done successfully, now you could download the results!"),
                     duration = 0, type = "message")
    
  })
  SSP.result = reactive({
    if (input$exp == "Sample_file") {
      pvalue = c();result = c()
      pvalue = lapply(sample_pvalue, function(x) manage.res(x))
      result = get_SSP(pvalue, sstargetGS$NSCLC)
      rownames(result) = NULL
      return(result)
    } else if (input$exp == "Tarca_datasets") {
      pvalue = c(); SSP.res = list(); newSSP.res = list(); GS = c();
      pvalue = lapply(Tarcadata_pvalue, function(x) manage.res(x))
      GS = TarcaGS;
      for (i in 1:length(pvalue)){
        temp = c();
        temp = word(names(pvalue)[i], 2, sep="\\.")
        targetGS = GS[names(GS) == temp][[1]] ;
        
        SSP.res[[i]] = SSP_calculation(pvalue[[i]], targetGS);  
        newSSP.res[[i]] = newSSPresult(SSP.res[[i]]);
        newSSP.res[[i]]$Datasets = names(pvalue)[[i]];
      }
      newSSP.result = do.call(rbind, newSSP.res);
      rownames(newSSP.result) = NULL
      return(newSSP.result)
    } else if (input$exp == "Browse" && !is.null(input$targetpath)) {
      data = resOutput()$SSP.result
      return(data)
      }
  })
  
  # pt4 <- reactive({
  #     if (!input$roc) return(NULL)
  #   ggplot(resOutput()$ROC.result, aes(d = Type,m = Value, color = Methods))+ 
  #     geom_roc() + 
  #     style_roc(theme = theme_grey, xlab = "Specificity",ylab = "Sensitivity") + ggtitle("ROC") +
  #     ggsci::scale_color_lancet()
  # })
  
  ### Plots of Results
  output$text1 <-  renderText({
    t = paste("You have choose the methods: ")
    print(c(t, input$GSEAMethods,input$comp))
  })
  
  output$text2 <-  renderText({
    if (input$exp == "Sample_file" && input$targetpath == "Sample_file") {
      print(paste("The matrix is ", nrow(SSP.result()), "rows X", ncol(SSP.result()), "columns", sep = " "))
    } else if (input$exp == "Tarca_datasets" && input$targetpath == "Tarca_datasets") {
      print(paste("The matrix is ", nrow(SSP.result()), "rows X", ncol(SSP.result()), "columns! But here only shown the first 10 rows of all results!", sep = " "))
    } else if (input$exp == "Browse" && input$targetpath == "Browse" ){
      print(paste("The matrix is ", nrow(SSP.result()), "rows X", ncol(SSP.result()), "columns! But here only shown the first 10 rows of all results!", sep = " "))
    }
  })
  
  output$text3 <-  renderText({
    if(is.null(input$sn) && is.null(input$sp) && is.null(input$pr)) {
      print("Please choose at least one benchmark metric!")
    }
  })
  
  output$SSP = renderPrint({
    if (length(input$GSEAMethods) == 2){
      cat("Please choose at least 3 methods!")
    } else {
      if (input$exp == "Sample_file" && input$targetpath == "Sample_file") {
        SSP.result()
      } else if (input$exp == "Tarca_datasets" && input$targetpath == "Tarca_datasets") {
        head(SSP.result(),10)
      } else if (input$exp == "Browse" && input$targetpath == "Browse" ){
        head(SSP.result(),10)
      }}
    
  })
  pt1 <- reactive({
    if (!input$sn) return(NULL)
    data = SSP.result()
    senplot = ggplot(data, aes(Methods,Sensitivity, fill = Methods))+ 
      geom_boxplot(outlier.shape = NA,varwidth = TRUE)+
      geom_point(position = position_jitter(0.1)) +
      scale_y_sqrt()+
      xlab("Methods")+
      ylab("Sensitivity")
    senplot
  })
  pt2 <- reactive({
    if (!input$sp) return(NULL)
    data = SSP.result()
    speplot = ggplot(data, aes(Methods,Specificity, fill = Methods))+
      geom_boxplot(outlier.shape = NA,varwidth = TRUE)+
      geom_point(position = position_jitter(0.1)) +
      scale_y_sqrt()+
      xlab("Methods")+
      ylab("Specificity")
    speplot
  })
  pt3 <- reactive({
    if (!input$pr) return(NULL)
    data = SSP.result()
    preplot = ggplot(data, aes(Methods,Precision, fill = Methods))+
      geom_boxplot(outlier.shape = NA,varwidth = TRUE)+
      geom_point(position = position_jitter(0.1)) +
      scale_y_sqrt()+
      xlab("Methods")+
      ylab("Precision")
    preplot
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
    if (input$exp == "Sample_file" && input$targetpath == "Sample_file") {
      data = list("expdata" = expdata(),
                  "targetGS" = targetgs(),
                  "ReferenceGS" = pathway,
                  # "pvalue.result" = pval.res(),
                  "SSP.result" = SSP.result(),
                  "Sensitivity.plot" = pt1(),
                  "Specificity.plot" = pt2(),
                  "Precision.plot" = pt3())
    } else if(input$exp == "Tarca_datasets" && input$targetpath == "Tarca_datasets"){
      data = list("expdata" = expdata(),
                  "targetGS" = targetgs(),
                  "ReferenceGS" = pathway,
                  # "pvalue.result" = pval.res(),
                  "SSP.result" = SSP.result(),
                  "Sensitivity.plot" = pt1(),
                  "Specificity.plot" = pt2(),
                  "Precision.plot" = pt3())
    } else if(input$exp == "Browse" & !is.null(input$userdata)){
      data = list("expdata" = expdata(),
                  "targetGS" = targetgs(),
                  "ReferenceGS" = pathway,
                  # "pvalue.result" = pval.res(),
                  "SSP.result" = SSP.result(),
                  "Sensitivity.plot" = pt1(),
                  "Specificity.plot" = pt2(),
                  "Precision.plot" = pt3())
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