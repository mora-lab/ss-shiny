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
      cat("It's under building...Please wait!")
    } else if(input$exp == "browse" & !is.null(input$userdata)){
      data = readRDS(input$userdata$datapath)
    }
    return(data)
  })
  
  targetgs <- reactive({
    tardata = c()
    if(input$targetpath == "sf1"){
      tardata = sampleGS
    }
    if(input$targetpath == "td1"){
      cat("It's under building...Please wait!")
    }
    if(input$targetpath == "browse1" & !is.null(input$userdata1)){
      tardata = read_data(input$userdata1$datapath)
      #tardata = head(tardata)
    }
    return(tardata)
  })
  
  GSdata <- reactive({
    gsdata = list("KEGGpathlist" = msigKEGG, 
                  "KEGGGScollection" = KEGGgscollection)
    return(gsdata)
    
  })
  ####################################################
  #############  Output of Data Preview  #############
  ####################################################
  output$expTable <-  DT::renderDataTable(DT::datatable({
    expdata()[[1]]
  }))
  
  output$text1 <-  renderText({
    if(input$exp == "td"){
      paste("It's under building...Please wait!")
    }
  })
  
  output$Target <- renderPrint({
    if(input$targetpath == "sf1"){
      targetgs()
    } else if(input$targetpath == "td1"){
      targetgs()[[1]]
    } else if(input$targetpath == "browse1" & !is.null(input$userdata1)){
      targetgs()[[1]]
    }
    
  })
  
  runner <- observeEvent(input$submit,{
    resOutput()
  })
  ####################################################
  ################  Process Data   ###################
  ####################################################
  resOutput = eventReactive(input$submit, {
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
    # if ((nsampleGS))
    # temp2 = list();
    # temp2 = word(names(nsampleGS), 2, sep="\\.");
    
    for (i in 1:length(nsamplefile)){
      temp1 = list();
      t.res[[i]] = run_methods(nsamplefile[[i]],
                               pathway, 
                               GSEA.Methods = input$GSEAMethods,
                               pvalCombMethod =input$comp);
      p.res[[i]] = t.res[[i]]$summary.pvalue.result;
      score[[i]] = t.res[[i]]$score;
      pval.res[[i]] = t.res[[i]]$pvalue.result;
      names(pval.res)[i] = names(score)[i] = names(p.res)[i] = names(t.res)[i] = names(nsamplefile)[i];
      # temp1 = word(names(nsamplefile)[[i]], 2, sep="\\.")
      # t = nsampleGS[[grep(temp2, names(t))]]
      showNotification(paste("Computing benchmark metrics data!"), duration = 20, type = "message");
      SSP.res[[i]] = SSP_calculation(p.res[[i]], nsampleGS[[1]]);  
      newSSP.res[[i]] = newSSPresult(SSP.res[[i]]);
      # colnames(newSSP.res[[i]]) = c("sen")
      newSSP.res[[i]]$Datasets = names(p.res)[[i]];
    }
    newSSP.result = do.call(rbind, newSSP.res);
    newSSP.result$Methods = gsub("PVAL.", "", toupper(newSSP.result$Methods))
    rownames(newSSP.result) = NULL
    
    ROCdata = newSSP.result
    ROCdata$FPR = 1-ROCdata$specificity
    ROCdata = ROCdata[,c("Methods", "sensitivity", "FPR")]
    ROCdata = gather(ROCdata, "Type", "Value", -Methods)
    
    result = list("score" = score,
                  "pvalue.result" = pval.res,
                  "SSP.result" = newSSP.result,
                  "ROC.result" = ROCdata)
    
    removeNotification(note)
    showNotification(paste("Computation is done successfully, now you could download the results, but 
                           it will take some time!"), duration = 20, type = "message")
    return(result)
    
  })
  pt1 <- reactive({
    if (!input$sn) return(NULL)
    ggplot(resOutput()$SSP.result, aes(Methods, sensitivity, fill = Methods))+ ## sensitivity
      geom_boxplot()+
      #geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+
      geom_jitter(shape=16, position=position_jitter(0.2))+
      xlab("Methods")+
      ylab("Sensitivity")
  })
  pt2 <- reactive({
    if (!input$sp) return(NULL)
    ggplot(resOutput()$SSP.result, aes(Methods,specificity, fill = Methods))+  ## specificity
      geom_boxplot()+
      #geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+
      geom_jitter(shape=16, position=position_jitter(0.2))+
      xlab("Methods")+
      ylab("Specificity")
  })
  pt3 <- reactive({
    if (!input$pr) return(NULL)
    ggplot(resOutput()$SSP.result, aes(Methods,precision, fill = Methods))+  ## precision
      geom_boxplot() +
      #geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+
      geom_jitter(shape=16, position=position_jitter(0.2))+
      scale_y_sqrt()+
      xlab("Methods")+
      ylab("Precision")
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
    if (is.null(input$plot)) {
      paste("Please choose at least one benchmark metrics") 
    }
    if(!is.null(input$plot)){
      paste("Please wait a moment to get plots!")
    }
  })
  
  output$result = renderPrint({
    if(!is.null(input$sn) | !is.null(input$sp) | !is.null(input$pr) |!is.null(input$roc) ){
      resOutput()$SSP.result
    }
  })
  
  output$cplot <- renderPlot({
    ptlist <- list(pt1(),pt2(),pt3())#,pt4())
    # wtlist <- c(input$sn,input$sp,input$pr, input$roc)
    to_delete <- !sapply(ptlist,is.null)
    
    ptlist <- ptlist[to_delete] 
    # wtlist <- wtlist[to_delete]
    if (length(ptlist)==0) return(NULL)
    ggpubr::ggarrange(plotlist = ptlist)
  })
  
  
  ### Download
  output$downloadData <- downloadHandler(
    filename = function() { paste('result.RDS') },  #Sys.Date(), 
    content = function(con) { saveRDS(resOutput(), con)}
  )
  
})