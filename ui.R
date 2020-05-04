shinyUI(fluidPage(
  theme = "style.css",
  navbarPage(
    title = span("SS_Shiny", style = "color: steelblue"), 
    ### Data analysis and comparison
    tabPanel("Data analysis", 
             span(h2(strong("Data analysis"), style = "color : steelblue")),
             sidebarLayout(
               sidebarPanel(
                 h3(strong("Data preparation:")),
                 radioButtons(inputId = "exp",
                              label = strong("1. Upload disease-related RNA expression datasets, as the gold standard:"),
                              choices = list("Sample file" = "sf",
                                             "Tarca's datasets"="td",
                                             "Browse"="browse"),
                              selected = c("sf")),
                 fileInput(inputId = "userdata",
                           label = "Browse users' datasets: ",
                           accept = ".rds"),
                 # tags$hr(),
                 radioButtons(inputId = "targetpath",
                              label = strong("2. Upload disease/target pathways:"),
                              choices = list("Sample file" = "sf1",
                                          "Tarca's datasets"="td1",
                                          "Browse"="browse1"),
                              selected = "sf1"),
                 fileInput(inputId = "userdata1",
                           label = "Browse users' datasets: ",
                           accept = c(".rds",
                                      ".gmt")),
                 # tags$hr(),
                 h3(strong("Method selection")),
                 checkboxGroupInput(inputId = "GSEAMethods",
                                    label = strong("3. GSA Methods(at least choose 3 methods):"),
                                    selected = c("plage",
                                                 "zscore",
                                                 "ssgsea",
                                                 "gsva",
                                                 "GRAPE"),
                                    inline = T, width = NULL, 
                                    choiceNames = list(
                                      tags$a(href = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1261155/", "PLAGE"),
                                      tags$a(href = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2563693/", "ZSCORE"),
                                      tags$a(href = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2783335/", "SSGSEA"),
                                      tags$a(href = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3618321/", "GSVA"),
                                      tags$a(href = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5485588/", "GRAPE")),
                                    choiceValues = c("plage",
                                                     "zscore",
                                                     "ssgsea",
                                                     "gsva",
                                                     "GRAPE")),
                 # tags$hr(),
                 radioButtons("comp", label = strong("4. Method to combine P-Values"),
                              choices = list("sumlog((Fisher’s method)"="sumlog", 
                                             "sumz(Stouffer’s method)(Run with at least 4 methods)" = "sumz",
                                             "meanp(Run with at least 4 methods)" = "meanp"), 
                              selected = "sumlog"),
                 # tags$hr(),
                 strong("5. Benchmark metrics:"),
                 checkboxInput(inputId = "sn","Sensitivity comparison", value = T),
                 checkboxInput(inputId = "sp","Specificity comparison"),
                 checkboxInput(inputId = "pr","Precision comparison"),
                 checkboxInput(inputId = "roc","ROC Curve"),
                 actionButton(inputId = "submit", label="Submit", icon("fas fa-magic"))
                 ),
									mainPanel(
									  tabsetPanel(type = "tab",
									              tabPanel(h3("Preview of datasets"), 
									                       h4(strong("1. Expression Data Preview(Only one dataset shown here)")),
									                       shiny::textOutput("text1"),
									                       DT::dataTableOutput("expTable"),
									                       h4(strong("2. Target pathway Preview(Only one list of target pathway here)")),
									                       verbatimTextOutput("Target")), # Show file contents.
									              tabPanel(h3("Preview of Results"),
									                       h4(strong("Results of Sensitivity, specificity and precision:")),
									                       shiny::textOutput("text2"),
									                       verbatimTextOutput("result"), # Show file contents.
									                       h4(strong("Comparison plot")),
									                       shiny::textOutput("text3"),
									                       shiny::plotOutput("cplot", width = "100%", height = "500px")),
									              tabPanel(h3("Download"),          
									                       downloadButton('downloadData', 'Download result RData'))
									              )
									  )
									)),
    ### Benchmark Study Results
    tabPanel("Example: Benchmark study",
             span(h2(strong("Results of SS-Benchmark Study"), style = "color : steelblue")),
             sidebarLayout(
               sidebarPanel(
                  includeMarkdown("WWW/mdfiles/Intro.md")
                 ),
               
               mainPanel(
                 tabsetPanel(type = "tab",
                             tabPanel(h3("Datasets Overview"),
                                      h4(strong("The following are the transcriptomic datasets chosen for this study:")),
                                      tags$img(src="image/data.png", align="center")),
                             tabPanel(h3("SS_Benchmark inputs"),
                                      h4(strong("A screenshot of all the options selected to perform these plots：")),
                                      tags$img(src="image/screen.png", align="center")),
                             tabPanel(h3("Comparison Plots"),
                                      h4(strong("Comparison plots:")),
                                      tags$img(src="image/allplots.png", align="center",
                                               width = "125%", height = "1250px"),
                                      p("Fig : A comparison of sensitivity, specificity, precision and ROC plot of 8 gene set analysis methods."))
                                     )))),
	  ### Contact
	  tabPanel("Help", 
	         span(h2(strong("SS_Shiny Help Center"), style = "color : steelblue")),
	         includeMarkdown("WWW/mdfiles/Help.md"))
    )))
