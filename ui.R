library(shiny)
library(shinyjs)
library(visNetwork)
library(shinythemes)
library(shinycssloaders)
library(shinyWidgets)
library(DT)
library(ggplot2)
library(car)
library(nortest)
library(tseries)
library(RcmdrMisc)
library(lmtest)
library(slickR)


datos <-read.csv("www/dataset.csv",dec = ",")
datosEmbriones <- read.csv("www/CPM_embriones_media_filter.txt", sep = "\t", row.names = 1 )
library(shiny)

source("Components/Comp_embryos.R")
source("Components/Comp_Network.R")
shinyUI(fluidPage(title = "BIO114",theme = shinytheme("cyborg"),setBackgroundImage(src = "wallpaper.jpg", shinydashboard = F),tags$head( includeScript("https://d3js.org/d3.v3.min.js") ),
                  tags$style(HTML("td {
  border: 1px solid black;
}")),
                  includeScript(path ="http://d3js.org/d3.v3.min.js"),
                  navbarPage(title = div(),
                             #"Biologia Molecular y Biotecnologia (BIO-114)",style="color:white;text-align:center"
                             tabPanel(div(img(src="LogoGrupo.svg",width="200px",height="86px"))),
                             
                             tabPanel("Embryos",
                                      tabsetPanel(
                                      tabPanel("Expression View",br(),sidebarLayout(
                                        
                                        sidebarPanel(
                                          
                                          textInput("seqID",p("ID gene database:",style="color:white; text-align:center")),
                                          imageOutput("legendtHeatE",height = "150px"),
                                          br(),
                                          textInput("keggID",p("Kegg ID:",style="color:white; text-align:center")),
                                          tableOutput("tableKeggSeq"),
                                          
                                          width = 3
                                        ),
                                        mainPanel(
                                          fluidRow(
                                            uiOutput("descriptionSeq"),br()
                                          ),
                                          fluidRow(
                                            br(),htmlOutput("heatMapEmbryos"),br()
                                          ),
                                          fluidRow(
                                            useShinyjs(),
                                            
                                            br(),column(plotOutput("barplotEmbryos"),width = 6),br(),
                                            br(),column(fluidRow(tableOutput("tableDiffEmbryos"),br(),tableOutput("tableGO")),width = 6),br()
                                            
                                          )
                                          
                                        ))),
                                      tabPanel("Correlation",br(),sidebarLayout(
                                        sidebarPanel(
                                          
                                          numericInput("CorteCorE",value = 0,max = 1, min = -1, step = 0.1,p("Cutoff correlation value:",style="color:white; text-align:center")),
                                          numericInput("CortePvalueE",value = 0.05,min = 0, max = 0.06,step = 0.01,p("Cutoff p-value value:",style="color:white; text-align:center")),
                                          br(),
                                          actionButton("buttonCorE", "Generate correlation table"),
                                          width = 3
                                        ),mainPanel(
                                          DT::dataTableOutput("tableCorE") %>% withSpinner(color="#0dc5c1")
                                        ))),
                                      tabPanel("General",br(),sidebarLayout(
                                        sidebarPanel(
                                          textInput("goETotal",p("Gene Onthology ID:",style="color:white; text-align:center")),
                                          selectInput("moduleETotal", h5("Module selection"), 
                                                      choices = moduleList(), selected = 1),
                                          checkboxGroupInput("checkGroupE", 
                                                             h5("Select diferential condition"), 
                                                             choices = list("PC-ES1 Zygotic embryos" = "PC-ES1,Zygotic embryos", 
                                                                            "EC-ES2 Zygotic embryos" = "EC-ES2,Zygotic embryos", 
                                                                            "C-ES3 Zygotic embryos" = "C-ES3,Zygotic embryos",
                                                                            "PC-ES1 Somatic embryos" = "PC-ES1,Somatic embryos",
                                                                            "EC-ES2 Somatic embryos" = "EC-ES2,Somatic embryos",
                                                                            "C-ES3 Somatic embryos" = "C-ES3,Somatic embryos")),
                                          br(),
                                          actionButton("buttonGeneralE", "Generate table"),
                                          br(),
                                          hr(),
                                          numericInput("clusterHeatE","HCluster cutree",value = 5),
                                          br(),
                                          fileInput("fileHeatmapEmbryos", "Choose csv File", accept = ".csv"),
                                          downloadButton("downloadHeatEmbryos", "Download heatmap"),
                                          width = 3
                                        ),mainPanel(
                                          DT::dataTableOutput("tableGeneralE") %>% withSpinner(color="#0dc5c1")
                                        ))),
                                      tabPanel("GraphTool",
                                               sidebarPanel(br(),
                                                              fileInput("fileNetEmbryos", "Choose csv File", accept = ".csv"),
                                                              textInput("seqID_E_N",p("ID gene database:",style="color:white; text-align:center")),
                                                              textInput("geneName_E_N",p("Node name:",style="color:white; text-align:center")),
                                                              actionButton("buttonAddNetwork", "Add Gene"),
                                                              br(),br(),
                                                              actionButton("buttonDeleteNetwork", "Delete Gene"),
                                                              br(),
                                                              hr(),
                                                              br(),
                                                              numericInput("CorteCorNetE",value = 0,max = 1, min = -1, step = 0.1,p("Cutoff correlation value:",style="color:white; text-align:center")),
                                                              numericInput("CortePvalueNetE",value = 0.05,min = 0, max = 0.06,step = 0.01,p("Cutoff p-value value:",style="color:white; text-align:center")),
                                                              actionButton("buttonFilterNetwork", "Filter"),
                                                              br(),
                                                              hr(),
                                                              br(),
                                                              actionButton("buttonEmptyNetwork", "Delete Network")
                                                            ),
                                               mainPanel(
                                                 visNetworkOutput("network"),
                                                 style = "height:700px;background-color:#FCFBF5;"
                                               ),
                                               br(),br(),
                                               )
                                      
                                      ),br(),br(),
                                      
                             ),
                             tabPanel("Needles",
                                      
                                      tabsetPanel(
                                        tabPanel("Expression View",br(),sidebarLayout(
                                          
                                          sidebarPanel(
                                            
                                            textInput("seqID_N",p("ID gene database:",style="color:white; text-align:center")),
                                            imageOutput("legendtHeatN",height = "150px"),
                                            br(),
                                            textInput("keggIDN",p("Kegg ID:",style="color:white; text-align:center")),
                                            tableOutput("tableKeggSeqN"),
                                            
                                            width = 3
                                          ),
                                          mainPanel(
                                            fluidRow(
                                              uiOutput("descriptionSeqN"),br()
                                            ),
                                            fluidRow(
                                              column(h5(p("May needles",style="color:white;text-align:center"),style="background-color:black;border-radius: 8px"),
                                                     imageOutput("heatNeedles_M"),width = 6),
                                              column(h5(p("November needles",style="color:white;text-align:center"),style="background-color:black;border-radius: 8px"),
                                                     imageOutput("heatNeedles_N"),width = 6)
                                              
                                            ),
                                            fluidRow(
                                              useShinyjs(),
                                              
                                              br(),column(plotOutput("barplotNeedles"),width = 6),br(),
                                              br(),column(fluidRow(tableOutput("tableDiffNeedles"),br(),tableOutput("tableGON")),width = 6),br()
                                              
                                            )
                                            
                                          ))),
                                        tabPanel("Correlation",br(),sidebarLayout(
                                          sidebarPanel(
                                            br(),
                                            actionButton("buttonCorN", "Generate correlation table"),
                                            width = 3
                                          ),mainPanel(
                                            DT::dataTableOutput("tableCorN") %>% withSpinner(color="#0dc5c1")
                                          ))),
                                        tabPanel("General",br(),sidebarLayout(
                                          sidebarPanel(
                                            textInput("goNTotal",p("Gene Onthology ID:",style="color:white; text-align:center")),
                                            selectInput("moduleNTotal", h5("Module selection"), 
                                                        choices = moduleList(), selected = 1),
                                            checkboxGroupInput("checkGroupN", 
                                                               h5("Select diferential condition"), 
                                                               choices = list("May > November" = "M>N", 
                                                                              "November > May" = "N>M", 
                                                                              "Whorl 0 > Whorl 1" = "W0>W1",
                                                                              "Whorl 1 > Whorl 0" = "W1>W0",
                                                                              "Whorl 3 > Whorl 1" = "W3>W1",
                                                                              "Whorl 1 > Whorl 3" = "W1>W3",
                                                                              "Whorl 3 > Whorl 0" = "W3>W0",
                                                                              "Whorl 0 > Whorl 3" = "W0>W3")),
                                            br(),
                                            actionButton("buttonGeneralN", "Generate table"),
                                            br(),
                                            hr(),
                                            numericInput("clusterHeatN","HCluster cutree",value = 5),
                                            br(),
                                            fileInput("fileHeatmapNeedles", "Choose csv File", accept = ".csv"),
                                            downloadButton("downloadHeatNeedles", "Download heatmap"),
                                            width = 3
                                          ),mainPanel(
                                            DT::dataTableOutput("tableGeneralN") %>% withSpinner(color="#0dc5c1")
                                          ))),
                                        tabPanel("GraphTool",
                                                 sidebarPanel(br(),
                                                              fileInput("fileNetNeedles", "Choose csv File", accept = ".csv"),
                                                              textInput("seqID_N_N",p("ID gene database:",style="color:white; text-align:center")),
                                                              textInput("geneName_N_N",p("Node name:",style="color:white; text-align:center")),
                                                              actionButton("buttonAddNetworkN", "Add Gene"),
                                                              br(),br(),
                                                              actionButton("buttonDeleteNetworkN", "Delete Gene"),
                                                              br(),
                                                              hr(),
                                                              br(),
                                                              numericInput("CorteCorNetN",value = 0,max = 1, min = -1, step = 0.1,p("Cutoff correlation value:",style="color:white; text-align:center")),
                                                              numericInput("CortePvalueNetN",value = 0.05,min = 0, max = 0.06,step = 0.01,p("Cutoff p-value value:",style="color:white; text-align:center")),
                                                              actionButton("buttonFilterNetworkN", "Filter"),
                                                              br(),
                                                              hr(),
                                                              br(),
                                                              actionButton("buttonEmptyNetworkN", "Delete Network")
                                                 ),
                                                 mainPanel(
                                                   visNetworkOutput("networkN"),
                                                   style = "height:700px;background-color:#FCFBF5;"
                                                 ),
                                                 br(),br(),
                                        )
                                        
                                      ),br(),br(),
                             )
                             
                  ), class='flex-center'
                  
                  
))