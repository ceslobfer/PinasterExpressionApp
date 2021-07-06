library(shiny)
library(shinyjs)
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
shinyUI(fluidPage(title = "BIO114",theme = shinytheme("cyborg"),setBackgroundImage(src = "wallpaper.jpg", shinydashboard = F),
                  tags$style(HTML("td {
  border: 1px solid black;
}")),
                  includeScript(path ="http://d3js.org/d3.v3.min.js"),
                  navbarPage(title = div(),
                             #"Biologia Molecular y Biotecnologia (BIO-114)",style="color:white;text-align:center"
                             tabPanel(div(img(src="LogoGrupo.svg",width="200px",height="86px")),
                                      fluidRow(column(div(),width = 2),
                                        column(
                                          tags$style(HTML("
      @import url('https://fonts.googleapis.com/css2?family=Yusei+Magic&display=swap');
      .slick-active {
  background: transparent;
}
.slick-track{
  background: transparent;
}
.slick-current{
  background: transparent;
}
.slick-slide{
  background: transparent;
}
.slick-cloned{
  background: transparent;
}
.slick-slide + .slick-cloned{
  background: transparent;
}
.slick-active + .slick-active {
  background: transparent;
}

.slick-active + .slick-active + .slick-active {
  background: transparent;
}")
                                        ),
                                          slickROutput("slickr", width="auto") %>% withSpinner(color="#0dc5c1"), width = 8),
                                      column(div(),width = 2),
                                      )),
                             
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
                                      tabPanel("GO search",br(),sidebarLayout(
                                        sidebarPanel(
                                          textInput("goE",p("Gene Onthology ID:",style="color:white; text-align:center")),
                                        ),mainPanel(
                                          DT::dataTableOutput("searchGO_E") %>% withSpinner(color="#0dc5c1")
                                        )
                                      ))
                                      
                                      ),br(),br(),
                                      
                             ),
                             tabPanel(div("Roots",style="text-align: center;"),
                                      tabsetPanel(
                                        tabPanel("Expression View",br(),
                                      sidebarLayout(
                                        
                                        sidebarPanel(
                                          useShinyjs(),
                                          textInput("seqID_R",p("ID gene database:",style="color:white; text-align:center")),
                                          imageOutput("legendtHeatR",height = "150px"),br(),
                                          div(tags$img(src="Root_legend.png",width="310px",height="320px")),
                                          
                                          width = 3
                                        ),
                                        mainPanel(
                                          fluidRow(
                                            uiOutput("descriptionSeq_R"),br()
                                          ),
                                          fluidRow(
                                            
                                            column(fluidRow(
                                              h5(p("Ammonium condition (3mM)",style="color:white;text-align:center"),style="background-color:black;border-radius: 8px"),
                                              imageOutput("heatRoots_N")),width = 4,style="text-align: center;"),
  
                                            column(fluidRow(
                                              useShinyjs(),
                                              h5(p("Differential expression",style="color:white;text-align:center"),style="background-color:black;border-radius: 8px"),
                                              tableOutput("diffRootTable"),br(),
                                              tableOutput("tableGOR"),br(),
                                              plotOutput("barplotRoots"), align="center"),br(),br(),
                                              width = 4,style="text-align: center;"),
                                            column(fluidRow(
                                              h5(p("Control condition",style="color:white;text-align:center"),style="background-color:black;border-radius: 8px"),
                                              imageOutput("heatRoots_C")),width = 4,style="text-align: center;")
                                          )
                                        ))),tabPanel("Correlation",br(),
                                                     sidebarPanel(
                                                       actionButton("buttonCorR", "Generate correlation table"),
                                                       width = 3
                                                     ),mainPanel(
                                                       DT::dataTableOutput("tableCorR") %>% withSpinner(color="#0dc5c1")
                                                     ))
                                      ,tabPanel("GO search",br(),
                                                sidebarPanel(
                                                  textInput("goR",p("Gene Onthology ID:",style="color:white; text-align:center")),
                                                ),mainPanel(
                                                  DT::dataTableOutput("searchGO_R") %>% withSpinner(color="#0dc5c1")
                                                )
                                      )
                                      # ,tabPanel("BLAST search",br(),
                                      #            sidebarPanel(
                                      #              textAreaInput("sequenceNUCL_r",p("Nucleotide sequences:",style="color:white; text-align:center")),width = 12
                                      #            ),mainPanel(br(),
                                      #              DT::dataTableOutput("bastSearch") %>% withSpinner(color="#0dc5c1"),
                                      #              br(),
                                      #            )
                                      # )
                                      )
                                      
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
                                        tabPanel("GO search",br(),sidebarLayout(
                                          sidebarPanel(
                                            textInput("goN",p("Gene Onthology ID:",style="color:white; text-align:center")),
                                          ),mainPanel(
                                            DT::dataTableOutput("searchGO_N") %>% withSpinner(color="#0dc5c1")
                                          )
                                        ))
                                        
                                      ),br(),br(),
                             ),
                             tabPanel("Step 4",
                                      
                                      fluidRow(column(width=2),
                                               column(
                                                   h4(p("Final model",style="color:black;text-align:center")),
                                                   width=8,style="background-color:lavender;border-radius: 10px")),
                                      br(),
                                      
                                      fluidRow(column(width=2, icon("hand-point-right","fa-5x"),align="center"),
                                               column(
                                                   p("Now we are going to build the final model, for this you will have to 
                                     select the independent variables you want to include in the model
                                     and especially select for which of them you want to include some non-linearity 
                                     in the model (this is related to the power transformations made to the 
                                     variables independent in a previous step)",style="color:black;text-align:justify"),
                                                   br(),
                                                   p("In this step it is very important to achieve a high adjusted coefficient of 
                                     determination and also to make the parameters of the model statistically significant, 
                                     for that reason we are going to test the following hypothesis:",style="color:black;text-align:justify"),
                                                   p('$$H_0: ~ \\beta_i = 0$$',style="color:black;border:1px solid black;background-color:white"),
                                                   width=8,style="background-color:lavender;border-radius: 10px")),
                                      br(),
                                      fluidRow(column(width=2),
                                               column(
                                                   p("Go ahead! let's build our final model",style="color:black; text-align:center"),width = 8,style="background-color:papayawhip; border-radius: 10px"
                                               )),
                                      hr(),
                                      sidebarLayout(
                                          
                                          sidebarPanel(
                                              
                                              
                                              fluidRow(column(
                                                  checkboxGroupInput("variablesincluidas",p("Select the independent variables you would like to include:",style="color:coral;text-align:justify"),choices = c("Projected population"=1,"Thefts"=3,"Traffic accidents"=5,"Homicides"=7,"School deserters"=9,"Sports venues"=11,"Extortions"=13),
                                                                     selected = c("Projected population"=1,"Thefts"=3,"Traffic accidents"=5,"Homicides"=7,"School deserters"=9,"Sports venues"=11,"Extortions"=13)),width = 6),
                                                  column(
                                                      checkboxGroupInput("incluirtrans",p("Do you want to include the transformation for this variable?",style="color:navy;text-align:justify"),choices = c("Include"=2,"Include"=4,"Include"=6,"Include"=8,"Include"=10,"Include"=12,"Include"=14)),
                                                      width=6)),
                                              br(),
                                              uiOutput("Anothermessage"),
                                              br(),
                                              hr(),
                                              p("However, you could have included some variable that is not really important for the model, so we are going to use another criterion (stepwise regression) to complement your choice",style="color:black;text-align:justify"),
                                              br(),
                                              selectInput("direccion",p("Choose the selection direction",style="color:dodgerblue"),choices = c("Backward"="backward","Forward"="forward")),
                                              selectInput("criterio",p("Criterion to be used",style="color:dodgerblue"),choices = c("AIC","BIC")),
                                              p("Read more about stepwise regression here →", a(href="https://en.wikipedia.org/wiki/Stepwise_regression",icon("wikipedia-w"),target="_blank"),style="color:black;text-align:justify")
                                              
                                              
                                          ),
                                          mainPanel(
                                              
                                              h3(p(strong('Model summary',style="color:salmon")),align="center"),
                                              
                                              tags$head(tags$style("#finalmodel{height: 400px;width:auto;border: 1px solid black; background-color: lavender}")),
                                              tags$head(tags$style("#ModeloBack{height: 400px;width:auto;border: 1px solid black; background-color: lavender}")),
                                              
                                              fluidRow(column(verbatimTextOutput("finalmodel"),width = 8),
                                                       
                                                       column(
                                                           uiOutput("Significancy2"),
                                                           br(),
                                                           uiOutput("FinalAlarma")
                                                           
                                                           ,width = 4)),
                                              hr(),
                                              h3(p(strong('Final model',style="color:salmon")),align="center"),
                                              fluidRow(column(verbatimTextOutput("ModeloBack"),width=8),
                                                       column(uiOutput("Determinacionfinal"),width=4))
                                              
                                              
                                          )
                                          
                                      )),
                             tabPanel("Step 5",
                                      
                                      fluidRow(column(width=2),
                                               column(
                                                   h4(p("Assumptions for model residuals",style="color:black;text-align:center")),
                                                   width=8,style="background-color:lavender;border-radius: 10px")),
                                      br(),
                                      
                                      fluidRow(column(width=2, icon("hand-point-right","fa-5x"),align="center"),
                                               column(
                                                   p("In addition to the assumptions that we have already discussed so far, our final model has three other assumptions associated, specifically in the residuals. These are:",style="color:black;text-align:justify"),
                                                   p("1. Normality:", '\\( H_0:~e_i \\sim Normal(\\mu,\\sigma) \\)',br(),br(
                                                       "2. Homoscedasticity", '\\( H_0: \\sigma_{e_i}^2 = \\sigma^2 \\)'),br(
                                                           "3. No auto-correlation", '\\( H_0: cor(e_i,e_{i+k}) = 0 \\)'),style="color:black;text-align:center;background-color:white;padding:15px;border:1px solid black"),br(),
                                                   
                                                   width=8,style="background-color:lavender;border-radius: 10px")),
                                      br(),
                                      fluidRow(column(width=2),
                                               column(
                                                   p("Let's do it. You are going to find some graphical and analytical tests in order to conclude about the previous hypothesis",style="color:black;text-align:center"),
                                                   width=8,style="background-color:papayawhip;border-radius: 10px")
                                      ),
                                      
                                      hr(),
                                      
                                      
                                      tabsetPanel(
                                          
                                          tabPanel("Normality",
                                                   
                                                   br(),
                                                   navlistPanel(widths=c(2,9),fluid = T,well = T,
                                                                
                                                                tabPanel("Graphical test",
                                                                         
                                                                         fluidRow(column(width=1),column(  
                                                                             fluidRow(
                                                                                 column(br(),plotOutput("Histograma2"),br(),width=4,style="border:1px solid black"),
                                                                                 column(br(),plotOutput("Boxplot2"),br(),width=4,style="border: 1px solid black;border-left: none"),
                                                                                 column(br(),plotOutput("qqPlot2"),br(),width=4,style="border:1px solid black;border-left:none")
                                                                                 
                                                                             ), width=11))),
                                                                
                                                                
                                                                tabPanel("Analytical test",
                                                                         
                                                                         selectInput("PruebaAnalitica2",p("Please select the test you want to try:",style="color:black; text-align:center"),choices=c("Shapiro-Wilk"=1,"Anderson-Darling"=2,"Cramer-von Mises"=3,"Kolmogorov-Smirnov"=4,"Jarque-Bera"=5)),
                                                                         uiOutput("ReadMore2"),
                                                                         
                                                                         
                                                                         fluidRow(
                                                                             
                                                                             tags$head(tags$style("#Conclusion12{color: navy;
                                                                         font-size: 15px;
                                                                         font-style: italic;
                                                                         font-weight: bold;
                                                                         text-align: center
                                                                         }")),
                                                                             tags$head(tags$style("#Prueba2{height: 155px; border: 1px solid black; background-color: lavender}")),
                                                                             column(verbatimTextOutput("Prueba2"),
                                                                                    br(),width = 6),
                                                                             column(br(),
                                                                                    p("Remember we are looking for a p-value greater than 0.05 (for a confidence level of 95%), so:",style="color:black"),
                                                                                    br(),
                                                                                    textOutput("Conclusion12"),
                                                                                    br(),width = 6,style="background-color:lavender;border-left:8px solid blue")
                                                                             
                                                                         )
                                                                         
                                                                ))),
                                          tabPanel("Constant variance",
                                                   
                                                   br(),
                                                   navlistPanel(widths=c(2,10),
                                                                
                                                                tabPanel("Graphical test",
                                                                         
                                                                         fluidRow(column(width=1),column(  
                                                                             fluidRow(
                                                                                 column(br(),plotOutput("FittedVsResiduals"),br(),width=11,style="border:1px solid black")
                                                                                 
                                                                             ), width=11))),
                                                                tabPanel("Analytical test",
                                                                         
                                                                         selectInput("PruebaAnalitica3",p("Please select the test you want to try:",style="color:black; text-align:center"),choices=c("Breusch-Pagan"=1,"Non-Constant Variance test"=2)),
                                                                         uiOutput("ReadMore3"),
                                                                         
                                                                         
                                                                         fluidRow(
                                                                             
                                                                             tags$head(tags$style("#Conclusion13{color: navy;
                                                                                         font-size: 15px;
                                                                                         font-style: italic;
                                                                                         font-weight: bold;
                                                                                         text-align: center
                                                                                         }")),
                                                                             tags$head(tags$style("#Prueba3{height: 140px; border: 1px solid black; background-color: lavender}")),
                                                                             column(verbatimTextOutput("Prueba3"),
                                                                                    br(),width = 6),
                                                                             column(br(),
                                                                                    p("Remember we are looking for a p-value greater than 0.05 (for a confidence level of 95%), so:",style="color:black"),
                                                                                    br(),
                                                                                    textOutput("Conclusion13"),
                                                                                    br(),width = 6,style="background-color:lavender;border-left:8px solid blue")
                                                                             
                                                                         )
                                                                         
                                                                ))),
                                          tabPanel("No auto-correlation",
                                                   
                                                   br(),
                                                   navlistPanel(widths=c(2,10),
                                                                
                                                                tabPanel("Graphical test",
                                                                         
                                                                         fluidRow(column(width=1),column(  
                                                                             fluidRow(
                                                                                 column(br(),plotOutput("ACF"),br(),width=4,style="border:1px solid black"),
                                                                                 column(br(),plotOutput("PACF"),br(),width=4,style="border: 1px solid black;border-left: none"),
                                                                                 column(br(),plotOutput("ResVsIndex"),br(),width=4,style="border:1px solid black;border-left:none")
                                                                                 
                                                                             ), width=11))
                                                                         
                                                                ),
                                                                tabPanel("Analytical test",
                                                                         
                                                                         selectInput("PruebaAnalitica4",p("Please select the test you want to try:",style="color:black; text-align:center"),choices=c("Durbin-Watson"=1,"Breusch-Godfrey"=2)),
                                                                         uiOutput("ReadMore4"),
                                                                         
                                                                         
                                                                         fluidRow(
                                                                             
                                                                             tags$head(tags$style("#Conclusion14{color: navy;
                                                                                font-size: 15px;
                                                                                font-style: italic;
                                                                                font-weight: bold;
                                                                                text-align: center
                                                                                }")),
                                                                             tags$head(tags$style("#Prueba4{height: 155px; border: 1px solid black; background-color: lavender}")),
                                                                             column(verbatimTextOutput("Prueba4"),
                                                                                    br(),width = 6),
                                                                             column(br(),
                                                                                    p("Remember we are looking for a p-value greater than 0.05 (for a confidence level of 95%), so:",style="color:black"),
                                                                                    br(),
                                                                                    textOutput("Conclusion14"),
                                                                                    br(),width = 6,style="background-color:lavender;border-left:8px solid blue")
                                                                             
                                                                         )
                                                                         
                                                                )))
                                          
                                      )),
                             tabPanel(icon("pencil-alt"),
                                      
                                      fluidRow(column(width=2),
                                               column(
                                                   h4(p("Let's review what we have learned",style="color:black;text-align:center")),
                                                   width=8,style="background-color:#BFF7BB;border-radius: 10px")),
                                      
                                      br(),
                                      fluidRow(column(width=2,icon("hand-point-right","fa-5x"),align="center"),
                                               column(
                                                   
                                                   p("This is the final section, here you will find a set of questions for 
                                     everything you have learned so far, there are 3 levels of difficulty.
                                      Try to answer correctly, if you make mistakes it 
                                     does not matter, after sending your answers you will find feedback regardless 
                                     of whether you did it right or wrong. You will also find a glossary in case 
                                     we have not dealt with a concept in a theoretical way."),
                                                   
                                                   width=8,style="background-color:#BFF7BB;border-radius: 10px")),
                                      
                                      br(),
                                      fluidRow(column(width=2),
                                               column(
                                                   p("Go for it!",style="color:black;text-align:center"),
                                                   width=8,style="background-color:lavender;border-radius: 10px")),
                                      hr(),
                                      navlistPanel(widths=c(2,10),
                                                   tabPanel("Level 1",
                                                            
                                                            column(width=12,
                                                                   h3("Level 1 questions (Theory)"),
                                                                   p("1. The following function represents a deterministic relationship between a response variable and a regression variable:"),
                                                                   fluidRow(column(width=4,
                                                                                   radioButtons("Question1","",choices=c('a. \\( Y = \\beta_0 + \\beta_1X \\)'=1,'b. \\( Y = \\beta_0 + \\beta_1X + e\\)'=2,'c. \\( \\hat Y = \\hat \\beta_0 + \\hat \\beta_1X \\)'=3,'d. \\( \\hat Y = \\hat \\beta_0 \\)'=4),selected = "_None")),
                                                                            column(width=8,uiOutput("Answer1"))),
                                                                   br(),
                                                                   p("2. With a regression model, statistical inference can be made as long as:"),
                                                                   fluidRow(column(width=4,
                                                                                   radioButtons("Question2","",choices=c('a. \\( R^2_{adj} \\Rightarrow 1 \\)'=1,'b. The model is statistically valid'=2,'c. a and b are corrects'=3, 'd. None of the above'=4),selected = "_none")),
                                                                            column(width=8,uiOutput("Answer2"))),
                                                                   br(),
                                                                   p("3. In order to study the safety problem in rural areas, it is decided to study the relationship between 
                                                school deserters and homicides. We run a simple linear regression model between these two variables, 
                                                we also validate the model. From the next output, what can we conclude?"),
                                                                   fluidRow(column(width=4,p(strong("Durbin-Watson test"),br(strong("data: model")),br(strong("DW = 2.8462, p-value = 0.008587")),style="font-family:courier;background-color:#F6F6F6;padding: 10px"))),
                                                                   fluidRow(column(width=4,
                                                                                   radioButtons("Question3","",choices=c('a. \\( Corr(e_i,e_j) \\neq 0 \\) for at least one \\(i \\neq j \\)'=1,'b. \\( \\sigma^2 \\) is not constant'=2,'c. \\( Corr(e_i,e_j) \\neq 0 ~  \\forall ~ i \\neq j \\)'=3, 'd. We cannot conclude, it is necessary other tests'=4),selected = "_none")),
                                                                            column(width=8,uiOutput("Answer3"))),
                                                                   br(),
                                                                   p("4. Data were collected on the number of sport venues in rural areas and the projected population. This information is useful to know the investment in sport and culture of each municipality based on its population. 
                                                Based on the following result of the adjustment of the linear regression model, conclude:"),
                                                                   fluidRow(column(width=4,p(strong("Analysis of Variance Table"),br(strong("Response: Sports venues")),br(strong("............. Df.. Sum Sq")),strong("ProjPopulation 1 . 8861835"),br(strong("Residuals ... 33 . 308965")),style="font-family:courier;background-color:#F6F6F6;padding: 10px"))),
                                                                   fluidRow(column(width=4,
                                                                                   radioButtons("Question4","",choices=c('a. The correlation between the number of sports venues and the projected population is 0.966'=1,'b. The correlation between the number of sports venues and the projected population is 0.983 '=2,'c. The correlation between the number of sports venues and the projected population is 0.287'=3, 'd. The correlation between the number of sports venues and the projected population is 0.536'=4),selected = "_none")),
                                                                            column(width=8,uiOutput("Answer4"))),
                                                                   br(),
                                                                   p("5. How would you interpret the property of the regression line that says that: \\( \\sum_{i=1}^n y_i = \\sum_{i=1}^n \\hat y_i \\)?"),
                                                                   fluidRow(column(width=4,
                                                                                   radioButtons("Question5","",choices=c('a. Neither the global adjustment nor the adjustment in each point of the response variable are perfect'=1,'b. The adjustment in each point is perfect, although the global adjustment of the response variable is not'=2,'c. The global adjustment and the adjustment at each point of the response variable are perfect'=3, 'd. The global adjustment of the response variable is perfect, although the adjustment at each point is not. '=4),selected = "_none")),
                                                                            column(width=8,uiOutput("Answer5"))),
                                                                   br(),
                                                                   p("6. What is the difference between the R squared and the adjusted R squared."),
                                                                   fluidRow(column(width=4,
                                                                                   radioButtons("Question6","",choices=c('a. The adjusted R squared increases the more variables there are while the R squared penalizes the inclusion of variables if they do not provide information.'=1,'b. The R squared increases the more variables there are while the adjusted R squared penalizes the inclusion of variables if they do not provide information.'=2,'c. The adjusted R squared represents the quality of the adjustment while the squared R indicates how good the predictions will be.'=3, 'd. There is no practical difference, both are measures of the performance of a model but were developed by different authors.'=4),selected = "_none")),
                                                                            column(width=8,uiOutput("Answer6"))),
                                                                   
                                                                   
                                                                   style="border: 10px solid transparent;border-image: url(border.png) 30 round")
                                                            
                                                   ),
                                                   tabPanel("Level 2",
                                                            
                                                            column(width=12,
                                                                   h3("Level 2 questions (Graphics analysis)"),
                                                                   p("1. The following diagram represents a scatter plot between a response variable and a regression variable. "),
                                                                   tags$img(src="Scatter1.png",width="300px",height="260px"),
                                                                   p("From the graph above one could conclude that:"),
                                                                   fluidRow(column(width=4,
                                                                                   radioButtons("Question7","",choices=c('a. The relationship between the response variable and the independent one is weak and it is not reasonable to adjust a linear regression model. '=1,'b. The relationship between the response variable and the independent one is non-linear and it is reasonable to adjust a linear regression model. '=2,'c. The relationship between the response variable and the independent one is strong and it is not reasonable to adjust a linear regression model.'=3,'d. The relationship between the response variable and the independent one is strong and it is reasonable to adjust a linear regression model. '=4),selected = "_None")),
                                                                            column(width=8,uiOutput("Answer7"))),
                                                                   br(),
                                                                   p("2. The following graph represents the residuals VS the adjusted values of a linear regression model."),
                                                                   tags$img(src="Scatter2.png",width="300px",height="260px"),
                                                                   p("From the graph above one could conclude that:"),
                                                                   fluidRow(column(width=4,
                                                                                   radioButtons("Question8","",choices=c('a. The assumption of normality of errors is fulfilled for this model and it is not necessary to transform Y.'=1,'b. The assumption of constant variance of the errors is not fulfilled for this model and it is necessary to transform Y.'=2,'c. The assumption of constant variance of the errors is fulfilled for this model and it is not necessary to transform to Y.'=3, 'd. The assumption of normality of errors is not fulfilled for this model and it is necessary to transform Y.'=4),selected = "_None")),
                                                                            column(width=8,uiOutput("Answer8"))),
                                                                   br(),
                                                                   p("3. The following graph represents the qqNorm of the residuals of a linear regression model."),
                                                                   tags$img(src="qqNorm.png",width="300px",height="260px"),
                                                                   p("From the graph above one could conclude that:"),
                                                                   fluidRow(column(width=4,
                                                                                   radioButtons("Question9","",choices=c('a. The assumption of normality of errors is fulfilled in this model.'=1,'b. The assumption of constant variance of errors is not fulfilled in this model. '=2,'c. The assumption of constant variance of the errors is fulfilled for this model.'=3, 'd. The assumption of normality of errors is not fulfilled for this model.'=4),selected = "_None")),
                                                                            column(width=8,uiOutput("Answer9"))),
                                                                   br(),
                                                                   p("4. The following graph represents a scatter plot between two variables"),
                                                                   tags$img(src="Scatter3.png",width="300px",height="260px"),
                                                                   p("From the graph above one could conclude that:"),
                                                                   fluidRow(column(width=4,
                                                                                   radioButtons("Question10","",choices=c('a. The relationship is linear and positive'=1,'b. The relationship is non-linear and negative '=2,'c. The relationship is non-linear and positive '=3, 'd. There is no relationship'=4),selected = "_None")),
                                                                            column(width=8,uiOutput("Answer10"))),
                                                                   br(),
                                                                   p("5. The following graph represents a scatter plot between two variables"),
                                                                   tags$img(src="Scatter4.png",width="300px",height="260px"),
                                                                   p("From the graph above one could conclude that:"),
                                                                   fluidRow(column(width=4,
                                                                                   radioButtons("Question11","",choices=c('a. The relationship is linear and positive'=1,'b. The relationship is non-linear and negative '=2,'c. The relationship is non-linear and positive '=3, 'd. There is no relationship'=4),selected = "_None")),
                                                                            column(width=8,uiOutput("Answer11"))),
                                                                   style="border: 10px solid transparent;border-image: url(border3.png) 30 round"
                                                            )
                                                            
                                                   ),
                                                   tabPanel("Level 3",
                                                            
                                                            column(width=12,
                                                                   h3("Level 3 questions (Problem analysis)"),
                                                                   p("1. What is the interpretation of the coefficient \\( \\hat \\beta_3 \\) corresponding to the variable traffic accidents in the following model."),
                                                                   tags$img(src="summary1.png",width="450px",height="240px"),
                                                                   fluidRow(column(width=4,
                                                                                   radioButtons("Question12","",choices=c('a. 67.668% of personal injuries originate in a traffic accident. '=1,'b. 0.67668% of personal injuries originate in a traffic accident. '=2,'c.For each traffic accident, the number of reported personal injuries increases by 0.6768, as long as the other variables are zero or constant. '=3,'d. None of the above'=4),selected = "_None")),
                                                                            column(width=8,br(),uiOutput("Answer12"))),
                                                                   br(),
                                                                   p("2. The following image presents the VIFS of each of the independent variables within the model, with that information we can determine that:"),
                                                                   tags$img(src="summary2.png",width="450px",height="120px"),
                                                                   fluidRow(column(width=4,
                                                                                   radioButtons("Question13","",choices=c('a. The variable "Sports venues" must be eliminated from the model because it has the lowest VIF.'=1,'b. Projected population is the variable that presents multicollinearity and must be irremediably eliminated from the model.'=2,'c. There are two variables with multicollinearity problems, so one of them should be removed from the model or grouped in an indicator.'=3, 'd. The range of values of the VIF is not so large, so the estimation of the model would not be affected by multicollinearity problems.'=4),selected = "_None")),
                                                                            column(width=8,br(),uiOutput("Answer13"))),
                                                                   br(),
                                                                   p("3. According to the following result it can be stated that: "),
                                                                   tags$img(src="summary3.png",width="450px",height="90px"),
                                                                   fluidRow(column(width=4,
                                                                                   radioButtons("Question14","",choices=c('a. The variables are linearly related since the coefficient associated with the regressor variable is significant.'=1,'b. The variables are not related linearly since the coefficient associated with the regressor variable is significant. '=2,'c. The variables are related linearly since the coefficient associated with the regressor variable is not significant.'=3, 'd. The variables are not related linearly since the coefficient associated with the regressor variable is not significant.'=4),selected = "_None")),
                                                                            column(width=8,br(),uiOutput("Answer14"))),
                                                                   br(),
                                                                   p("4. From the following output we can conclude that: "),
                                                                   tags$img(src="summary4.png",width="450px",height="90px"),
                                                                   fluidRow(column(width=4,
                                                                                   radioButtons("Question15","",choices=c('a. The model has good goodness of fit and explains 42.56% the variability of the response variable around its mean.'=1,'b. The model does not have a good godness of fit since it only explains 0.7303% of the variability of the response variable around its mean.'=2,'c. The model has a good goodness of fit and explains 71.32% the variability of the response variable around its mean.'=3, 'd. The model has no statistical significance because it has a very small p-value.'=4),selected = "_None")),
                                                                            column(width=8,br(),uiOutput("Answer15"))),
                                                                   br(),
                                                                   p("5. From the following model we can conclude that:"),
                                                                   tags$img(src="summary5.png",width="450px",height="300px"),
                                                                   fluidRow(column(width=4,
                                                                                   radioButtons("Question16","",choices=c('a. The variables contained in the database completely explain the safety problem in the rural areas of Antioquia.'=1,'b. There are variables that probably do not help to explain the safety problem in rural areas of Antioquia such as traffic accidents and homicides.'=2,'c. There are variables that do not have a linear relationship with the response variable, so other types of relationships should be studied before stating that they do not help explain the safety problem.'=3, 'd.  It should be considered to make a model with only two explanatory variables, in this case projected population and thefts. '=4),selected = "_None")),
                                                                            column(width=8,br(),uiOutput("Answer16"))),
                                                                   style="border: 10px solid transparent;border-image: url(border2.png) 30 round"
                                                            )
                                                            
                                                   ),
                                                   tabPanel(icon("book"),
                                                            
                                                            column(width=12,
                                                                   h3("Glossary"),
                                                                   br(),
                                                                   p(strong("Deterministic relationship:"),"Two variables are related in a deterministic way as long as exist an exact relationship between them. e.g you will pay $8 for each hour playing a video game, it would be always the same price and the more hours you play the more you will pay, increasing exactly 8$ per hour."),
                                                                   br(),
                                                                   p(strong("Statistical relationship:"),"This type of relationship means that one variable influences the behavior of another. However, there is no exact relationship between them because there will always be random factors that also influence. e.g The more hours you exercise, the more weight you lose, however, you do not know in what proportion, because it is different in each organism."),
                                                                   br(),
                                                                   p(strong("Coefficient of determination or R squared:"), "Represents the proportion of the variability of the response variable that is predicted from the independent variable, this measure increases with every variable include in the model. This is a measure of how well the regression predictions approximate the real data points."),
                                                                   br(),
                                                                   p(strong("Adjusted Coefficient of determination:"), "this is also a measure of godness of fit of the model, however this measure penalizes the inclusion of variables if they do not provide additional information."),
                                                                   br(),
                                                                   p(strong("Correlation coefficient"), "is a measure of the linear relationship between two variables, taking values form -1 to +1, being -1 a perfect negative linear relationship and +1 a perfect positive linear relationship. This coefficient is the root squared of R squared"),
                                                                   style="border: 10px solid transparent;border-image: url(border4.png) 30 round"
                                                            )
                                                            
                                                            
                                                   ))
                                      
                             )
                             
                  ), class='flex-center'
                  
                  
))