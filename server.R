library(shiny)
library(shinyjs)
library(shinythemes)
library(shinycssloaders)
library(DT)
library(ggplot2)
library(car)
library(nortest)
library(tseries)
library(RcmdrMisc)
library(lmtest)

datos <-read.csv("www/dataset.csv",dec = ",")
datosEmbriones <- read.csv("www/CPM_embriones_media_filter.txt", sep = "\t", row.names = 1 )
tablaSeqGo <- read.csv("www/seq_gos.txt", sep = "\t",col.names = c("seqName","GoID"), stringsAsFactors = F)
tablaSeqDes <- read.csv("www/seq_description.txt", sep = "\t",col.names = c("seq","description"), row.names = 1,stringsAsFactors = F)
RM_data <- read.csv("www/CPM_raizMicrodiseccion_RM_mean.txt", sep = "\t", row.names = 1, stringsAsFactors = F )
tablaSeqDesRoots <- read.csv("www/seq_description_raiz.txt", sep = "\t",col.names = c("seq","description"), row.names = 1,stringsAsFactors = F)
tablaDiffRoots <-read.csv("www/microdiseccion_diffExpression.txt", sep = "\t", col.names = c("seqName","Tissue","Condition","Log Fold Change"), stringsAsFactors = F )
tablaDiffEmbryos <-read.csv("www/embryo_differential.txt", sep = "\t", col.names = c("seqName","Comparative","Type","Log Fold Change"), stringsAsFactors = F )
tablaSeqGoR <- read.csv("www/seq_gos_raiz.txt", sep = "\t",col.names = c("seqName","GoID"),stringsAsFactors = F)
keggTable <- read.csv("www/kegg_total_transcriptome.txt", col.names = c("seqName","KeggID"),sep = "\t", stringsAsFactors = F )
expressionNeedles <- read.csv("www/Needles_CPM_mean.txt", sep = "\t", row.names = 1 )
tablaDiffNeedles <-read.csv("www/needles_differential.txt", sep = "\t", col.names = c("seqName","Condition","Log Fold Change"), stringsAsFactors = F )

source("Components/Comp_Tables.R")
source("Components/Comp_embryos.R")
source("Components/Comp_needle.R")
source("Components/Comp_barplot.R")
source("Components/Comp_corTable.R")
source("Components/Comp_roots.R")
source("Components/Comp_blast.R")
shinyServer(function(input, output) {
    
    output$RawData <- DT::renderDataTable(
        DT::datatable({
            datos
        },
        options = list(lengthMenu=list(c(5,15,20),c('5','15','20')),pageLength=10,
                       initComplete = JS(
                           "function(settings, json) {",
                           "$(this.api().table().header()).css({'background-color': '1c1b1b', 'color': '1c1b1b'});",
                           "}"),
                       columnDefs=list(list(className='dt-center',targets="_all"))
        ),
        filter = "top",
        selection = 'multiple',
        style = 'bootstrap',
        class = 'cell-border stripe',
        rownames = FALSE,
        colnames = c("Subregion","Municipality","Projected population","Thefts","Traffic accidents","Homicides","School deserters","Sports venues","Extortions","Personal injuries")
        ))
    
    variabletrans <- reactive({
        
        attach(datos)
        if(input$Transformacion == 0){Lesiones.t=log(LesionesPer)}else{Lesiones.t=LesionesPer^as.double(input$Transformacion)}
        if(input$Transformacion ==1){titulo="Personal Injuries"}else{titulo="Transformed Personal Injuries"}  
        Resultado <- cbind(Lesiones.t,titulo)
        
        
    })
    
    variabletrans2 <- reactive({
        
        attach(datos)
        if(input$Transformacion == 0){Lesiones.t=log(LesionesPer)}else{Lesiones.t=LesionesPer^as.double(input$Transformacion)}
        if(input$Transformacion ==1){titulo="Personal Injuries"}else{titulo="Transformed Personal Injuries"}  
        Lesiones.t
    })
    
    output$Histograma <- renderPlot({
        
        ggplot(NULL,aes(as.double(variabletrans()[,1])))+geom_histogram(bins=nclass.Sturges(as.double(variabletrans()[,1])),color="white",
                                                                        fill="seagreen1",aes(y=..density..),lwd=0.8)+geom_density(color="seagreen4",
                                                                                                                                  alpha=0.3,fill="seagreen4",lty=1)+
            labs(title = paste(variabletrans()[1,2], "\n histogram"),x=variabletrans()[1,2],y="Density")+
            theme(plot.title = element_text(color="navy", size=15, face="bold.italic",hjust=0.5),
                  axis.title.x = element_text(color="navy", size=13, face="bold"),
                  axis.title.y = element_text(color="navy", size=13, face="bold"))
        
        
    })
    
    output$Boxplot <- renderPlot({
        
        ggplot(NULL,aes(x=0,y=as.double(variabletrans()[,1])))+geom_boxplot(color="black",fill="skyblue",alpha=0.5)+ stat_summary(fun.y=mean, colour="darkred", geom="point",shape=18, size=3)+
            labs(title = paste(variabletrans()[1,2], "\n boxplot"),x="",y=variabletrans()[1,2])+
            theme(plot.title = element_text(color="navy", size=15, face="bold.italic",hjust=0.5),
                  axis.title.x = element_text(),
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.title.y = element_text(color="navy", size=13, face="bold"))
        
        
    })
    
    output$qqPlot <- renderPlot({
        
        par(font.main=4,font.lab=2,col.main="navy",col.lab="navy",cex.lab=1.1)
        qqPlot(as.double(variabletrans()[,1]),grid=F,xlab="",ylab="")
        u <-par("usr")
        rect(u[1], u[3], u[2], u[4], col="#EBE9E9", border=TRUE)
        grid(NULL,NULL,col="white",lty=1)
        par(new=TRUE)
        qqPlot(as.double(variabletrans()[,1]),col="coral",pch=16,id=T,lwd=1.9,col.lines = "black",grid = F, main = paste(variabletrans()[1,2], "\n Q-Q plot"),xlab="Normal quantiles",ylab=variabletrans()[1,2])
        box(col="white")
        
        
    })
    
#Blast search Roots
    # output$bastSearch <- DT::renderDataTable({
    #     DT::datatable(blast_raiz(input$sequenceNUCL_r),options = list(initComplete = JS(
    #         "function(settings, json) {",
    #         "$(this.api().table().header()).css({'background-color': #1c1b1b', 'color': '#1c1b1b'});",
    #         "}")),
    #         filter = "top",
    #         selection = 'multiple',
    #         style = 'bootstrap',
    #         class = 'cell-border stripe',
    #         rownames = FALSE,
    #     )
    # })
    
#Slider imagenes
    output$slickr <- renderSlickR({
        imgs <- list.files("www/fotosSlider/", pattern=".jpg", full.names = TRUE)
        
        
        slick <- slickR(imgs,height = 600,
               width = '100%')
        slick <- slick + settings(dots = TRUE, autoplay = TRUE, autoplaySpeed = 3000)
        slick
    })
    
    newCorTable <- eventReactive(input$buttonCorE, {
        if (input$seqID %in% row.names(datosEmbriones)){
            tableCor <- seqCorEmbryos(input$seqID)
            tableCor
        }
    })

    
    output$tableCorE <- DT::renderDataTable({
        table <- newCorTable()
        DT::datatable(table,extensions = 'Buttons',options = list(initComplete = JS(
            "function(settings, json) {",
            "$(this.api().table().header()).css({'background-color': '#1c1b1b', 'color': '#1c1b1b'});",
            "}"),dom = 'Bfrtip',pageLength = length(row.names(table)),buttons = c('copy', 'csv', 'excel', 'pdf', 'print')),
            filter = "top",
            selection = 'multiple',
            style = 'bootstrap',
            class = 'cell-border stripe',
            rownames = FALSE,
        )
    })
    
    newCorTableN <- eventReactive(input$buttonCorN, {
        if (input$seqID_N %in% row.names(expressionNeedles)){
            tableCor <- seqCorNeedles(input$seqID_N)
            tableCor
        }
    })
    
    
    output$tableCorN <- DT::renderDataTable({
        table <- newCorTableN()
        DT::datatable(table,extensions = 'Buttons',options = list(initComplete = JS(
            "function(settings, json) {",
            "$(this.api().table().header()).css({'background-color': '#1c1b1b', 'color': '#1c1b1b'});",
            "}"),dom = 'Bfrtip',pageLength = length(row.names(table)),buttons = c('copy', 'csv', 'excel', 'pdf', 'print')),
            filter = "top",
            selection = 'multiple',
            style = 'bootstrap',
            class = 'cell-border stripe',
            rownames = FALSE,
        )
    })
    
    newCorTableR <- eventReactive(input$buttonCorR, {
        tablaCorResR <- NULL
        if (input$seqID_R %in% row.names(RM_data)){
            tablaCorResR<- seqCorRoots(input$seqID_R)
            
        }
        tablaCorResR
    })
    
    output$tableCorR <- DT::renderDataTable({
        if (input$seqID_R %in% row.names(RM_data)){
            table <- newCorTableR()
            DT::datatable(table,extensions = 'Buttons',options = list(initComplete = JS(
                "function(settings, json) {",
                "$(this.api().table().header()).css({'background-color': '#1c1b1b', 'color': '#1c1b1b'});",
                "}"),dom = 'Bfrtip',pageLength = length(row.names(table)),buttons = c('copy', 'csv', 'excel', 'pdf', 'print')),
                filter = "top",
                selection = 'multiple',
                style = 'bootstrap',
                class = 'cell-border stripe',
                rownames = FALSE,
            )
        }
    })
    
    output$searchGO_E <- DT::renderDataTable({
        table <- searchGoSeq(input$goE)
        DT::datatable(table,extensions = 'Buttons',options = list(initComplete = JS(
            "function(settings, json) {",
            "$(this.api().table().header()).css({'background-color': '#1c1b1b', 'color': '#1c1b1b'});",
            "}"),dom = 'Bfrtip',pageLength = length(row.names(table)),buttons = c('copy', 'csv', 'excel', 'pdf', 'print')),
            filter = "top",
            selection = 'multiple',
            style = 'bootstrap',
            class = 'cell-border stripe',
            rownames = FALSE,
        )
    })

    output$searchGO_N <- DT::renderDataTable({
        table <- searchGoSeqN(input$goN)
        DT::datatable(table,extensions = 'Buttons',options = list(initComplete = JS(
            "function(settings, json) {",
            "$(this.api().table().header()).css({'background-color': '#1c1b1b', 'color': '#1c1b1b'});",
            "}"),dom = 'Bfrtip',pageLength = length(row.names(table)),buttons = c('copy', 'csv', 'excel', 'pdf', 'print')),
            filter = "top",
            selection = 'multiple',
            style = 'bootstrap',
            class = 'cell-border stripe',
            rownames = FALSE,
        )
    })
    
    output$searchGO_R <- DT::renderDataTable({
        table<-searchGoSeqR(input$goR)
        DT::datatable(table,extensions = 'Buttons',options = list(initComplete = JS(
            "function(settings, json) {",
            "$(this.api().table().header()).css({'background-color': '#1c1b1b', 'color': '#1c1b1b'});",
            "}"),dom = 'Bfrtip',pageLength = length(row.names(table)),buttons = c('copy', 'csv', 'excel', 'pdf', 'print')),
            filter = "top",
            selection = 'multiple',
            style = 'bootstrap',
            class = 'cell-border stripe',
            rownames = FALSE,
        )
    })
    
    output$tableKeggSeq <- renderTable({
        if (input$keggID %in% keggTable$KeggID) {
            keggTable[which(keggTable$KeggID == input$keggID & keggTable$seqName %in% row.names(datosEmbriones)),]
        }else{
            if (input$seqID %in% row.names(datosEmbriones) & input$seqID %in% keggTable$seqName) {
                keggTable[which(keggTable$seqName == input$seqID),]
            }
        }
        
    })
    output$tableKeggSeqN <- renderTable({
        if (input$keggIDN %in% keggTable$KeggID) {
            keggTable[which(keggTable$KeggID == input$keggIDN & keggTable$seqName %in% row.names(expressionNeedles)),]
        }else{
            if (input$seqID_N %in% row.names(expressionNeedles) & input$seqID_N %in% keggTable$seqName) {
                keggTable[which(keggTable$seqName == input$seqID_N),]
            }
        }
        
    })
    
    output$tableGO <- renderTable({
        if (input$seqID %in% row.names(datosEmbriones) &&
            input$seqID %in% tablaSeqGo$seqName) {
            goTableDes(input$seqID)
        }
        
    })
    output$tableGOR <- renderTable({
        if (input$seqID_R %in% row.names(RM_data) &&
            input$seqID_R %in% tablaSeqGoR$seqName) {
            goTableDesR(input$seqID_R)
        }
        
    })
    output$tableGON <- renderTable({
        if (input$seqID_N %in% row.names(expressionNeedles) &&
            input$seqID_N %in% tablaSeqGo$seqName) {
            goTableDesN(input$seqID_N)
        }
        
    })
    
    output$descriptionSeq <- renderUI({
        if (input$seqID %in% row.names(datosEmbriones)) {
        column(
            h5(p(tablaSeqDes[input$seqID,],style="color:white;text-align:center")),
            width=12,style="background-color:black;border-radius: 8px")
            }
    })
    output$descriptionSeq_R <- renderUI({
        if (input$seqID_R %in% row.names(RM_data)) {
            column(
                h5(p(tablaSeqDesRoots[input$seqID_R,],style="color:white;text-align:center")),
                width=12,style="background-color:black;border-radius: 8px")
        }
    })
    output$descriptionSeqN <- renderUI({
        if (input$seqID_N %in% row.names(expressionNeedles)) {
            column(
                h5(p(tablaSeqDes[input$seqID_N,],style="color:white;text-align:center")),
                width=12,style="background-color:black;border-radius: 8px")
        }
    })
#Barplot Embriones    
    hide("barplotEmbryos")
    output$barplotEmbryos <- renderPlot({
        res <- plot.new()
        if (input$seqID %in% row.names(datosEmbriones)) {
            res <- barplotEmbryosFunc(input$seqID)
        }
        res
    })

    observeEvent(input$seqID,{
        if (input$seqID %in% row.names(datosEmbriones)) {
            show("barplotEmbryos")
        }else{
            hide("barplotEmbryos")
        }
        
    })
#Barplot Raiz
    hide("barplotRoots")
    output$barplotRoots <- renderPlot({
        res <- plot.new()
        if (input$seqID_R %in% row.names(RM_data)) {
            res <- barplotRootsFunc(input$seqID_R)
        }
        res
    })
    
    observeEvent(input$seqID_R,{
        if (input$seqID_R %in% row.names(RM_data)) {
            show("barplotRoots")
        }else{
            hide("barplotRoots")
        }
        
    })
    #Barplot Needles
    hide("barplotNeedles")
    output$barplotNeedles<- renderPlot({
        res <- plot.new()
        if (input$seqID_N %in% row.names(expressionNeedles)) {
            res <- barplotNeedlesFunc(input$seqID_N)
        }
        res
    })
    
    observeEvent(input$seqID_N,{
        if (input$seqID_N %in% row.names(expressionNeedles)) {
            show("barplotNeedles")
        }else{
            hide("barplotNeedles")
        }
        
    })
#Legend plot Roots
    hide("legendtHeatR")
    output$legendtHeatR <- renderImage({
        if (input$seqID_R %in% row.names(RM_data)) {
            legendColorsR(input$seqID_R)
        }
        list(src="www/legendR.svg",width = "310px",
             height = "180px")
    },deleteFile = TRUE)
    
    observeEvent(input$seqID_R,{
        if (input$seqID_R %in% row.names(RM_data)) {
            show("legendtHeatR")
        }else{
            hide("legendtHeatR")
        }
        
    })

    #Legend plot Roots
    hide("legendtHeatE")
    output$legendtHeatE <- renderImage({
        if (input$seqID %in% row.names(datosEmbriones)) {
            legendColorsE(input$seqID)
        }
        list(src="www/legendE.svg",width = "310px",
             height = "180px")
    },deleteFile = TRUE)
    
    observeEvent(input$seqID,{
        if (input$seqID %in% row.names(datosEmbriones)) {
            show("legendtHeatE")
        }else{
            hide("legendtHeatE")
        }
        
    })    
    
    #Legend plot Needles
    hide("legendtHeatN")
    output$legendtHeatN <- renderImage({
        if (input$seqID_N %in% row.names(expressionNeedles)) {
            legendColorsN(input$seqID_N)
        }
        list(src="www/legendN.svg",width = "310px",
             height = "180px")
    },deleteFile = TRUE)
    
    observeEvent(input$seqID_N,{
        if (input$seqID_N %in% row.names(expressionNeedles)) {
            show("legendtHeatN")
        }else{
            hide("legendtHeatN")
        }
        
    })
    
    output$diffRootTable <- renderTable({
        if (input$seqID_R %in% tablaDiffRoots$seqName) {
            tableDiffExpressionRoots(input$seqID_R)
        }
        
    })
    output$tableDiffEmbryos <- renderTable({
        if (input$seqID %in% tablaDiffEmbryos$seqName) {
            tableDiffExpressionEmbryos(input$seqID)
        }
        
    })
    output$tableDiffNeedles <- renderTable({
        if (input$seqID_N %in% tablaDiffNeedles$seqName) {
            tableDiffExpressionNeedles(input$seqID_N)
        }
        
    })
    
    output$heatMapEmbryos <- renderText({
        
        heatMapEmbry(input$seqID)
            
    })

    #Roots Heatmaps
    hide("heatRoots_N")
    output$heatRoots_N <- renderImage({
        if (input$seqID_R %in% row.names(RM_data)) {
            heatMapRootsN(input$seqID_R)
        }
        list(src="www/raiz_heat_modN.svg",width = "280px",
             height = "660px",style="margin-right:-40px")
    },deleteFile = TRUE)
    
    observeEvent(input$seqID_R,{
        if (input$seqID_R %in% row.names(RM_data)) {
            show("heatRoots_N")
        }else{
            hide("heatRoots_N")
        }
        
    })
    
    hide("heatRoots_C")
    output$heatRoots_C <- renderImage({
        if (input$seqID_R %in% row.names(RM_data)) {
            heatMapRootsC(input$seqID_R)
        }
        list(src="www/raiz_heat_modC.svg",width = "280px",
             height = "660px",style="margin-right:-40px")
    },deleteFile = TRUE)
    
    observeEvent(input$seqID_R,{
        if (input$seqID_R %in% row.names(RM_data)) {
            show("heatRoots_C")
        }else{
            hide("heatRoots_C")
        }
        
    })

    #Needles heatmaps
    hide("heatNeedles_N")
    output$heatNeedles_N <- renderImage({
        if (input$seqID_N %in% row.names(expressionNeedles)) {
            heatMapNeedlesN(input$seqID_N)
        }
        list(src="www/needle_N.svg",width = "700px",
             height = "530px",style="margin-left:-100px;margin-top:-100px;margin-bottom:-50px")
    },deleteFile = TRUE)
    
    observeEvent(input$seqID_N,{
        if (input$seqID_N %in% row.names(expressionNeedles)) {
            show("heatNeedles_N")
        }else{
            hide("heatNeedles_N")
        }
        
    })
    
    hide("heatNeedles_M")
    output$heatNeedles_M <- renderImage({
        if (input$seqID_N %in% row.names(expressionNeedles)) {
            heatMapNeedlesM(input$seqID_N)
        }
        list(src="www/needle_M.svg",width = "700px",
             height = "530px",style="margin-left:-100px;margin-top:-100px;margin-bottom:-50px")
    },deleteFile = TRUE)
    
    observeEvent(input$seqID_N,{
        if (input$seqID_N %in% row.names(expressionNeedles)) {
            show("heatNeedles_M")
        }else{
            hide("heatNeedles_M")
        }
        
    })
    
    output$Conclusion1 <- renderText({
        
        if(testanalitico()$p.value < 0.05){mensaje="We reject the null hypothesis, the variable is not normally distributed"}else{mensaje="We keep the null hypothesis, the variable is normally distributed"}
        mensaje
        
    })
    
    output$ReadMore <- renderUI({
        
        
        if(input$PruebaAnalitica == 1){
            
            More <- p("Read more about this test here → ", a(href="https://en.wikipedia.org/wiki/Shapiro%E2%80%93Wilk_test", icon("wikipedia-w"),target="_blank"),style="color:black;text-align:center")
            
        }
        else
        {
            
            if(input$PruebaAnalitica == 2){
                
                More <- p("Read more about this test here → ", a(href="https://en.wikipedia.org/wiki/Anderson%E2%80%93Darling_test", icon("wikipedia-w"),target="_blank"),style="color:black;text-align:center")
                
            }
            else
            {
                
                if(input$PruebaAnalitica == 3){
                    
                    More <- p("Read more about this test here → ", a(href="https://en.wikipedia.org/wiki/Cram%C3%A9r%E2%80%93von_Mises_criterion", icon("wikipedia-w"),target="_blank"),style="color:black;text-align:center")
                    
                }
                else
                {
                    
                    if(input$PruebaAnalitica ==4){
                        
                        More <- p("Read more about this test here → ", a(href="https://en.wikipedia.org/wiki/Kolmogorov%E2%80%93Smirnov_test", icon("wikipedia-w"),target="_blank"),style="color:black;text-align:center")
                        
                        
                    }
                    else
                    {
                        
                        More <- p("Read more about this test here → ", a(href="https://en.wikipedia.org/wiki/Jarque%E2%80%93Bera_test", icon("wikipedia-w"),target="_blank"),style="color:black;text-align:center")
                        
                        
                    }
                    
                }
                
            }
            
        }
        
        More
        
    })
    
    Trans1 <- reactive({
        
        attach(datos)
        if(input$Transformacion1 == 0){ProjectedPopulation.t=log(ProjectedPopulation)}else{ProjectedPopulation.t=ProjectedPopulation^as.double(input$Transformacion1)}
        if(input$Transformacion1 ==1){titulo="Projected Population"}else{
            if(input$Transformacion1==0){titulo="Projected Population Logarithm"}else{titulo="Transformed Projected Population"}}  
        Resultado1 <- cbind(ProjectedPopulation.t,titulo)
        
    })
    
    Trans2 <- reactive({
        
        attach(datos)
        if(input$Transformacion2 == 0){Thefts.t=log(Thefts)}else{Thefts.t=Thefts^as.double(input$Transformacion2)}
        if(input$Transformacion2 ==1){titulo="Thefts"}else{
            if(input$Transformacion2==0){titulo="Thefts Logarithm"}else{titulo="Transformed Thefts"}}  
        Resultado2 <- cbind(Thefts.t,titulo)
        
    })
    
    Trans3 <- reactive({
        
        attach(datos)
        if(input$Transformacion3 == 0){TrafAccid.t=log(TrafAccid)}else{TrafAccid.t=TrafAccid^as.double(input$Transformacion3)}
        if(input$Transformacion3 ==1){titulo="Traffic accidents"}else{
            if(input$Transformacion3==0){titulo="Traffic accidents Logarithm"}else{titulo="Transformed Traffic accidents"}}  
        Resultado3 <- cbind(TrafAccid.t,titulo)
        
    })
    
    Trans4 <- reactive({
        
        attach(datos)
        if(input$Transformacion4 == 0){Homicides.t=log(Homicides)}else{Homicides.t=Homicides^as.double(input$Transformacion4)}
        if(input$Transformacion4 ==1){titulo="Homicides"}else{
            if(input$Transformacion4 ==0){titulo="Homicides Logarithm"}else{titulo="Transformed Homicides"}}  
        Resultado4 <- cbind(Homicides.t,titulo)
        
    })
    
    Trans5 <- reactive({
        
        attach(datos)
        if(input$Transformacion5 == 0){SchoolDes.t=log(SchoolDes)}else{SchoolDes.t=SchoolDes^as.double(input$Transformacion5)}
        if(input$Transformacion5 ==1){titulo="School deserters"}else{
            if(input$Transformacion5 ==0){titulo="School deserters Logarithm"}else{titulo="Transformed School deserters"}}  
        Resultado5 <- cbind(SchoolDes.t,titulo)
        
    })
    
    Trans6 <- reactive({
        
        attach(datos)
        if(input$Transformacion6 == 0){SportsScenari.t=log(SportsScenari)}else{SportsScenari.t=SportsScenari^as.double(input$Transformacion6)}
        if(input$Transformacion6 ==1){titulo="Sports venues"}else{
            if(input$Transformacion6 ==0){titulo="Sports venues Logarithm"}else{titulo="Transformed Sports venues"}}  
        Resultado6 <- cbind(SportsScenari.t,titulo)
        
    })
    
    Trans7 <- reactive({
        
        attach(datos)
        if(input$Transformacion7 == 0){Extortions.t=log(Extortions)}else{Extortions.t=Extortions^as.double(input$Transformacion7)}
        if(input$Transformacion7 ==1){titulo="Extortions"}else{
            if(input$Transformacion7 ==0){titulo="Extortions Logarithm"}else{titulo="Transformed Extortions"}}  
        Resultado7 <- cbind(Extortions.t,titulo)
        
    })
    
    output$Dispersion1 <- renderPlot({
        
        ggplot(NULL,aes(x=as.double(Trans1()[,1]),y=as.double(variabletrans()[,1])))+
            geom_point(shape=18,color="blue",size=3)+
            geom_smooth(method = lm,linetype="dashed",color="black",fill="seagreen3")+
            labs(title = paste("\n",variabletrans()[1,2], "vs", Trans1()[1,2], "\n"),x=paste(Trans1()[1,2]),y=paste(variabletrans()[1,2],"\n"))+
            theme(plot.title = element_text(color="navy", size=18, face="bold.italic",hjust=0.5),
                  axis.title.x = element_text(color="navy", size=15, face="bold"),
                  axis.text.x = element_text(size=13),
                  axis.title.y = element_text(color="navy", size=15, face="bold"),
                  axis.text.y = element_text(size = 13))
        
        
    })
    
    output$correlacion1 <- renderText({
        
        corr <- cor(as.double(Trans1()[,1]),as.double(variabletrans()[,1]))
        
        mensaje <- paste("Cor = ", corr)
        mensaje
        
    })
    
    output$Dispersion2 <- renderPlot({
        
        ggplot(NULL,aes(x=as.double(Trans2()[,1]),y=as.double(variabletrans()[,1])))+
            geom_point(shape=18,color="blue",size=3)+
            geom_smooth(method = lm,linetype="dashed",color="black",fill="peachpuff")+
            labs(title = paste("\n",variabletrans()[1,2], "vs", Trans2()[1,2], "\n"),x=paste(Trans2()[1,2]),y=paste(variabletrans()[1,2],"\n"))+
            theme(plot.title = element_text(color="navy", size=18, face="bold.italic",hjust=0.5),
                  axis.title.x = element_text(color="navy", size=15, face="bold"),
                  axis.text.x = element_text(size=13),
                  axis.title.y = element_text(color="navy", size=15, face="bold"),
                  axis.text.y = element_text(size = 13))
        
        
    })
    
    output$correlacion2 <- renderText({
        
        corr <- cor(as.double(Trans2()[,1]),as.double(variabletrans()[,1]))
        
        mensaje <- paste("Cor = ", corr)
        mensaje
        
    })
    
    output$Dispersion3 <- renderPlot({
        
        ggplot(NULL,aes(x=as.double(Trans3()[,1]),y=as.double(variabletrans()[,1])))+
            geom_point(shape=18,color="blue",size=3)+
            geom_smooth(method = lm,linetype="dashed",color="black",fill="cyan")+
            labs(title = paste("\n",variabletrans()[1,2], "vs", Trans3()[1,2], "\n"),x=paste(Trans3()[1,2]),y=paste(variabletrans()[1,2],"\n"))+
            theme(plot.title = element_text(color="navy", size=18, face="bold.italic",hjust=0.5),
                  axis.title.x = element_text(color="navy", size=15, face="bold"),
                  axis.text.x = element_text(size=13),
                  axis.title.y = element_text(color="navy", size=15, face="bold"),
                  axis.text.y = element_text(size = 13))
        
    })
    
    output$correlacion3 <- renderText({
        
        corr <- cor(as.double(Trans3()[,1]),as.double(variabletrans()[,1]))
        
        mensaje <- paste("Cor = ", corr)
        mensaje
        
    })
    
    output$Dispersion4 <- renderPlot({
        
        ggplot(NULL,aes(x=as.double(Trans4()[,1]),y=as.double(variabletrans()[,1])))+
            geom_point(shape=18,color="blue",size=3)+
            geom_smooth(method = lm,linetype="dashed",color="black",fill="violet")+
            labs(title = paste("\n",variabletrans()[1,2], "vs", Trans4()[1,2], "\n"),x=paste(Trans4()[1,2]),y=paste(variabletrans()[1,2],"\n"))+
            theme(plot.title = element_text(color="navy", size=18, face="bold.italic",hjust=0.5),
                  axis.title.x = element_text(color="navy", size=15, face="bold"),
                  axis.text.x = element_text(size=13),
                  axis.title.y = element_text(color="navy", size=15, face="bold"),
                  axis.text.y = element_text(size = 13))
        
        
    })
    
    output$correlacion4 <- renderText({
        
        corr <- cor(as.double(Trans4()[,1]),as.double(variabletrans()[,1]))
        
        mensaje <- paste("Cor = ", corr)
        mensaje
        
    })
    
    output$Dispersion5 <- renderPlot({
        
        ggplot(NULL,aes(x=as.double(Trans5()[,1]),y=as.double(variabletrans()[,1])))+
            geom_point(shape=18,color="blue",size=3)+
            geom_smooth(method = lm,linetype="dashed",color="black",fill="sandybrown")+
            labs(title = paste("\n",variabletrans()[1,2], "vs", Trans5()[1,2], "\n"),x=paste(Trans5()[1,2]),y=paste(variabletrans()[1,2],"\n"))+
            theme(plot.title = element_text(color="navy", size=18, face="bold.italic",hjust=0.5),
                  axis.title.x = element_text(color="navy", size=15, face="bold"),
                  axis.text.x = element_text(size=13),
                  axis.title.y = element_text(color="navy", size=15, face="bold"),
                  axis.text.y = element_text(size = 13))
        
    })
    
    output$correlacion5 <- renderText({
        
        corr <- cor(as.double(Trans5()[,1]),as.double(variabletrans()[,1]))
        
        mensaje <- paste("Cor = ", corr)
        mensaje
        
    })
    
    output$Dispersion6 <- renderPlot({
        
        ggplot(NULL,aes(x=as.double(Trans6()[,1]),y=as.double(variabletrans()[,1])))+
            geom_point(shape=18,color="blue",size=3)+
            geom_smooth(method = lm,linetype="dashed",color="black",fill="hotpink")+
            labs(title = paste("\n",variabletrans()[1,2], "vs", Trans6()[1,2], "\n"),x=paste(Trans6()[1,2]),y=paste(variabletrans()[1,2],"\n"))+
            theme(plot.title = element_text(color="navy", size=18, face="bold.italic",hjust=0.5),
                  axis.title.x = element_text(color="navy", size=15, face="bold"),
                  axis.text.x = element_text(size=13),
                  axis.title.y = element_text(color="navy", size=15, face="bold"),
                  axis.text.y = element_text(size = 13))
        
    })
    
    output$correlacion6 <- renderText({
        
        corr <- cor(as.double(Trans6()[,1]),as.double(variabletrans()[,1]))
        
        mensaje <- paste("Cor = ", corr)
        mensaje 
        
    })
    
    output$Dispersion7 <- renderPlot({
        
        ggplot(NULL,aes(x=as.double(Trans7()[,1]),y=as.double(variabletrans()[,1])))+
            geom_point(shape=18,color="blue",size=3)+
            geom_smooth(method = lm,linetype="dashed",color="black",fill="yellow")+
            labs(title = paste("\n",variabletrans()[1,2], "vs", Trans7()[1,2], "\n"),x=paste(Trans7()[1,2]),y=paste(variabletrans()[1,2],"\n"))+
            theme(plot.title = element_text(color="navy", size=18, face="bold.italic",hjust=0.5),
                  axis.title.x = element_text(color="navy", size=15, face="bold"),
                  axis.text.x = element_text(size=13),
                  axis.title.y = element_text(color="navy", size=15, face="bold"),
                  axis.text.y = element_text(size = 13)) 
        
    })
    
    output$correlacion7 <- renderText({
        
        corr <- cor(as.double(Trans7()[,1]),as.double(variabletrans()[,1]))
        
        mensaje <- paste("Cor = ", corr)
        mensaje 
        
    })
    
    Modelcuanti <- reactive({
        
        cuantis  <- datos[,c(3:10)]
        xnam <- paste0(colnames(cuantis[as.double(input$variablescuantis)]))
        
        if(length(xnam)==0)
        {
            fmla <- as.formula(paste(deparse(substitute(variabletrans2())), "~ 1"))
            nombrefmla <- paste(variabletrans()[1,2], "~ 1")
        }
        else
        {
            fmla <- as.formula(paste(deparse(substitute(variabletrans2())), "~",paste(xnam,collapse = " + ")))
            nombrefmla <- paste(variabletrans()[1,2], "~",paste(xnam,collapse = " + "))
        }
        Model1 <- lm(fmla)
        Model1[["call"]] <- nombrefmla
        Model1
        
    })
    
    output$Model1 <- renderPrint({
        
        summary(Modelcuanti())
        
    })
    
    output$VIF <- renderPrint({
        
        if(length(Modelcuanti()$coefficients)<=2){
            mensaje="The model must have at least two explanatory variables to execute this function"
            mensaje
        }
        else
        {
            vif(Modelcuanti())
        }
        
    })
    
    output$Alarm <- renderText({
        
        if(length(Modelcuanti()$coefficients)<=2){
            alerta = "The model does not has enough independent variables"
            alerta
        }
        else
        {
            listadevifs <- vif(Modelcuanti())
            
            mensaje="There are no multicollinearity problems"
            nombres <- vector(mode = "numeric",length = 7)
            
            for(i in 1:length(listadevifs)){
                
                if(listadevifs[[i]]>5){
                    
                    mensaje="There are multicollinearity problems in the following variables:"
                    nombres[i] = i
                }
            }
            
            variablesconproblemas <- paste(names(listadevifs[nombres]),collapse = ", ")
            
            
            if(nombres[1]==0 & nombres[2]==0 & nombres[3]==0 & nombres[4]==0 & nombres[5]==0 & nombres[6]==0 & nombres[7]==0){
                mensaje
            }
            else
            {
                paste(mensaje,variablesconproblemas,". You should keep only one of these.") 
            }
        }
        
    })
    
    output$Determination <- renderText({
        
        valoresp <- summary(Modelcuanti())$coefficients[,4]
        
        mensaje1 = "All the model parameters are significant for a confidence level of 95%"
        nombresbetas <- vector(mode = "numeric",length = 6)
        
        for(i in 1:length(valoresp)){
            
            if(valoresp[[i]]>0.05){
                mensaje1="There are parameters which are not significant with a confidence level of 95%, that is, there are relationships between variables that are not important, these parameters correspond to:"
                nombresbetas[i]=i
            }
        }
        
        betasnosignificativos <- paste(names(valoresp[nombresbetas]),collapse = ", ")
        paste(mensaje1, betasnosignificativos)
        
    })
    
    output$AdjustedDetermination <- renderText({
        
        Rajustado <- summary(Modelcuanti())$adj.r.squared
        
        paste("Whit the current model you got an adjusted R squared of",Rajustado, ". But remember, if there are multicollinearity problems this is not a reliable result, because the estimation can be quite imprecise, and the variances and IC will be too broad")
        
    })
    
    Modelofinalfinal <- reactive({
        
        todasvariables <- cbind(datos[,c(3:9)],as.double(Trans1()[,1]),as.double(Trans2()[,1]),as.double(Trans3()[,1]),as.double(Trans4()[,1]),as.double(Trans5()[,1]),as.double(Trans6()[,1]),as.double(Trans7()[,1]))
        todasvariables <- todasvariables[, c(matrix(1:ncol(todasvariables), nrow = 2, byrow = T))]
        
        nombrestodas <- c("Projected population","Transformed Projected population","Thefts","Transformed Thefts","Traffic accidents","Transformed Traffic accidents","Homicides","Transformed Homicides","School deserters","Transformed School deserters","Sports venues", "Transformed Sports venues","Extortions", "Transformed Extortions")
        
        variables <- c(input$variablesincluidas,input$incluirtrans)
        sort(variables,decreasing = F)
        
        xnames <- paste0(colnames(todasvariables[as.double(variables)]))
        
        if(length(xnames)==0)
        {
            fmla2 <- as.formula(paste(deparse(substitute(variabletrans2())), "~ 1"))
            nombrefmla2 <- paste(variabletrans()[1,2], "~ 1")
        }
        else
        {
            fmla2 <- as.formula(paste(deparse(substitute(variabletrans2())), "~",paste(xnames,collapse = " + ")))
            nombrefmla2 <- paste(variabletrans()[1,2], "~",paste(nombrestodas[as.double(variables)],collapse = " + "))
        }
        Model2 <- lm(fmla2)
        Model2[["call"]] <- nombrefmla2
        names(Model2$coefficients) <- c("Intercept",nombrestodas[as.double(variables)])
        Model2
        
    })
    
    output$finalmodel <- renderPrint({
        
        summary(Modelofinalfinal())
        
    })
    
    output$Significancy2 <- renderUI({
        
        valoresp2 <- summary(Modelofinalfinal())$coefficients[,4]
        
        mensaje2 = "All the model parameters are significant for a confidence level of 95%"
        nombresbetas2 <- vector(mode = "numeric",length = 15)
        
        for(i in 1:length(valoresp2)){
            
            if(valoresp2[[i]]>0.05){
                mensaje2="There are parameters which are not significant with a confidence level of 95%, that is, there are relationships between variables that are not important, these parameters correspond to:"
                nombresbetas2[i]=i
            }
        }
        
        betasnosignificativos2 <- paste(names(valoresp2[nombresbetas2]),collapse = ", ")
        p(paste(mensaje2, betasnosignificativos2),style="padding:25px;background-color:lavender;border-left:8px solid blue;border-top: 1px solid black;border-right:1px solid black;border-bottom: 1px solid black;color:black;text-align:center")
        
        
    })
    
    output$Anothermessage <- renderUI({
        
        p("Those variables whose betas are not significant should be eliminated from the model, try to get them out one by one prioritizing those whose betas have a higher p-value (Pr(>|t|))",style="padding:25px;background-color:papayawhip;border-left:8px solid coral;border-top: 1px solid black;border-right:1px solid black;border-bottom: 1px solid black;color:black;text-align:center")
        
        
    })
    
    output$FinalAlarma <- renderUI({
        
        multicoli <- summary(Modelofinalfinal())$coefficients[,4]
        contador = 0
        
        for(i in 1:length(names(multicoli))){
            
            if(names(multicoli[i])=="Projected population"| names(multicoli[i])=="School deserters")
            {
                contador = contador + 1
            }
        }
        
        if(contador >= 2)
        {
            mensaje = "There are multicollinearity problems, remember the analysis of the previous section"
        }
        else
        {
            mensaje = "It seems that there are no problems of multicollinearity, you are doing great so far"
        }
        
        p(mensaje, style="padding:25px;background-color:papayawhip;border-left:8px solid coral;border-top: 1px solid black;border-right:1px solid black;border-bottom: 1px solid black;color:black;text-align:center" )
        
    })
    
    selecciondevariables <- reactive({
        
        todasvariables2 <- cbind(datos[,c(3:9)],as.double(Trans1()[,1]),as.double(Trans2()[,1]),as.double(Trans3()[,1]),as.double(Trans4()[,1]),as.double(Trans5()[,1]),as.double(Trans6()[,1]),as.double(Trans7()[,1]))
        todasvariables2 <- todasvariables2[, c(matrix(1:ncol(todasvariables2), nrow = 2, byrow = T))]
        
        nombrestodas2 <- c("Projected population","Transformed Projected population","Thefts","Transformed Thefts","Traffic accidents","Transformed Traffic accidents","Homicides","Transformed Homicides","School deserters","Transformed School deserters","Sports venues", "Transformed Sports venues","Extortions", "Transformed Extortions")
        
        variables <- c(input$variablesincluidas,input$incluirtrans)
        sort(variables,decreasing = F)
        
        xnames2 <- paste0(colnames(todasvariables2[as.double(variables)]))
        
        if(length(xnames2)==0)
        {
            fmla3 <- as.formula(paste(deparse(substitute(variabletrans2())), "~ 1"))
        }
        else
        {
            fmla3 <- as.formula(paste(deparse(substitute(variabletrans2())), "~",paste(xnames2,collapse = " + ")))
            
        }
        
        Model3 <- lm(fmla3)
        Model3
        
        
        Modelbase <- stepwise(Model3,direction = input$direccion,criterion = input$criterio)
        Modelbase
        
    })
    
    Modelazofinal <- reactive({
        
        formulita <- selecciondevariables()$terms
        
        if(any(grepl("Trans1", names(selecciondevariables()$coefficients))==TRUE)){
            if(any(grepl("ProjectedPopulation", names(selecciondevariables()$coefficients))==TRUE)){
                
            }else{
                formulita <- update.formula(formulita,~ . + ProjectedPopulation)
            }
        }
        
        if(any(grepl("Trans2", names(selecciondevariables()$coefficients))==TRUE)){
            if(any(grepl("Thefts", names(selecciondevariables()$coefficients))==TRUE)){
                
            }else{
                formulita <- update.formula(formulita,~ . + Thefts)
            }
        }
        
        if(any(grepl("Trans3", names(selecciondevariables()$coefficients))==TRUE)){
            if(any(grepl("TrafAccid", names(selecciondevariables()$coefficients))==TRUE)){
                
            }else{
                formulita <- update.formula(formulita,~ . + TrafAccid)
            }
        }
        
        if(any(grepl("Trans4", names(selecciondevariables()$coefficients))==TRUE)){
            if(any(grepl("Homicides", names(selecciondevariables()$coefficients))==TRUE)){
                
            }else{
                formulita <- update.formula(formulita,~ . + Homicides)
            }
        }
        
        if(any(grepl("Trans5", names(selecciondevariables()$coefficients))==TRUE)){
            if(any(grepl("SchoolDes", names(selecciondevariables()$coefficients))==TRUE)){
                
            }else{
                formulita <- update.formula(formulita,~ . + SchoolDes)
            }
        }
        
        if(any(grepl("Trans6", names(selecciondevariables()$coefficients))==TRUE)){
            if(any(grepl("SportsScenari", names(selecciondevariables()$coefficients))==TRUE)){
                
            }else{
                formulita <- update.formula(formulita,~ . + SportsScenari)
            }
        }
        
        if(any(grepl("Trans7", names(selecciondevariables()$coefficients))==TRUE)){
            if(any(grepl("Extortions", names(selecciondevariables()$coefficients))==TRUE)){
                
            }else{
                formulita <- update.formula(formulita,~ . + Extortions)
            }
        }
        
        Model4 <- lm(formulita)
        Model4[["call"]] <- selecciondevariables()$call
        Model4[["call"]] <- gsub("[(]","",Model4[["call"]])
        Model4[["call"]] <- gsub("[)]","",Model4[["call"]])
        Model4[["call"]] <- gsub("variabletrans2",variabletrans()[1,2],Model4[["call"]])
        Model4[["call"]] <- paste0(Model4[["call"]][1],"(formula = ",Model4[["call"]][2],")")
        
        names(Model4$coefficients) <- gsub("[.]","",names(Model4$coefficients))
        names(Model4$coefficients) <- gsub("[,]","",names(Model4$coefficients))
        names(Model4$coefficients) <- gsub("[(]","",names(Model4$coefficients))
        names(Model4$coefficients) <- gsub("[)]","",names(Model4$coefficients))
        names(Model4$coefficients) <- gsub("[[]","",names(Model4$coefficients))
        names(Model4$coefficients) <- gsub("[]]","",names(Model4$coefficients))
        names(Model4$coefficients) <- gsub("asdoubleTrans1 1","Transformed Projected Population",names(Model4$coefficients))
        names(Model4$coefficients) <- gsub("asdoubleTrans2 1","Transformed Thefts",names(Model4$coefficients))
        names(Model4$coefficients) <- gsub("asdoubleTrans3 1","Transformed TrafAccid",names(Model4$coefficients))
        names(Model4$coefficients) <- gsub("asdoubleTrans4 1","Transformed Homicides",names(Model4$coefficients))
        names(Model4$coefficients) <- gsub("asdoubleTrans5 1","Transformed SchoolDes",names(Model4$coefficients))
        names(Model4$coefficients) <- gsub("asdoubleTrans6 1","Transformed SportsScenari",names(Model4$coefficients))
        names(Model4$coefficients) <- gsub("asdoubleTrans3 1","Transformed Extortions",names(Model4$coefficients))
        
        Model4
        
        
    })
    
    output$ModeloBack <- renderPrint({
        
        invisible(capture.output(modelo <- summary(Modelazofinal())))
        modelo
        
    })
    
    output$Determinacionfinal <- renderUI({
        
        coeficiente = summary(Modelazofinal())$adj.r.squared
        p(paste("With the final model you built, you get an adjusted square R of:",coeficiente),style="padding:25px;background-color:papayawhip;border-left:8px solid coral;border-top: 1px solid black;border-right:1px solid black;border-bottom: 1px solid black;color:black;text-align:center" )
        
    })
    
    output$Histograma2 <- renderPlot({
        
        ggplot(NULL,aes(Modelazofinal()$residuals))+geom_histogram(bins=nclass.Sturges(Modelazofinal()$residuals),color="white",
                                                                   fill="seagreen1",aes(y=..density..),lwd=0.8)+geom_density(color="seagreen4",
                                                                                                                             alpha=0.3,fill="seagreen4",lty=1)+
            labs(title = "Residuals \n histogram",x="Residuals",y="Density")+
            theme(plot.title = element_text(color="navy", size=15, face="bold.italic",hjust=0.5),
                  axis.title.x = element_text(color="navy", size=13, face="bold"),
                  axis.title.y = element_text(color="navy", size=13, face="bold"))
        
        
    })
    
    output$Boxplot2 <- renderPlot({
        
        ggplot(NULL,aes(x=0,y=Modelazofinal()$residuals))+geom_boxplot(color="black",fill="skyblue",alpha=0.5)+ stat_summary(fun.y=mean, colour="darkred", geom="point",shape=18, size=3)+
            labs(title = " Residuals \n boxplot",x="",y="Residuals")+
            theme(plot.title = element_text(color="navy", size=15, face="bold.italic",hjust=0.5),
                  axis.title.x = element_text(),
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.title.y = element_text(color="navy", size=13, face="bold"))
        
        
    })
    
    output$qqPlot2 <- renderPlot({
        
        par(font.main=4,font.lab=2,col.main="navy",col.lab="navy",cex.lab=1.1)
        qqPlot(Modelazofinal()$residuals,grid=F,xlab="",ylab="")
        u <-par("usr")
        rect(u[1], u[3], u[2], u[4], col="#EBE9E9", border=TRUE)
        grid(NULL,NULL,col="white",lty=1)
        par(new=TRUE)
        qqPlot(Modelazofinal()$residuals,col="coral",pch=16,id=T,lwd=1.9,col.lines = "black",grid = F, main = "Residuals \n Q-Q plot",xlab="Normal quantiles",ylab="Residuals")
        box(col="white")
        
        
    })
    
    testanalitico2 <- reactive({
        
        if(input$PruebaAnalitica2 == 1){
            
            prueba <- shapiro.test(Modelazofinal()$residuals)
            
        }
        else
        {
            
            if(input$PruebaAnalitica2 ==2){
                
                prueba <- ad.test(Modelazofinal()$residuals)
                
            }
            else
            {
                
                if(input$PruebaAnalitica2 == 3){
                    
                    prueba <- cvm.test(Modelazofinal()$residuals)
                    
                }
                else
                {
                    
                    if(input$PruebaAnalitica2 == 4){
                        
                        prueba <- lillie.test(Modelazofinal()$residuals)
                        
                    }
                    else
                    {
                        
                        prueba <- jarque.bera.test(Modelazofinal()$residuals)
                        
                    }
                    
                }
                
            }
            
            
        }
        
        prueba$data.name <- "Model Residuals"
        prueba
        
    })
    
    output$Prueba2 <- renderPrint({
        
        testanalitico2()
        
    })
    
    output$Conclusion12 <- renderText({
        
        if(testanalitico2()$p.value < 0.05){mensaje="We reject the null hypothesis, the residuals are not normally distributed"}else{mensaje="We keep the null hypothesis, the residuals are normally distributed"}
        mensaje
        
    })
    
    output$ReadMore2 <- renderUI({
        
        
        if(input$PruebaAnalitica2 == 1){
            
            More <- p("Read more about this test here → ", a(href="https://en.wikipedia.org/wiki/Shapiro%E2%80%93Wilk_test", icon("wikipedia-w"),target="_blank"),style="color:black")
            
        }
        else
        {
            
            if(input$PruebaAnalitica2 == 2){
                
                More <- p("Read more about this test here → ", a(href="https://en.wikipedia.org/wiki/Anderson%E2%80%93Darling_test", icon("wikipedia-w"),target="_blank"),style="color:black")
                
            }
            else
            {
                
                if(input$PruebaAnalitica2 == 3){
                    
                    More <- p("Read more about this test here → ", a(href="https://en.wikipedia.org/wiki/Cram%C3%A9r%E2%80%93von_Mises_criterion", icon("wikipedia-w"),target="_blank"),style="color:black")
                    
                }
                else
                {
                    
                    if(input$PruebaAnalitica2 ==4){
                        
                        More <- p("Read more about this test here → ", a(href="https://en.wikipedia.org/wiki/Kolmogorov%E2%80%93Smirnov_test", icon("wikipedia-w"),target="_blank"),style="color:black")
                        
                        
                    }
                    else
                    {
                        
                        More <- p("Read more about this test here → ", a(href="https://en.wikipedia.org/wiki/Jarque%E2%80%93Bera_test", icon("wikipedia-w"),target="_blank"),style="color:black")
                        
                        
                    }
                    
                }
                
            }
            
        }
        
        More
        
    })
    
    output$FittedVsResiduals <- renderPlot({
        
        ggplot(NULL,aes(x=Modelazofinal()$fitted.values,y=Modelazofinal()$residuals))+
            geom_point(shape=18,color="blue",size=3)+
            geom_smooth(method = lm,linetype="dashed",color="black",fill="violet")+
            labs(title = "Model fitted values vs residuals",x="Model fitted values",y="Model residuals")+
            theme(plot.title = element_text(color="navy", size=18, face="bold.italic",hjust=0.5),
                  axis.title.x = element_text(color="navy", size=15, face="bold"),
                  axis.text.x = element_text(size=13),
                  axis.title.y = element_text(color="navy", size=15, face="bold"),
                  axis.text.y = element_text(size = 13))
        
    })
    
    formulit <- reactive({
        
        formulita <- selecciondevariables()$terms
        
        if(any(grepl("Trans1", names(selecciondevariables()$coefficients))==TRUE)){
            if(any(grepl("ProjectedPopulation", names(selecciondevariables()$coefficients))==TRUE)){
                
            }else{
                formulita <- update.formula(formulita,~ . + ProjectedPopulation)
            }
        }
        
        if(any(grepl("Trans2", names(selecciondevariables()$coefficients))==TRUE)){
            if(any(grepl("Thefts", names(selecciondevariables()$coefficients))==TRUE)){
                
            }else{
                formulita <- update.formula(formulita,~ . + Thefts)
            }
        }
        
        if(any(grepl("Trans3", names(selecciondevariables()$coefficients))==TRUE)){
            if(any(grepl("TrafAccid", names(selecciondevariables()$coefficients))==TRUE)){
                
            }else{
                formulita <- update.formula(formulita,~ . + TrafAccid)
            }
        }
        
        if(any(grepl("Trans4", names(selecciondevariables()$coefficients))==TRUE)){
            if(any(grepl("Homicides", names(selecciondevariables()$coefficients))==TRUE)){
                
            }else{
                formulita <- update.formula(formulita,~ . + Homicides)
            }
        }
        
        if(any(grepl("Trans5", names(selecciondevariables()$coefficients))==TRUE)){
            if(any(grepl("SchoolDes", names(selecciondevariables()$coefficients))==TRUE)){
                
            }else{
                formulita <- update.formula(formulita,~ . + SchoolDes)
            }
        }
        
        if(any(grepl("Trans6", names(selecciondevariables()$coefficients))==TRUE)){
            if(any(grepl("SportsScenari", names(selecciondevariables()$coefficients))==TRUE)){
                
            }else{
                formulita <- update.formula(formulita,~ . + SportsScenari)
            }
        }
        
        if(any(grepl("Trans7", names(selecciondevariables()$coefficients))==TRUE)){
            if(any(grepl("Extortions", names(selecciondevariables()$coefficients))==TRUE)){
                
            }else{
                formulita <- update.formula(formulita,~ . + Extortions)
            }
        }
        
        formulita
        
    })
    
    testanalitico3 <- reactive({
        
        
        if(input$PruebaAnalitica3 == 1){
            
            modelo <- Modelazofinal()
            prueba <- bptest(modelo)
            prueba$data.name <- "Model Residuals"
            
        }
        else
        {
            modelo <- lm(formulit())    
            prueba <- ncvTest(modelo)
            
        }
        
        prueba
        
    })
    
    output$Prueba3 <- renderPrint({
        
        testanalitico3()
        
    })
    
    output$Conclusion13 <- renderText({
        
        if(input$PruebaAnalitica3 == 1){
            
            if(testanalitico3()$p.value < 0.05){mensaje="We reject the null hypothesis, the residuals are not homoscedastic"}else{mensaje="We keep the null hypothesis, the residuals are homoscedastic"}
            mensaje}
        else
        {
            if(testanalitico3()$p < 0.05){mensaje="We reject the null hypothesis, the residuals are not homoscedastic"}else{mensaje="We keep the null hypothesis, the residuals are homoscedastic"}
            mensaje}
        
        
    })
    
    output$ReadMore3 <- renderUI({
        
        More <- p("Read more about this test here → ", a(href="https://en.wikipedia.org/wiki/Breusch%E2%80%93Pagan_test", icon("wikipedia-w"),target="_blank"),style="color:black")
        More
        
    })
    
    output$ACF <- renderPlot({
        
        
        par(font.main=4,font.lab=2,col.main="navy",col.lab="navy",cex.lab=1.1)
        acf(Modelazofinal()$residuals,ylim=c(-1,1),main="")
        u <-par("usr")
        rect(u[1], u[3], u[2], u[4], col="#EBE9E9", border=TRUE)
        grid(NULL,NULL,col="white",lty=1)
        par(new=TRUE)
        acf(Modelazofinal()$residuals,ylim=c(-1,1),main="ACF de los residuales del modelo")
        box(col="white")
        
    })
    
    output$PACF <- renderPlot({
        
        par(font.main=4,font.lab=2,col.main="navy",col.lab="navy",cex.lab=1.1)
        pacf(Modelazofinal()$residuals,ylim=c(-1,1),main="")
        u <-par("usr")
        rect(u[1], u[3], u[2], u[4], col="#EBE9E9", border=TRUE)
        grid(NULL,NULL,col="white",lty=1)
        par(new=TRUE)
        pacf(Modelazofinal()$residuals,ylim=c(-1,1),main="ACF de los residuales del modelo")
        box(col="white")
        
    })
    
    output$ResVsIndex <- renderPlot({
        
        ggplot(NULL,aes(x=seq(1, 118, 1),y=Modelazofinal()$residuals))+
            geom_point(shape=18,color="blue",size=3)+
            geom_smooth(method = lm,linetype="dashed",color="black",fill="violet")+
            labs(title = "Residuals vs Index",x="Index",y="Model residuals")+
            theme(plot.title = element_text(color="navy", size=18, face="bold.italic",hjust=0.5),
                  axis.title.x = element_text(color="navy", size=15, face="bold"),
                  axis.text.x = element_text(size=13),
                  axis.title.y = element_text(color="navy", size=15, face="bold"),
                  axis.text.y = element_text(size = 13))
        
    })
    
    testanalitico4 <- reactive({
        
        if(input$PruebaAnalitica4 == 1){
            
            prueba <- dwtest(Modelazofinal(),alternative = "two.sided")
            
        }
        else
        {
            
            prueba <- bgtest(Modelazofinal())
            
        } 
        
        prueba$data.name <- "Model Residuals"   
        prueba
        
    })
    
    output$Prueba4 <- renderPrint({
        
        testanalitico4()
        
    })
    
    output$Conclusion14 <- renderText({
        
        if(testanalitico4()$p.value < 0.05){mensaje="We reject the null hypothesis, the residuals are auto-correlated"}else{mensaje="We keep the null hypothesis, the residuals are not auto-correlated"}
        mensaje
        
    })
    
    output$ReadMore4 <- renderUI({
        
        
        if(input$PruebaAnalitica4 == 1){
            
            More <- p("Read more about this test here → ", a(href="https://en.wikipedia.org/wiki/Durbin%E2%80%93Watson_statistic", icon("wikipedia-w"),target="_blank"),style="color:black")
            
        }
        else
        {
            
            More <- p("Read more about this test here → ", a(href="https://en.wikipedia.org/wiki/Breusch%E2%80%93Godfrey_test", icon("wikipedia-w"),target="_blank"),style="color:black")
        }
        
        More
        
    })
    
    output$Answer1 <- renderUI({
        
        actionID <- input$Question1
        if(!is.null(actionID)){
            
            if(input$Question1 == 1){mensaje = "It is right! this equation represents a straight line, a function of two variables that are related in a perfect or deterministic way."}
            if(input$Question1 == 2){mensaje = "It is not correct! note that in this equation there is a term that represents random factors (e). These random factors make the relationship between X and Y not perfect or deterministic."}
            if(input$Question1 == 3){mensaje = "It is not correct! this equation represents an adjustment, note that the equation includes the adjusted response variable and not the original."}
            if(input$Question1 == 4){mensaje = "This relationship includes an adjusted variable and a constant factor, therefore it is not even a relationship between two variables."}
            
            if(input$Question1 == 1){p(mensaje,icon("smile-wink","fa-2x"),style="background-color:#BFF7BB;padding:15px;text-align:justify;border-left: 8px solid green")} 
            else
            {
                p(mensaje,icon("sad-cry","fa-2x"),style="background-color:#FFA8A8;padding:15px;text-align:justify;border-left: 8px solid red")
            }
            
        }
        else
        {""}
        
    })
    
    output$Answer2 <- renderUI({
        
        actionID <- input$Question2
        if(!is.null(actionID)){
            
            if(input$Question2 == 1){mensaje = "It is not correct! This is an important factor to make inferences in these models but, it is the only one?"}
            if(input$Question2 == 2){mensaje = "It is not correct! This is an important factor to make inferences in these models but, it is the only one?"}
            if(input$Question2 == 3){mensaje = "It is correct! We need the variability of the response variable to be explained by the independent variables (option a) and also that the statistical assumptions be fulfilled (option b)."}
            if(input$Question2 == 4){mensaje = "It is not correct, just try again!"}
            
            if(input$Question2 == 3){p(mensaje,icon("smile-wink","fa-2x"),style="background-color:#BFF7BB;padding:15px;text-align:justify;border-left: 8px solid green")} 
            else
            {
                p(mensaje,icon("sad-cry","fa-2x"),style="background-color:#FFA8A8;padding:15px;text-align:justify;border-left: 8px solid red")
            }
            
        }
        else
        {""}
        
    })
    
    output$Answer3 <- renderUI({
        
        actionID <- input$Question3
        if(!is.null(actionID)){
            
            if(input$Question3 == 1){mensaje = "It is right! For a null hypothesis where the residuals of the model are not autocorrelated, the Durbin Watson independence test shows that there is at least a lag with a significant autocorrelation coefficient. "}
            if(input$Question3 == 2){mensaje = "It is not correct! We cannot know this using the Durbin Watson test which purpose is proving independence of residuals"}
            if(input$Question3 == 3){mensaje = "It is not correct! We cannot know if every autocorrelation coefficient is significant, to prove this we should do graphical tests like ACF"}
            if(input$Question3 == 4){mensaje = "It is not correct! In fact we can conclude about a null hypothesis where the residuals of the model are not autocorrelated and we can do it using the p - value"}
            
            if(input$Question3 == 1){p(mensaje,icon("smile-wink","fa-2x"),style="background-color:#BFF7BB;padding:15px;text-align:justify;border-left: 8px solid green")} 
            else
            {
                p(mensaje,icon("sad-cry","fa-2x"),style="background-color:#FFA8A8;padding:15px;text-align:justify;border-left: 8px solid red")
            }
            
        }
        else
        {""}
        
    })
    
    output$Answer4 <- renderUI({
        
        actionID <- input$Question4
        if(!is.null(actionID)){
            
            if(input$Question4 == 1){mensaje = "It is not right! Remember that the root of the R squared is the coefficient of autocorrelation, and remember that the r squared is calculated in the following way: "}
            if(input$Question4 == 2){mensaje = "It is correct! You made the right calculations. Was it just luck? Check then the concepts of R squared and coefficient of autocorrelation in the section Glossary"}
            if(input$Question4 == 3){mensaje = "It is not correct! Remember that the root of the R squared is the coefficient of autocorrelation, and remember that the r squared is calculated in the following way: "}
            if(input$Question4 == 4){mensaje = "It is not correct! Remember that the root of the R squared is the coefficient of autocorrelation, and remember that the r squared is calculated in the following way: "}
            
            if(input$Question4 == 2){p(mensaje,icon("smile-wink","fa-2x"),style="background-color:#BFF7BB;padding:15px;text-align:justify;border-left: 8px solid green")} 
            else
            {
                p(mensaje,withMathJax('\\( R^2 = \\frac{SSR}{SSR+SSE} \\)'),icon("sad-cry","fa-2x"),style="background-color:#FFA8A8;padding:15px;text-align:justify;border-left: 8px solid red")
            }
            
        }
        else
        {""}
        
    })
    
    output$Answer5 <- renderUI({
        
        actionID <- input$Question5
        if(!is.null(actionID)){
            
            if(input$Question5 == 1){mensaje = "It is not right! This mathematical expression tells us that the sum of the real values and the adjusted values are equal, therefore at a general level the adjustment could be considered perfect, but we must understand the point-to-point adjustment in other way"}
            if(input$Question5 == 2){mensaje = "It is not correct! If the point-to-point adjustment was perfect there were no terms associated with randomness and uncontrollable factors"}
            if(input$Question5 == 3){mensaje = "It is not correct! If the point-to-point adjustment was perfect there were no terms associated with randomness and uncontrollable factors"}
            if(input$Question5 == 4){mensaje = "It is correct! This mathematical expression tells us that the sum of the real values and the adjusted values are equal, therefore at a general level the adjustment could be considered perfect, even when in each point there is an error percentage"}
            
            if(input$Question5 == 4){p(mensaje,icon("smile-wink","fa-2x"),style="background-color:#BFF7BB;padding:15px;text-align:justify;border-left: 8px solid green")} 
            else
            {
                p(mensaje,icon("sad-cry","fa-2x"),style="background-color:#FFA8A8;padding:15px;text-align:justify;border-left: 8px solid red")
            }
            
        }
        else
        {""}
        
    })
    
    output$Answer6 <- renderUI({
        
        actionID <- input$Question6
        if(!is.null(actionID)){
            
            if(input$Question6 == 1){mensaje = "It is not correct! Actually the r squared increases as more variables have the model."}
            if(input$Question6 == 2){mensaje = "That is right! That is exactly the answer, but was it just luck? then review the glossary section to better understand the concepts of R squared and adjusted R squared."}
            if(input$Question6 == 3){mensaje = "It is not correct! This does not make much sense, review the glossary section to better understand the concepts of R squared and adjusted R squared."}
            if(input$Question6 == 4){mensaje = "It is not correct! This does not make much sense, review the glossary section to better understand the concepts of R squared and adjusted R squared."}
            
            if(input$Question6 == 2){p(mensaje,icon("smile-wink","fa-2x"),style="background-color:#BFF7BB;padding:15px;text-align:justify;border-left: 8px solid green")} 
            else
            {
                p(mensaje,icon("sad-cry","fa-2x"),style="background-color:#FFA8A8;padding:15px;text-align:justify;border-left: 8px solid red")
            }
            
        }
        else
        {""}
        
    })
    
    output$Answer7 <- renderUI({
        
        actionID <- input$Question7
        if(!is.null(actionID)){
            
            if(input$Question7 == 1){mensaje = "It is not correct! Actually the points of the diagram are not so scattered and we can see a not so weak linear relationship."}
            if(input$Question7 == 2){mensaje = "That is not right! In the diagram we can see how the points follow the shape of a straight line."}
            if(input$Question7 == 3){mensaje = "It is not correct! The relationship seems to be strong so it would be worthwhile to perform a regression model between these variables."}
            if(input$Question7 == 4){mensaje = "It is right! The relationship seems to be strong so it would be worthwhile to perform a regression model between these variables."}
            
            if(input$Question7 == 4){p(mensaje,icon("smile-wink","fa-2x"),style="background-color:#BFF7BB;padding:15px;text-align:justify;border-left: 8px solid green")} 
            else
            {
                p(mensaje,icon("sad-cry","fa-2x"),style="background-color:#FFA8A8;padding:15px;text-align:justify;border-left: 8px solid red")
            }
            
        }
        else
        {""}
        
    })
    
    output$Answer8 <- renderUI({
        
        actionID <- input$Question8
        if(!is.null(actionID)){
            
            if(input$Question8 == 1){mensaje = "It is not correct! Graphing the residuals of the model and its adjusted values does not tell me anything about the assumption of normality."}
            if(input$Question8 == 2){mensaje = "That is right! We can see how the points are not uniformly distributed so the variance is not constant and a solution could be to transform the response variable."}
            if(input$Question8 == 3){mensaje = "It is not correct! We can see how the points are not uniformly distributed so the variance is not constant."}
            if(input$Question8 == 4){mensaje = "It is not correct! Graphing the residuals of the model and its adjusted values does not tell me anything about the assumption of normality."}
            
            if(input$Question8 == 2){p(mensaje,icon("smile-wink","fa-2x"),style="background-color:#BFF7BB;padding:15px;text-align:justify;border-left: 8px solid green")} 
            else
            {
                p(mensaje,icon("sad-cry","fa-2x"),style="background-color:#FFA8A8;padding:15px;text-align:justify;border-left: 8px solid red")
            }
            
        }
        else
        {""}
        
    })
    
    output$Answer9 <- renderUI({
        
        actionID <- input$Question9
        if(!is.null(actionID)){
            
            if(input$Question9 == 1){mensaje = "It is not correct! We can see how there are many points outside the confidence bands so there is no normality."}
            if(input$Question9 == 2){mensaje = "That is not correct! This graph is not useful to conclude about the assumption of constant variance."}
            if(input$Question9 == 3){mensaje = "It is not correct! This graph is not useful to conclude about the assumption of constant variance."}
            if(input$Question9 == 4){mensaje = "It is correct! We can see how there are many points outside the confidence bands so there is no normality."}
            
            if(input$Question9 == 4){p(mensaje,icon("smile-wink","fa-2x"),style="background-color:#BFF7BB;padding:15px;text-align:justify;border-left: 8px solid green")} 
            else
            {
                p(mensaje,icon("sad-cry","fa-2x"),style="background-color:#FFA8A8;padding:15px;text-align:justify;border-left: 8px solid red")
            }
            
        }
        else
        {""}
        
    })
    
    output$Answer10 <- renderUI({
        
        actionID <- input$Question10
        if(!is.null(actionID)){
            
            if(input$Question10 == 1){mensaje = "It is not correct! Indeed, the relationship is positive, but is it linear? It does not look like a straight line, does it?"}
            if(input$Question10 == 2){mensaje = "That is not correct! Indeed, the relationship is non-linear, but is it negative? It does not seem to decrease, does it?"}
            if(input$Question10 == 3){mensaje = "It is correct! We can see the positive relationship with a functional form that is not linear."}
            if(input$Question10 == 4){mensaje = "It is not correct! Can not you see any functional form?"}
            
            if(input$Question10 == 3){p(mensaje,icon("smile-wink","fa-2x"),style="background-color:#BFF7BB;padding:15px;text-align:justify;border-left: 8px solid green")} 
            else
            {
                p(mensaje,icon("sad-cry","fa-2x"),style="background-color:#FFA8A8;padding:15px;text-align:justify;border-left: 8px solid red")
            }
            
        }
        else
        {""}
        
    })
    
    output$Answer11 <- renderUI({
        
        actionID <- input$Question11
        if(!is.null(actionID)){
            
            if(input$Question11 == 1){mensaje = "It is not correct! Can you see any functional form?"}
            if(input$Question11 == 2){mensaje = "That is not correct! Can you see any functional form?"}
            if(input$Question11 == 3){mensaje = "It is not correct! Can you see any functional form?"}
            if(input$Question11 == 4){mensaje = "It is correct! There is no association between the variables, neither positive, nor negative, nor with any functional form."}
            
            if(input$Question11 == 4){p(mensaje,icon("smile-wink","fa-2x"),style="background-color:#BFF7BB;padding:15px;text-align:justify;border-left: 8px solid green")} 
            else
            {
                p(mensaje,icon("sad-cry","fa-2x"),style="background-color:#FFA8A8;padding:15px;text-align:justify;border-left: 8px solid red")
            }
            
        }
        else
        {""}
        
    })
    
    output$Answer12 <- renderUI({
        
        actionID <- input$Question12
        if(!is.null(actionID)){
            
            if(input$Question12 == 1){mensaje = "It is not correct! 0.67668 is the value of the parameter associated with the traffic accidents variable and should not be interpreted as a percentage in this case."}
            if(input$Question12 == 2){mensaje = "That is not correct! 0.67668 is the value of the parameter associated with the traffic accidents variable and should not be interpreted as a percentage in this case."}
            if(input$Question12 == 3){mensaje = "That is right!"}
            if(input$Question12 == 4){mensaje = "It is not correct! Think again"}
            
            if(input$Question12 == 3){p(mensaje,icon("smile-wink","fa-2x"),style="background-color:#BFF7BB;padding:15px;text-align:justify;border-left: 8px solid green")} 
            else
            {
                p(mensaje,icon("sad-cry","fa-2x"),style="background-color:#FFA8A8;padding:15px;text-align:justify;border-left: 8px solid red")
            }
            
        }
        else
        {""}
        
    })
    
    output$Answer13 <- renderUI({
        
        actionID <- input$Question13
        if(!is.null(actionID)){
            
            if(input$Question13 == 1){mensaje = "It is not correct! Remember that the problem of multicollinearity is greater the larger the value of the VIF."}
            if(input$Question13 == 2){mensaje = "That is not correct! It is true that this variable has the greatest multicollinearity problem but alternatives must be considered before eliminating it."}
            if(input$Question13 == 3){mensaje = "It is correct! It is important to try to make indicators before eliminating variables to group all the information that both variables can provide."}
            if(input$Question13 == 4){mensaje = "It is not correct! This statement does not make much sense, review step 3 to remember this concept."}
            
            if(input$Question13 == 3){p(mensaje,icon("smile-wink","fa-2x"),style="background-color:#BFF7BB;padding:15px;text-align:justify;border-left: 8px solid green")} 
            else
            {
                p(mensaje,icon("sad-cry","fa-2x"),style="background-color:#FFA8A8;padding:15px;text-align:justify;border-left: 8px solid red")
            }
            
        }
        else
        {""}
        
    })
    
    output$Answer14 <- renderUI({
        
        actionID <- input$Question14
        if(!is.null(actionID)){
            
            if(input$Question14 == 1){mensaje = "It is correct! With this small p value we can conclude that the parameter that represents the linear relationship (slope of the regression line) between the variables is significant."}
            if(input$Question14 == 2){mensaje = "That is not correct! With this small p value we can conclude that the parameter that represents the linear relationship (slope of the regression line) between the variables is significant."}
            if(input$Question14 == 3){mensaje = "It is not correct! Are you sure that the parameter is not significant?"}
            if(input$Question14 == 4){mensaje = "It is not correct! First, look at the p value, are you sure that the parameter is not significant?"}
            
            if(input$Question14 == 1){p(mensaje,icon("smile-wink","fa-2x"),style="background-color:#BFF7BB;padding:15px;text-align:justify;border-left: 8px solid green")} 
            else
            {
                p(mensaje,icon("sad-cry","fa-2x"),style="background-color:#FFA8A8;padding:15px;text-align:justify;border-left: 8px solid red")
            }
            
        }
        else
        {""}
        
    })
    
    output$Answer15 <- renderUI({
        
        actionID <- input$Question15
        if(!is.null(actionID)){
            
            if(input$Question15 == 1){mensaje = "It is not correct! We can not conclude on the goodness of fit having only one test statistic (F in this case)."}
            if(input$Question15 == 2){mensaje = "That is not correct! The multiple R squared is not given directly in terms of percentage."}
            if(input$Question15 == 3){mensaje = "It is correct! The adjusted R squared converted to percentage tells us the variability of the response variable that is explained by the independent variables."}
            if(input$Question15 == 4){mensaje = "It is not correct! The presented p-value does not show the statistical validity of the model."}
            
            if(input$Question15 == 3){p(mensaje,icon("smile-wink","fa-2x"),style="background-color:#BFF7BB;padding:15px;text-align:justify;border-left: 8px solid green")} 
            else
            {
                p(mensaje,icon("sad-cry","fa-2x"),style="background-color:#FFA8A8;padding:15px;text-align:justify;border-left: 8px solid red")
            }
            
        }
        else
        {""}
        
    })
    
    output$Answer16 <- renderUI({
        
        actionID <- input$Question16
        if(!is.null(actionID)){
            
            if(input$Question16 == 1){mensaje = "It is not correct! There is no deterministic relationship with the response variable, in addition the adjusted R squared and R squared coefficients do not exceed 80% of variability explained."}
            if(input$Question16 == 2){mensaje = "That is not correct! Before affirming this, another type of relationship should be considered, not just linear."}
            if(input$Question16 == 3){mensaje = "It is right! it seems that there are variables whose linear relationship with the response variable is not significant but may have other types of relationships."}
            if(input$Question16 == 4){mensaje = "It is not correct! Sports venues and extortions are important variables that should not be excluded"}
            
            if(input$Question16 == 3){p(mensaje,icon("smile-wink","fa-2x"),style="background-color:#BFF7BB;padding:15px;text-align:justify;border-left: 8px solid green")} 
            else
            {
                p(mensaje,icon("sad-cry","fa-2x"),style="background-color:#FFA8A8;padding:15px;text-align:justify;border-left: 8px solid red")
            }
            
        }
        else
        {""}
        
    })
    
})