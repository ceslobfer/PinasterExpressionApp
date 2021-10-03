library(shiny)
library(shinyjs)
library(visNetwork)
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
tablaSeqGo2 <- read.csv("www/seq_go_embrios.txt", sep = "\t",col.names = c("seqName","GoID"),stringsAsFactors = F)
tablaSeqDes <- read.csv("www/seq_description.txt", sep = "\t",col.names = c("seq","description"), row.names = 1,stringsAsFactors = F)
RM_data <- read.csv("www/CPM_raizMicrodiseccion_RM_mean.txt", sep = "\t", row.names = 1, stringsAsFactors = F )
tablaSeqDesRoots <- read.csv("www/seq_description_raiz.txt", sep = "\t",col.names = c("seq","description"), row.names = 1,stringsAsFactors = F)
tablaDiffRoots <-read.csv("www/microdiseccion_diffExpression.txt", sep = "\t", col.names = c("seqName","Tissue","Condition","Log Fold Change"), stringsAsFactors = F )
tablaDiffEmbryos <-read.csv("www/embryo_differential.txt", sep = "\t", col.names = c("seqName","Comparative","Type","Log Fold Change"), stringsAsFactors = F )
tablaSeqGoR <- read.csv("www/seq_gos_raiz.txt", sep = "\t",col.names = c("seqName","GoID"),stringsAsFactors = F)
keggTable <- read.csv("www/kegg_total_transcriptome.txt", col.names = c("seqName","KeggID"),sep = "\t", stringsAsFactors = F )
expressionNeedles <- read.csv("www/Needles_CPM_mean.txt", sep = "\t", row.names = 1 )
tablaDiffNeedles <-read.csv("www/needles_differential.txt", sep = "\t", col.names = c("seqName","Condition","Log Fold Change"), stringsAsFactors = F )
module_dataE2 <- read.csv("www/seq_module_file_Embryos.txt", sep = "\t",col.names = c("seqName","module"), stringsAsFactors = F )
module_dataN2 <- read.csv("www/needles_module.txt", sep = "\t",col.names = c("seqName","module"), stringsAsFactors = F )
hubsNeedles <- read.csv("www/hubGenes_needles.txt", sep = "\t",col.names = c("seqName","module"),stringsAsFactors = F)
goPerSeq <- read.csv("www/seq_go_embrios_join.txt", sep = "\t",col.names = c("seqName","go"),row.names = 1,stringsAsFactors = F)

source("Components/Comp_Tables.R")
source("Components/Comp_embryos.R")
source("Components/Comp_needle.R")
source("Components/Comp_barplot.R")
source("Components/Comp_corTable.R")
source("Components/Comp_roots.R")
source("Components/Comp_blast.R")
source("Components/Comp_heatmaps.R")
source("Components/Comp_Network.R")
shinyServer(function(input, output) {
    
    #Network gene Needles
    values <- reactiveValues()
    values$list <- list()
    edgeNodes <- reactiveValues()
    edgeNodes$list <- list()
    observeEvent(input$buttonAddNetworkN, {
        if (input$seqID_N_N %in% row.names(expressionNeedles)) {
            if (input$geneName_N_N != "") {
                values$list[[input$seqID_N_N]] <- list("Name" = input$geneName_N_N,"expression" = expressionNeedles[input$seqID_N_N,])
            }else{
                values$list[[input$seqID_N_N]] <- list("Name" = input$seqID_N_N,"expression" = expressionNeedles[input$seqID_N_N,])
            }
            
            edgeNodes$list <- corGenerationNeedles(values$list,0.05,0)
        }
    })
    
    observeEvent(input$buttonDeleteNetwork, {
        if (input$seqID_N_N %in% names(values$list)) {
            values$list[[input$seqID_N_N]] <- NULL
            edgeNodes$list <- corGenerationNeedles(values$list,0.05,0)
        }
    })
    
    observeEvent(input$fileNetNeedles,{
        inFile <- input$fileNetNeedles
        listSeq <- read.csv(inFile$datapath, header = F,stringsAsFactors = F)
        if (ncol(listSeq) == 1) {
            for (i in c(1:nrow(listSeq))) {
                seqName <- listSeq$V1[i]
                if (seqName %in% row.names(expressionNeedles)) {
                    values$list[[seqName]] <- list("Name" = seqName,"expression" = expressionNeedles[seqName,])
                }
            }
        }
        if (ncol(listSeq) == 2) {
            for (i in c(1:nrow(listSeq))) {
                seqName <- listSeq$V1[i]
                nameGene <- listSeq$V2[i]
                if (seqName %in% row.names(expressionNeedles)) {
                    values$list[[seqName]] <- list("Name" = nameGene,"expression" = expressionNeedles[seqName,])
                }
            }
        }
        print(values$list)
        edgeNodes$list <- corGenerationNeedles(values$list,0.05,0)
    })
    
    observeEvent(input$buttonEmptyNetworkN, {
        for (seqID in names(values$list)) {
            values$list[[seqID]] <- NULL
        }
        edgeNodes$list <- corGenerationNeedles(values$list,0.05,0)
    })
    observeEvent(input$buttonFilterNetworkN,{
        edgeNodes$list <- corGenerationNeedles(values$list,input$CortePvalueNetN,input$CorteCorNetN)
    })
    
    output$networkN <- renderVisNetwork({
        resNet <- edgeNodes$list
        if (length(edgeNodes$list) != 0) {
            nodes <- resNet[[1]]
            edges <- resNet[[2]]
            visNetwork(nodes, edges)  %>% 
                visOptions(highlightNearest = TRUE, height ="650px",nodesIdSelection = TRUE) %>% 
                visNodes(borderWidth = 2,physics = FALSE, font = list(size = 17, face="arial black"),
                         shadow = list(enabled = TRUE, size = 8), shape = "ellipse")%>% 
                visInteraction(navigationButtons = TRUE) %>% 
                visIgraphLayout(physics = FALSE, randomSeed = 12) %>% 
                visEdges(shadow = TRUE,width = 4,length = 400, smooth = FALSE,physics = FALSE)
            
        }
        
    })
#Network gene Embryos
    values <- reactiveValues()
    values$list <- list()
    edgeNodes <- reactiveValues()
    edgeNodes$list <- list()
    observeEvent(input$buttonAddNetwork, {
        if (input$seqID_E_N %in% row.names(datosEmbriones)) {
            if (input$geneName_E_N != "") {
                values$list[[input$seqID_E_N]] <- list("Name" = input$geneName_E_N,"expression" = datosEmbriones[input$seqID_E_N,])
            }else{
                values$list[[input$seqID_E_N]] <- list("Name" = input$seqID_E_N,"expression" = datosEmbriones[input$seqID_E_N,])
            }
        
            edgeNodes$list <- corGenerationEmbryos(values$list,0.05,0)
        }
    })
    
    observeEvent(input$buttonDeleteNetwork, {
        if (input$seqID_E_N %in% names(values$list)) {
            values$list[[input$seqID_E_N]] <- NULL
            edgeNodes$list <- corGenerationEmbryos(values$list,0.05,0)
        }
    })
    
    observeEvent(input$fileNetEmbryos,{
        inFile <- input$fileNetEmbryos
        listSeq <- read.csv(inFile$datapath, header = F,stringsAsFactors = F)
        if (ncol(listSeq) == 1) {
            for (i in c(1:nrow(listSeq))) {
                seqName <- listSeq$V1[i]
                if (seqName %in% row.names(datosEmbriones)) {
                    values$list[[seqName]] <- list("Name" = seqName,"expression" = datosEmbriones[seqName,])
                }
            }
        }
        if (ncol(listSeq) == 2) {
            for (i in c(1:nrow(listSeq))) {
                seqName <- listSeq$V1[i]
                nameGene <- listSeq$V2[i]
                if (seqName %in% row.names(datosEmbriones)) {
                    values$list[[seqName]] <- list("Name" = nameGene,"expression" = datosEmbriones[seqName,])
                }
            }
        }
        print(values$list)
        edgeNodes$list <- corGenerationEmbryos(values$list,0.05,0)
    })
    
    observeEvent(input$buttonEmptyNetwork, {
        for (seqID in names(values$list)) {
            values$list[[seqID]] <- NULL
        }
        edgeNodes$list <- corGenerationEmbryos(values$list,0.05,0)
    })
    observeEvent(input$buttonFilterNetwork,{
        edgeNodes$list <- corGenerationEmbryos(values$list,input$CortePvalueNetE,input$CorteCorNetE)
    })

    output$network <- renderVisNetwork({
        resNet <- edgeNodes$list
        if (length(edgeNodes$list) != 0) {
            nodes <- resNet[[1]]
            edges <- resNet[[2]]
            visNetwork(nodes, edges)  %>% 
                visOptions(highlightNearest = TRUE, height ="650px",nodesIdSelection = TRUE) %>% 
                visNodes(borderWidth = 2,physics = FALSE, font = list(size = 17, face="arial black"),
                         shadow = list(enabled = TRUE, size = 8), shape = "ellipse")%>% 
                visInteraction(navigationButtons = TRUE) %>% 
                visIgraphLayout(physics = FALSE, randomSeed = 12) %>% 
                visEdges(shadow = TRUE,width = 4,length = 400, smooth = FALSE,physics = FALSE)
                
        }
        
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
#General table Needles
    generalTableN <- eventReactive(input$buttonGeneralN, {
        seqList <- c()
        #Differential Condition
        if (!is.null(input$checkGroupN) && input$checkGroupN != "" ) {
            for (condition in input$checkGroupN) {
                condition <- condition
                searchCondition <- tablaDiffNeedles[which(tablaDiffNeedles$Condition== condition),1]
                if (length(seqList) == 0) {
                    seqList <- searchCondition
                }else{
                    seqList <- intersect(seqList, searchCondition)
                }
            }
        }
        
        #Go Condition
        if (!is.null(input$goNTotal)) {
            goSearch <- strsplit(input$goNTotal,split = ",")[[1]]
            for (go in goSearch) {
                seqsGo <- tablaSeqGo2[which(tablaSeqGo2$GoID == go),1]
                if (length(seqList) == 0) {
                    seqList <- seqsGo
                }else{
                    seqList <- intersect(seqList, seqsGo)
                }
            }
        }
        if (input$moduleNTotal != "Module") {
            seqsModule <- module_dataN2[which(module_dataN2$module==input$moduleNTotal),]$seqName
            if (length(seqList) == 0) {
                seqList <- seqsModule
            }else{
                seqList <- intersect(seqList, seqsModule)
            }
        }
        tablaSearch <- data.frame(stringsAsFactors = F)
        if (!is.null(seqList)) {
            
            for (seq in seqList) {
                
                if (seq %in% hubsNeedles$seqName) {
                    tablaSearch <- rbind.data.frame(tablaSearch,c(seq, tablaSeqDes[seq,],goPerSeq[seq,],"HUB"),stringsAsFactors = F)
                }else{
                    tablaSearch <- rbind.data.frame(tablaSearch,c(seq, tablaSeqDes[seq,],goPerSeq[seq,],"NO-H"),stringsAsFactors = F)
                }
                
            }
        }
        if (input$goNTotal == "" && input$moduleNTotal == "Module" && is.null(input$checkGroupN)) {
            seqList <- row.names(expressionNeedles)
            for (seq in seqList) {
                if (seq %in% hubsNeedles$seqName) {
                    tablaSearch <- rbind.data.frame(tablaSearch,c(seq, tablaSeqDes[seq,],goPerSeq[seq,],"HUB"),stringsAsFactors = F)
                }else{
                    tablaSearch <- rbind.data.frame(tablaSearch,c(seq, tablaSeqDes[seq,],goPerSeq[seq,],"NO-H"),stringsAsFactors = F)
                }
                
            }
        }
        
        
        tablaSearch <- tablaSearch[which(tablaSearch[,1] %in% row.names(expressionNeedles)),]
        colnames(tablaSearch) <- c("Seq ID","Seq Description","GO ID", "HUBGENE")
        tablaSearch
        
        
        
    })
    
    output$tableGeneralN <- DT::renderDataTable({
        table <- generalTableN()
        if (length(rownames(table)) > 0) {
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
    heatmapValuesN <- reactiveValues()
    observeEvent(input$fileHeatmapNeedles,{
        inFile <- input$fileHeatmapNeedles
        listSeq <- read.csv(inFile$datapath, header = T,stringsAsFactors = F)
        table <- expressionNeedles[listSeq[,1],]
        heatmapValuesN$table <- table
    })
    #Descarga heatmap Needles
    output$downloadHeatNeedles <- downloadHandler(
        filename = function(){paste("heatMapNeedles",'.png',sep='')},
        content = function(file){
            ggsave(filename = file, plot = heatMapNeedles(heatmapValuesN$table,input$clusterHeatN) )
        }
    )
#General table Embryos
    generalTableE <- eventReactive(input$buttonGeneralE, {
        seqList <- c()
        #Differential Condition
        if (!is.null(input$checkGroupE) && input$checkGroupE != "" ) {
            for (condition in input$checkGroupE) {
                condition <- strsplit(condition,split = ",")[[1]]
                searchCondition <- tablaDiffEmbryos[which(tablaDiffEmbryos$Comparative == condition[1] & tablaDiffEmbryos$Type == condition[2]),1]
                if (length(seqList) == 0) {
                    seqList <- searchCondition
                }else{
                    seqList <- intersect(seqList, searchCondition)
                }
            }
        }
        
        #Go Condition
        if (!is.null(input$goETotal)) {
            goSearch <- strsplit(input$goETotal,split = ",")[[1]]
            for (go in goSearch) {
                seqsGo <- tablaSeqGo2[which(tablaSeqGo2$GoID == go),1]
                if (length(seqList) == 0) {
                    seqList <- seqsGo
                }else{
                    seqList <- intersect(seqList, seqsGo)
                }
            }
        }
        if (input$moduleETotal != "Module") {
            seqsModule <- module_dataE2[which(module_dataE2$module==input$moduleETotal),]$seqName
            if (length(seqList) == 0) {
                seqList <- seqsModule
            }else{
                seqList <- intersect(seqList, seqsModule)
            }
        }
        tablaSearch <- data.frame(stringsAsFactors = F)
        if (!is.null(seqList)) {
            
            for (seq in seqList) {
                
                if (seq %in% hubsEmbrio$seqName) {
                    tablaSearch <- rbind.data.frame(tablaSearch,c(seq, tablaSeqDes[seq,],goPerSeq[seq,],"HUB"),stringsAsFactors = F)
                }else{
                    tablaSearch <- rbind.data.frame(tablaSearch,c(seq, tablaSeqDes[seq,],goPerSeq[seq,],"NO-H"),stringsAsFactors = F)
                }
                
            }
        }
        if (input$goETotal == "" && input$moduleETotal == "Module" && is.null(input$checkGroupE)) {
            seqList <- row.names(datosEmbriones)
                for (seq in seqList) {
                    
                    if (seq %in% hubsEmbrio$seqName) {
                        tablaSearch <- rbind.data.frame(tablaSearch,c(seq, tablaSeqDes[seq,],goPerSeq[seq,],"HUB"),stringsAsFactors = F)
                    }else{
                        tablaSearch <- rbind.data.frame(tablaSearch,c(seq, tablaSeqDes[seq,],goPerSeq[seq,],"NO-H"),stringsAsFactors = F)
                    }
                    
                }
        }
        
        
        tablaSearch <- tablaSearch[which(tablaSearch[,1] %in% row.names(datosEmbriones)),]
        colnames(tablaSearch) <- c("Seq ID","Seq Description","GO ID", "HUBGENE")
        tablaSearch
        
        

    })
    
    output$tableGeneralE <- DT::renderDataTable({
        table <- generalTableE()
        if (length(rownames(table)) > 0) {
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
    heatmapValues <- reactiveValues()
    observeEvent(input$fileHeatmapEmbryos,{
        inFile <- input$fileHeatmapEmbryos
        listSeq <- read.csv(inFile$datapath, header = T,stringsAsFactors = F)
        table <- datosEmbriones[listSeq[,1],]
        heatmapValues$table <- table
    })
    
#Descarga heatmap Embriones
output$downloadHeatEmbryos <- downloadHandler(
    filename = function(){paste("heatMapEmbryos",'.png',sep='')},
    content = function(file){
        ggsave(filename = file, plot = heatMapEmbryos(heatmapValues$table,input$clusterHeatE) )
        }
    )
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
            "$(this.api().table().header()).css({'background-color': '#1c1b1b', 'color': '#ffffff'});",
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
            "$(this.api().table().header()).css({'background-color': '#1c1b1b', 'color': '#ffffff'});",
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
                "$(this.api().table().header()).css({'background-color': '#1c1b1b', 'color': '#ffffff'});",
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
            "$(this.api().table().header()).css({'background-color': '#1c1b1b', 'color': '#ffffff'});",
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
            "$(this.api().table().header()).css({'background-color': '#1c1b1b', 'color': '#ffffff'});",
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
            "$(this.api().table().header()).css({'background-color': '#1c1b1b', 'color': '#ffffff'});",
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
    
    
    
})