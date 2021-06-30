library(shiny)

RC_data <- read.csv("www/CPM_raizMicrodiseccion_RC_mean.txt", sep = "\t", row.names = 1, stringsAsFactors = F )
RDV_data <- read.csv("www/CPM_raizMicrodiseccion_RDV_mean.txt", sep = "\t", row.names = 1, stringsAsFactors = F )
RDC_data <- read.csv("www/CPM_raizMicrodiseccion_RDC_mean.txt", sep = "\t", row.names = 1, stringsAsFactors = F )
RM_data <- read.csv("www/CPM_raizMicrodiseccion_RM_mean.txt", sep = "\t", row.names = 1, stringsAsFactors = F )
tablaDiff <-read.csv("www/microdiseccion_diffExpression.txt", sep = "\t", col.names = c("seqName","Tissue","Condition","Log Fold Change"), stringsAsFactors = F )


fun_color_range <- colorRampPalette(c("white","#fff700","red"))
rescaleColor <- function(rangeExpression,minExpression,expression){
  my_colors <- fun_color_range(100)
  value <- ((expression-minExpression)/rangeExpression)*100
  if (value < 1) {
    value <- 1
  }
  color <- my_colors[value]
  return(color)
}

tableDiffExpressionRoots <- function(seq){
  tablaRes <- data.frame()
  if (seq %in% tablaDiff$seqName) {
    tablaRes <-  tablaDiff[which(tablaDiff$seqName ==  seq),]
    row.names(tablaRes) <- tablaRes$Tissue
    tablaRes$seqName <- NULL

  }  
  return(tablaRes)
}

heatMapRootsC <- function(seq){
  res <- ""
  if(seq %in% c(row.names(RC_data),row.names(RM_data),row.names(RDC_data),row.names(RDV_data))){
    RC_C <- RC_data[seq,"RC_C"]
    RC_N <- RC_data[seq,"RC_N"]
    RDC_C <- RDC_data[seq,"RDC_C"]
    RDC_N <- RDC_data[seq,"RDC_N"]
    RDV_C <- RDV_data[seq,"RDV_C"]
    RDV_N <- RDV_data[seq,"RDV_N"]
    RM_C <- RM_data[seq,"RM_C"]
    RM_N <- RM_data[seq,"RM_N"]
    expression <- c(RC_C,RDC_C,RDV_C,RM_C,RC_N,RM_N,RDC_N,RDV_N)
    minExpression <- min(expression)
    maxExpression <- max(expression)
    rangeExpression <- maxExpression-minExpression
    
    my_colors <- fun_color_range(rangeExpression)
    svgRootFile_C  <- readLines("www/svg.svg")
    RM_C_Color <- rescaleColor(rangeExpression,minExpression,RM_C)
    svgRootFile_C  <- gsub(pattern = "#b180b1", replace = RM_C_Color, x = svgRootFile_C)
    RDV_C_Color <- rescaleColor(rangeExpression,minExpression,RDV_C)
    svgRootFile_C  <- gsub(pattern = "#3c4bcb", replace = RDV_C_Color, x = svgRootFile_C)
    RDC_C_Color <- rescaleColor(rangeExpression,minExpression,RDC_C)
    svgRootFile_C  <- gsub(pattern = "#24b44c", replace = RDC_C_Color, x = svgRootFile_C)
    RC_C_Color <- rescaleColor(rangeExpression,minExpression,RC_C)
    svgRootFile_C  <- gsub(pattern = "#7494bc", replace = RC_C_Color, x = svgRootFile_C)
    writeLines(svgRootFile_C, con="www/raiz_heat_modC.svg")
    #res <- div(tags$img(src="raiz_heat_modC.svg",width="280px",height="660px",style="margin-right:-40px"))
  }
  return(res)
}

heatMapRootsN <- function(seq){
  res<- ""
  if(seq %in% c(row.names(RC_data),row.names(RM_data),row.names(RDC_data),row.names(RDV_data))){
    RC_C <- RC_data[seq,"RC_C"]
    RC_N <- RC_data[seq,"RC_N"]
    RDC_C <- RDC_data[seq,"RDC_C"]
    RDC_N <- RDC_data[seq,"RDC_N"]
    RDV_C <- RDV_data[seq,"RDV_C"]
    RDV_N <- RDV_data[seq,"RDV_N"]
    RM_C <- RM_data[seq,"RM_C"]
    RM_N <- RM_data[seq,"RM_N"]
    expression <- c(RC_C,RDC_C,RDV_C,RM_C,RC_N,RM_N,RDC_N,RDV_N)
    minExpression <- min(expression)-1
    maxExpression <- max(expression)+1
    rangeExpression <- maxExpression-minExpression
    
    svgRootFile_N  <- readLines("www/svg.svg")
    RM_N_Color <- rescaleColor(rangeExpression,minExpression,RM_N)
    svgRootFile_N  <- gsub(pattern = "#b180b1", replace = RM_N_Color, x = svgRootFile_N)
    RDV_N_Color <- rescaleColor(rangeExpression,minExpression,RDV_N)
    svgRootFile_N  <- gsub(pattern = "#3c4bcb", replace = RDV_N_Color, x = svgRootFile_N)
    RDC_N_Color <- rescaleColor(rangeExpression,minExpression,RDC_N)
    svgRootFile_N  <- gsub(pattern = "#24b44c", replace = RDC_N_Color, x = svgRootFile_N)
    RC_N_Color <- rescaleColor(rangeExpression,minExpression,RC_N)
    svgRootFile_N  <- gsub(pattern = "#7494bc", replace = RC_N_Color, x = svgRootFile_N)
    writeLines(svgRootFile_N, con="www/raiz_heat_modN.svg")
    #res <- div(tags$img(src="raiz_heat_modN.svg",width="280px",height="660px",style="margin-right:-40px"))
  }
  
  return(res)
}

legendColorsR <- function(seq){
  RC_C <- RC_data[seq,"RC_C"]
  RC_N <- RC_data[seq,"RC_N"]
  RDC_C <- RDC_data[seq,"RDC_C"]
  RDC_N <- RDC_data[seq,"RDC_N"]
  RDV_C <- RDV_data[seq,"RDV_C"]
  RDV_N <- RDV_data[seq,"RDV_N"]
  RM_C <- RM_data[seq,"RM_C"]
  RM_N <- RM_data[seq,"RM_N"]
  expression <- c(RC_C,RDC_C,RDV_C,RM_C,RC_N,RM_N,RDC_N,RDV_N)
  minExpression <- min(expression)
  maxExpression <- max(expression)
  rangeExpression <- maxExpression-minExpression
  svg(filename="www/legendR.svg",width = 20, height = 8)
  plot(rep(0,100),col=fun_color_range(100),font.lab = 2,cex.lab=8,pch=15,cex=20,yaxt='n',ylab = "",xlab = "CPM",xaxt = 'n')
  text(7,0,round(minExpression,digits = 2), cex= 6)
  text(95,0,round(maxExpression,digits = 2), cex= 6)
  dev.off()
  #res <- div(tags$img(src="legendR.svg",width="300px",height="180px"))
  #return(res)
}

plot(rep(1,100),col=fun_color_range(100),pch=19,cex=3, yaxt='n',ylab = "",xlab = "CPM",xaxt = 'n')

axis(side=1,at=c(0,100),labels=c(1,100))
