
tablaDiff <-read.csv("www/needles_differential.txt", sep = "\t", col.names = c("seqName","Condition","Log Fold Change"), stringsAsFactors = F )
expressionNeedles <- read.csv("www/Needles_CPM_mean.txt", sep = "\t", row.names = 1 )


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

tableDiffExpressionNeedles <- function(seq){
  tablaRes <- data.frame()
  if (seq %in% tablaDiff$seqName) {
    tablaRes <-  tablaDiff[which(tablaDiff$seqName ==  seq),]
    row.names(tablaRes) <- tablaRes$Tissue
    tablaRes$seqName <- NULL
    
  }  
  return(tablaRes)
}


heatMapNeedlesN <- function(seq){
  res <- ""
  if(seq %in% c(row.names(expressionNeedles))){

    expression <- expressionNeedles[seq,]
    minExpression <- min(expression)
    maxExpression <- max(expression)
    rangeExpression <- maxExpression-minExpression
    
    my_colors <- fun_color_range(rangeExpression)
    svgNeedleFile_N  <- readLines("www/needle.svg")
    W0_N_Color <- rescaleColor(rangeExpression,minExpression,expressionNeedles[seq,c("N0")])
    svgNeedleFile_N  <- gsub(pattern = "#1bf224", replace = W0_N_Color, x = svgNeedleFile_N)
    W1_N_Color <- rescaleColor(rangeExpression,minExpression,expressionNeedles[seq,c("N1")])
    svgNeedleFile_N  <- gsub(pattern = "#3c4ccc", replace = W1_N_Color, x = svgNeedleFile_N)
    W3_N_Color <- rescaleColor(rangeExpression,minExpression,expressionNeedles[seq,c("N3")])
    svgNeedleFile_N  <- gsub(pattern = "#548536", replace = W3_N_Color, x = svgNeedleFile_N)
    writeLines(svgNeedleFile_N, con="www/needle_N.svg")
    #res <- div(tags$img(src="raiz_heat_modC.svg",width="280px",height="660px",style="margin-right:-40px"))
  }
  return(res)
}
heatMapNeedlesM <- function(seq){
  res <- ""
  if(seq %in% c(row.names(expressionNeedles))){
    
    expression <- expressionNeedles[seq,]
    minExpression <- min(expression)
    maxExpression <- max(expression)
    rangeExpression <- maxExpression-minExpression
    
    my_colors <- fun_color_range(rangeExpression)
    svgNeedleFile_N  <- readLines("www/needle.svg")
    W0_N_Color <- rescaleColor(rangeExpression,minExpression,expressionNeedles[seq,c("M0")])
    svgNeedleFile_N  <- gsub(pattern = "#1bf224", replace = W0_N_Color, x = svgNeedleFile_N)
    W1_N_Color <- rescaleColor(rangeExpression,minExpression,expressionNeedles[seq,c("M1")])
    svgNeedleFile_N  <- gsub(pattern = "#3c4ccc", replace = W1_N_Color, x = svgNeedleFile_N)
    W3_N_Color <- rescaleColor(rangeExpression,minExpression,expressionNeedles[seq,c("M3")])
    svgNeedleFile_N  <- gsub(pattern = "#548536", replace = W3_N_Color, x = svgNeedleFile_N)
    writeLines(svgNeedleFile_N, con="www/needle_M.svg")
    #res <- div(tags$img(src="raiz_heat_modC.svg",width="280px",height="660px",style="margin-right:-40px"))
  }
  return(res)
}