createHeatRect <- function(color, sample){
  
  resRect <- column(HTML(
    paste("<img src='",sample,".png' width='120'' height='120'>
          ",
    "<svg width='120' height='120'
       viewBox='0 0 120 120'
       xmlns='http://www.w3.org/2000/svg'>
    <rect style='fill:",color,"; position: absolute' stroke='black' width=120 height=120></rect>
  </svg>",sep="")),width = 2)
  return(resRect)
}
datosEmbriones <- read.csv("www/CPM_embriones_media_filter.txt", sep = "\t", row.names = 1 )

#seqID <- "ppt_212074"
fun_color_range <- colorRampPalette(c("white","#fff700","red"))
tablaDiffEmbryos <-read.csv("www/embryo_differential.txt", sep = "\t", col.names = c("seqName","Comparative","Type","Log Fold Change"), stringsAsFactors = F )

tableDiffExpressionEmbryos <- function(seq){
  tablaRes <- data.frame()
  if (seq %in% tablaDiffEmbryos$seqName) {
    tablaRes <-  tablaDiffEmbryos[which(tablaDiffEmbryos$seqName ==  seq),]
    row.names(tablaRes) <- tablaRes$Comparative
    tablaRes$seqName <- NULL
    
  }  
  return(tablaRes)
}
my_colors <- fun_color_range(100)
heatMapEmbry <- function(seq){
  res <- ""
  if(seq %in% row.names(datosEmbriones)){
    expression <- datosEmbriones[seq,]
    minExpression <- min(expression)
    maxExpression <- max(expression)
    rangeExpression <- maxExpression-minExpression
    for (sample in colnames(expression)) {
      sample_expression <- datosEmbriones[seq,sample]
      sample_expression <- ((sample_expression-minExpression)/rangeExpression)*100
      if (sample_expression < 1) {
        sample_expression <- 1
      }
      color <- my_colors[sample_expression]
      
      res <- paste(res,HTML(toString(createHeatRect(color,sample))),sep = "\n")
    }
  }
  return(res)
}

legendColorsE <- function(seq){
  expression <- datosEmbriones[seq,]
  minExpression <- min(expression)
  maxExpression <- max(expression)
  rangeExpression <- maxExpression-minExpression
  svg(filename="www/legendE.svg",width = 20, height = 8)
  plot(rep(0,100),col=fun_color_range(100),font.lab = 2,cex.lab=8,pch=15,cex=20,yaxt='n',ylab = "",xlab = "CPM",xaxt = 'n')
  text(7,0,round(minExpression,digits = 2), cex= 6)
  text(95,0,round(maxExpression,digits = 2), cex= 6)
  dev.off()
  #res <- div(tags$img(src="legendR.svg",width="300px",height="180px"))
  #return(res)
}
