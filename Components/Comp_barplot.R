library(ggplot2)
datosEmbriones <- read.csv("www/CPM_embriones_media_filter.txt", sep = "\t", row.names = 1 ,stringsAsFactors = F)
expressionNeedles <- read.csv("www/Needles_CPM_mean.txt", sep = "\t", row.names = 1 )

#seq <- "ppt_212074"

barplotEmbryosFunc <- function(seq){
  barData <- data.frame(stringsAsFactors = F)
  expression <- datosEmbriones[seq,]
  barData <- rbind(expression[,c("PC","EC","C")])
  colnames(barData) <- c("PC-ES1","EC-ES2","C-ES3")
  barData <- rbind(barData,as.numeric(expression[,c("ES1","ES2","ES3")]))
  row.names(barData)<-c("Cygot","Somatic")

  barplot(as.matrix(barData), main="CPM expression comparative",font.lab = 2,cex.lab=1.2,cex.names = 2, font = 2,
          xlab="Comparatives",ylab = "CPM expression",ylim = c(0,(max(expression) + (min(expression)/2))), density=c(30,30),angle =c(45) ,col=c("darkblue","red"),
          legend = rownames(barData),args.legend = list(text.font = 2),beside=TRUE)
}

RC_data <- read.csv("www/CPM_raizMicrodiseccion_RC_mean.txt", sep = "\t", row.names = 1, stringsAsFactors = F )
RDV_data <- read.csv("www/CPM_raizMicrodiseccion_RDV_mean.txt", sep = "\t", row.names = 1, stringsAsFactors = F )
RDC_data <- read.csv("www/CPM_raizMicrodiseccion_RDC_mean.txt", sep = "\t", row.names = 1, stringsAsFactors = F )
RM_data <- read.csv("www/CPM_raizMicrodiseccion_RM_mean.txt", sep = "\t", row.names = 1, stringsAsFactors = F )
#seq <- "pp_117813"
barplotRootsFunc <- function(seq){
  RC_C <- RC_data[seq,"RC_C"]
  RC_N <- RC_data[seq,"RC_N"]
  RDC_C <- RDC_data[seq,"RDC_C"]
  RDC_N <- RDC_data[seq,"RDC_N"]
  RDV_C <- RDV_data[seq,"RDV_C"]
  RDV_N <- RDV_data[seq,"RDV_N"]
  RM_C <- RM_data[seq,"RM_C"]
  RM_N <- RM_data[seq,"RM_N"]
  expression <- c(RC_C,RDC_C,RDV_C,RM_C,RC_N,RM_N,RDC_N,RDV_N)
  barData <- data.frame(stringsAsFactors = F)
  barData <- rbind.data.frame(c(RC_N,RM_N,RDC_N,RDV_N), stringsAsFactors = F)
  colnames(barData) <- c("RC_N-RC_C","RM_N-RM_C","RDC_N-RDC_C","RDV_N-RDV_C")
  barData <- rbind.data.frame(barData,c(RC_C,RM_C,RDC_C,RDV_C), stringsAsFactors = F)
  row.names(barData)<-c("Ammonium (3mM)","Control")
  rangeExpression <- 25*((max(expression)-min(expression))/100)
  barplot(as.matrix(barData), main="CPM expression comparative",font.lab = 2,cex.lab=1.2,cex.names = 0.6, font = 2,
          xlab="Comparatives",ylab = "CPM expression",ylim = c(0,(max(expression) + (rangeExpression))), density=c(30,30),angle =c(45) ,col=c("darkblue","red"),
          legend = rownames(barData),args.legend = list(text.font = 2),beside=TRUE)
}

barplotNeedlesFunc <- function(seq){
  expression <- expressionNeedles[seq,c("N0","M0","N1","M1","N3","M3")]
  barData <- data.frame(stringsAsFactors = F)
  barData <- rbind.data.frame(as.numeric(expressionNeedles[seq,c("N0","N1","N3")]), stringsAsFactors = F)
  colnames(barData) <- c("W0","W1","W3")
  barData <- rbind.data.frame(barData,as.numeric(expressionNeedles[seq,c("M0","M1","M3")]), stringsAsFactors = F)
  row.names(barData)<-c("November","May")
  rangeExpression <- 25*((max(expression)-min(expression))/100)
  barplot(as.matrix(barData), main="CPM expression comparative",font.lab = 2,cex.lab=1.2,cex.names = 0.6, font = 2,
          xlab="Comparatives",ylab = "CPM expression",ylim = c(0,(max(expression) + (rangeExpression))), density=c(30,30),angle =c(45) ,col=c("darkblue","red"),
          legend = rownames(barData),args.legend = list(text.font = 2),beside=TRUE)
}
