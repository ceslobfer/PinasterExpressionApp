datosEmbriones <- read.csv("www/CPM_embriones_media_filter.txt", sep = "\t", row.names = 1, stringsAsFactors = F )
tablaSeqDes <- read.csv("www/seq_description.txt", sep = "\t",col.names = c("seq","description"), row.names = 1,stringsAsFactors = F)
RC_data <- read.csv("www/CPM_raizMicrodiseccion_RC_mean.txt", sep = "\t", row.names = 1, stringsAsFactors = F )
RDV_data <- read.csv("www/CPM_raizMicrodiseccion_RDV_mean.txt", sep = "\t", row.names = 1, stringsAsFactors = F )
RDC_data <- read.csv("www/CPM_raizMicrodiseccion_RDC_mean.txt", sep = "\t", row.names = 1, stringsAsFactors = F )
RM_data <- read.csv("www/CPM_raizMicrodiseccion_RM_mean.txt", sep = "\t", row.names = 1, stringsAsFactors = F )
module_dataE <- read.csv("www/seq_module_file_Embryos.txt", sep = "\t",col.names = c("seqName","module"), stringsAsFactors = F )
tablaSeqDesRoots <- read.csv("www/seq_description_raiz.txt", sep = "\t",col.names = c("seq","description"), row.names = 1,stringsAsFactors = F)
hubsEmbrio <- read.csv("www/hubModule.txt", sep = "\t",col.names = c("seqName","module"),stringsAsFactors = F)

#seq <- "ppt_55700"
seqCorEmbryos <- function(seq){
  tablaCor <- data.frame(stringsAsFactors = F)
  seqQuery <- as.numeric(datosEmbriones[seq,])
  moduleSeq <- module_dataE[which(module_dataE$seqName == seq),2]
  print(moduleSeq)
  tablaModule <- module_dataE[which(module_dataE$module == moduleSeq),]
  for (seq2 in tablaModule$seqName) {
    if (seq2!=seq) {
      correlation <- cor.test(seqQuery,as.numeric(datosEmbriones[seq2,]))
      if (correlation$p.value < 0.05 & seq2 %in% hubsEmbrio$seqName){
        tablaCor <- rbind.data.frame(tablaCor,c(seq2, tablaSeqDes[seq2,],correlation$estimate,correlation$p.value,TRUE),stringsAsFactors = F)
      }
      if (correlation$p.value < 0.05 & !(seq2 %in% hubsEmbrio$seqName)){
        tablaCor <- rbind.data.frame(tablaCor,c(seq2, tablaSeqDes[seq2,],correlation$estimate,correlation$p.value,FALSE),stringsAsFactors = F)
      }
    }
  }
  
  colnames(tablaCor) <- c("Sequence.ID","Blast.description","Correlation.value","Correlation.pvalue","HUBGENE")
  tablaCor <- transform(tablaCor, Correlation.value = as.numeric(Correlation.value))
  tablaCor <- transform(tablaCor, Correlation.pvalue = as.numeric(Correlation.pvalue))
  tablaCor <- tablaCor[order(abs(tablaCor$Correlation.value), decreasing = TRUE),]
  return(tablaCor)
}



expressionRoots <- function(seq){
  RC_C <- RC_data[seq,"RC_C"]
  RC_N <- RC_data[seq,"RC_N"]
  RDC_C <- RDC_data[seq,"RDC_C"]
  RDC_N <- RDC_data[seq,"RDC_N"]
  RDV_C <- RDV_data[seq,"RDV_C"]
  RDV_N <- RDV_data[seq,"RDV_N"]
  RM_C <- RM_data[seq,"RM_C"]
  RM_N <- RM_data[seq,"RM_N"]
  expression <- c(RC_C,RDC_C,RDV_C,RM_C,RC_N,RM_N,RDC_N,RDV_N)
  return(expression)
}

#seq <- "pp_100138"
seqCorRoots <- function(seq){
  tablaCor <- data.frame(stringsAsFactors = F)
  seqQuery <- as.numeric(expressionRoots(seq))
  for (seq2 in row.names(RM_data)) {
    if (seq2!=seq) {
      correlation <- cor.test(seqQuery,as.numeric(expressionRoots(seq2)))
      if (correlation$p.value < 0.05){
        tablaCor <- rbind.data.frame(tablaCor,c(seq2, tablaSeqDesRoots[seq2,],correlation$estimate,correlation$p.value),stringsAsFactors = F)
      }
    }
  }
  colnames(tablaCor) <- c("Sequence.ID","Blast.description","Correlation.value","Correlation.pvalue")
  tablaCor <- transform(tablaCor, Correlation.value = as.numeric(Correlation.value))
  tablaCor <- transform(tablaCor, Correlation.pvalue = as.numeric(Correlation.pvalue))
  tablaCor <- tablaCor[order(abs(tablaCor$Correlation.value), decreasing = TRUE),]
  return(tablaCor)
}