library(shiny)
library(xtable)
module_dataE <- read.csv("www/seq_module_file_Embryos.txt", sep = "\t", row.names = 1,col.names = c("seqName","module"), stringsAsFactors = F )
tablaSeqDes <- read.csv("www/seq_description.txt", sep = "\t",col.names = c("seq","description"), row.names = 1,stringsAsFactors = F)
tablaDiffEmbryos <-read.csv("www/embryo_differential.txt", sep = "\t", col.names = c("seqName","Comparative","Type","Log Fold Change"), stringsAsFactors = F )

tablaDiffNeedles <-read.csv("www/needles_differential.txt", sep = "\t", col.names = c("seqName","Condition","Log Fold Change"), stringsAsFactors = F )
module_dataN2 <- read.csv("www/needles_module.txt", sep = "\t",col.names = c("seqName","module"), stringsAsFactors = F )


moduleList <- function(){
  resList <- list()
  resList[["Module selection"]] <- "Module"
  for (module in unique(module_dataE$module)) {
    resList[module]<-module
  }
  return(resList)
}

corGenerationEmbryos <- function(listGenes, pvalueFilter, corFilter){
  NameVector <- c()
  sourceV <- c()
  targetV <- c()
  corValuesV <- c()
  signCor <- c()
  colorGene <- c()
  labelGene <- c()
  toolTipGene <- c()
  groupGene <- c()
  for (seq1 in names(listGenes)) {
    #Node configuration
    if (!(seq1 %in% row.names(module_dataE))) {
      colorGene <-c(colorGene,"gray")
      groupGene <- c(groupGene,"gray")
    }else{
      colorGene <-c(colorGene,module_dataE[seq1,])
      groupGene <- c(groupGene,module_dataE[seq1,])
    }
    #print(paste0("<b>", listGenes[[seq1]][["Name"]],"</b>"))
    labelGene <- c(labelGene,listGenes[[seq1]][["Name"]])
    diffTable <- tablaDiffEmbryos[which(tablaDiffEmbryos$seqName ==  seq1),c(2,3,4)]
    htmlTable <- print(xtable(diffTable, align="llll"), include.rownames=FALSE,
          type="html")
    toolTipGene <- c(toolTipGene,paste0("<p><b>", seq1,"</b><br>",tablaSeqDes[seq1,],"<br>",
                                        htmlTable,"</p>"))
    for (seq2 in names(listGenes)) {
      pasteSeq <- ""
      if (seq1 > seq2) {
        pasteSeq <- paste(seq1,seq2,sep = "")
      }else{
        if (seq2 > seq1) {
          pasteSeq <- paste(seq2,seq1,sep = "")
          
        }
      }
      
      if (seq1 != seq2 && !(pasteSeq %in% NameVector)) {
        NameVector <- c(NameVector,pasteSeq)
        corTest <- cor.test(as.numeric(listGenes[[seq1]]$expression),as.numeric(listGenes[[seq2]]$expression))
        corValue <- corTest$estimate
        pvalue <- corTest$p.value
        if (pvalue < pvalueFilter && abs(as.numeric(corValue))>corFilter) {
          sourceV <- c(sourceV, seq1 )
          targetV <- c(targetV, seq2 )
          corValuesV <- c(corValuesV, as.character(round(as.numeric(corValue), digits = 3)))
          if (as.numeric(corValue > 0)) {
            signCor <- c(signCor, "red")
          }else{
            signCor <- c(signCor, "blue")
          }
          
        }
        

      }
      
    }
  }
  
  nodes <- data.frame(id = names(listGenes), 
                      color = list(background = colorGene, border = rep("balck",length(names(listGenes)))), 
                      label= labelGene,
                      title = toolTipGene)
  
  edges <- data.frame(from = sourceV, to = targetV,
                      label = corValuesV,
                      color = signCor)
  return(list(nodes, edges))
}

corGenerationNeedles <- function(listGenes, pvalueFilter, corFilter){
  NameVector <- c()
  sourceV <- c()
  targetV <- c()
  corValuesV <- c()
  signCor <- c()
  colorGene <- c()
  labelGene <- c()
  toolTipGene <- c()
  groupGene <- c()
  for (seq1 in names(listGenes)) {
    #Node configuration
    if (!(seq1 %in% row.names(module_dataN2))) {
      colorGene <-c(colorGene,"gray")
      groupGene <- c(groupGene,"gray")
    }else{
      colorGene <-c(colorGene,module_dataN2[seq1,])
      groupGene <- c(groupGene,module_dataN2[seq1,])
    }
    
    labelGene <- c(labelGene,listGenes[[seq1]][["Name"]])
    diffTable <- tablaDiffNeedles[which(tablaDiffNeedles$seqName ==  seq1),c(2,3)]
    htmlTable <- print(xtable(diffTable, align="lll"), include.rownames=FALSE,
                       type="html")
    toolTipGene <- c(toolTipGene,paste0("<p><b>", seq1,"</b><br>",tablaSeqDes[seq1,],"<br>",
                                        htmlTable,"</p>"))
    for (seq2 in names(listGenes)) {
      pasteSeq <- ""
      if (seq1 > seq2) {
        pasteSeq <- paste(seq1,seq2,sep = "")
      }else{
        if (seq2 > seq1) {
          pasteSeq <- paste(seq2,seq1,sep = "")
          
        }
      }
      if (seq1 != seq2 && !(pasteSeq %in% NameVector)) {
        NameVector <- c(NameVector,pasteSeq)
        corTest <- cor.test(as.numeric(listGenes[[seq1]]$expression),as.numeric(listGenes[[seq2]]$expression))
        corValue <- corTest$estimate
        pvalue <- corTest$p.value
        if (pvalue < pvalueFilter && abs(as.numeric(corValue))>corFilter) {
          sourceV <- c(sourceV, seq1 )
          targetV <- c(targetV, seq2 )
          corValuesV <- c(corValuesV, as.character(round(as.numeric(corValue), digits = 3)))
          if (as.numeric(corValue > 0)) {
            signCor <- c(signCor, "red")
          }else{
            signCor <- c(signCor, "blue")
          }
          
        }
        
        
      }
      
    }
  }
  
  nodes <- data.frame(id = names(listGenes), 
                      color = list(background = colorGene, border = rep("balck",length(names(listGenes)))), 
                      label= labelGene,
                      title = toolTipGene)
  
  edges <- data.frame(from = sourceV, to = targetV,
                      label = corValuesV,
                      color = signCor)
  return(list(nodes, edges))
}
