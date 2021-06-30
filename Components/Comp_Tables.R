
tablaGoDes <- read.csv("www/go_des_file.txt", sep = "\t",col.names = c("go","descripcion"), row.names = 1 ,stringsAsFactors = F)
tablaSeqGoR <- read.csv("www/seq_gos_raiz.txt", sep = "\t",col.names = c("seqName","GoID"),stringsAsFactors = F)
tablaSeqDes <- read.csv("www/seq_description.txt", sep = "\t",col.names = c("seq","description"), row.names = 1,stringsAsFactors = F)
tablaSeqGo <- read.csv("www/seq_go_embrios.txt", sep = "\t",col.names = c("seqName","GoID"),stringsAsFactors = F)
datosEmbriones <- read.csv("www/CPM_embriones_media_filter.txt", sep = "\t", row.names = 1 )
tablaSeqDesRoots <- read.csv("www/seq_description_raiz.txt", sep = "\t",col.names = c("seq","description"), row.names = 1,stringsAsFactors = F)
RM_data <- read.csv("www/CPM_raizMicrodiseccion_RM_mean.txt", sep = "\t", row.names = 1, stringsAsFactors = F )

#seqID <- "ppt_212074"
goTableDes <- function(seq){
  goTable <- data.frame(stringsAsFactors = F)
  tablaSearch <- data.frame(stringsAsFactors = F)
  if (seq %in% tablaSeqGo$seqName) {
    tableSearch <- tablaSeqGo[which(tablaSeqGo$seqName == seq),]
  }
  for (go in tableSearch$GoID) {
    goTable <- rbind(goTable,c(go,tablaGoDes[go,]),stringsAsFactors = F)
    colnames(goTable) <- c("GO id", "GO descriptor")
  }
  
  return(goTable)
}
goTableDesR <- function(seq){
  
  goTable <- data.frame(stringsAsFactors = F)
  tablaSearch <- data.frame(stringsAsFactors = F)
  if (seq %in% tablaSeqGoR$seqName) {
    tableSearch <- tablaSeqGoR[which(tablaSeqGoR$seqName == seq),]
    for (go in tableSearch$GoID) {
      goTable <- rbind(goTable,c(go,tablaGoDes[go,]),stringsAsFactors = F)
    }
    colnames(goTable) <- c("GO id", "GO descriptor")
  }
  return(goTable)
}

#go <- "GO:0006032"
searchGoSeq <- function(go){
  res <- NULL
  tablaSearch <- data.frame(stringsAsFactors = F)
  if (go %in% tablaSeqGo$GoID) {
    tablaSearch <- tablaSeqGo[which(tablaSeqGo$GoID == go),]

    for (seq in tablaSearch$seqName) {
      if (seq %in% row.names(datosEmbriones)) {
        res <- rbind.data.frame(res,c(seq,tablaSeqDes[seq,]),stringsAsFactors = F)
      }
      
    }
    colnames(res) <- c("Sequence ID", "NCBI description blast")
  }

  return(res)
}


searchGoSeqR <- function(go){
  res <- NULL
  tablaSearch <- data.frame(stringsAsFactors = F)
  if (go %in% tablaSeqGoR$GoID) {
    tablaSearch <- tablaSeqGoR[which(tablaSeqGoR$GoID == go),]
    
    for (seq in tablaSearch$seqName) {
      if (seq %in% row.names(RM_data)) {
        res <- rbind.data.frame(res,c(seq,tablaSeqDesRoots[seq,]),stringsAsFactors = F)
      }
      
    }
    colnames(res) <- c("Sequence ID", "NCBI description blast")
  }
  
  return(res)
}

