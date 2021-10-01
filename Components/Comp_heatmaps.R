library("pheatmap")

heatMapEmbryos <- function(table,clusterCut){
  table[table==0] <- Inf
  logMinValue <- log(min(table))
  table <- log10(table)
  table <- table + abs(logMinValue)
  table[table==Inf] <- 0
  table <- table[,c("ES1","ES2","ES3","PC","EC","C")]
  if (length(row.names(table))<=55) {
    pheatmap(table, cutree_rows = clusterCut,clustering_method = "ward.D2", cluster_cols = F)
  }else{
    pheatmap(table, cutree_rows = clusterCut,clustering_method = "ward.D2", cluster_cols = F,show_rownames = F)
  }

}

heatMapNeedles <- function(table,clusterCut){
  table[table==0] <- Inf
  logMinValue <- log(min(table))
  table <- log10(table)
  table <- table + abs(logMinValue)
  table[table==Inf] <- 0
  if (length(row.names(table))<=55) {
    pheatmap(table, cutree_rows = clusterCut,clustering_method = "ward.D2", cluster_cols = F)
  }else{
    pheatmap(table, cutree_rows = clusterCut,clustering_method = "ward.D2", cluster_cols = F,show_rownames = F)
  }
  
}

