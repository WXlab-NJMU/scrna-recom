rm(list = ls())
setwd("/home/useryk1/project/scRNA/huadongxu/data_20220331_analysis/DE_every_cluster")
library(Seurat)
library(plyr)
library(reshape2)
library(dplyr)
library(tidyverse)

sce.big <- readRDS("../celltype/sce.big_celltype.rds")
Idents(object = sce.big) <- "my_celltype_res0.3"
for( i in unique(sce.big@meta.data$my_celltype_res0.3)[!unique(sce.big@meta.data$my_celltype_res0.3) %in% "Leydig precursor cell"]){
  
  e<-gsub(" ", "_",i)
  
  markers_df<-diff_funtion(sce.big,i,"S","HP")
}

diff_funtion <- function(seuset,celltype,group1,group2){
  
  subpbmc <- subset(x = seuset,idents=celltype)
  
  markers_df <- FindMarkers(object = subpbmc, ident.1 = group1, ident.2 = group2,group.by = "group",min.pct = 0, logfc.threshold = 0)
  
  markers_df$gene_id <- rownames(markers_df)
  
  markers_df$change = as.factor(ifelse(markers_df$p_val_adj <= 0.05 & abs(markers_df$avg_logFC) > 1,
                                       ifelse(markers_df$avg_logFC > 1 , 'UP', 'DOWN' ), 'NOT' ) )
  
  e<-gsub(" ", "_",celltype)
  #print(e)
  write.table(markers_df,file=paste0(group1,"_vs_",group2,"_",e,'_markers.txt'),sep="\t",row.names=F,quote=F)
  
  this_tile <- paste0( 'The number of up gene is ', nrow(markers_df[ markers_df$change =='UP', ]),  "  ",
                       'The number of down gene is ', nrow(markers_df[ markers_df$change =='DOWN', ]) )
  volcano = ggplot(data = markers_df, aes( x = avg_logFC , y = -log10(p_val_adj), color = change)) +
    geom_point( alpha = 0.4, size = 1.75) +
    theme_set( theme_set( theme_bw( base_size = 15 ) ) ) +
    xlab( "avg_logFC" ) + ylab( "-log10 p-value" ) + theme( plot.title = element_text( size = 10, hjust = 0.5)) +
    scale_colour_manual( values = c('blue','black','red') )+ggtitle( this_tile ) 
  
  ggsave(volcano,filename=paste0(group1,"_vs_",group2,"_",e,"_volcano.pdf"))
  
  return(markers_df)
}
