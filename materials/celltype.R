rm(list = ls())
setwd("/home/useryk1/project/scRNA/huadongxu/data_20220331_analysis/celltype")
library(Seurat)
library(plyr)
library(reshape2)
library(dplyr)
library(tidyverse)
##肿瘤常用分类marker
celltype_marker=c(
  "EPCAM",#上皮细胞 epithelial
  "PECAM1",#内皮细胞 endothelial
  "COL3A1",#成纤维细胞 fibroblasts
  "CD163","AIF1",#髓系细胞 myeloid
  "CD79A",#B细胞
  "JCHAIN",#浆细胞 plasma cell
  "CD3D","CD8A","CD4",#T细胞
  "GNLY","NKG7",#NK细胞
  "PTPRC"#免疫细胞
)
VlnPlot(sce.big,features = celltype_marker,pt.size = 0,ncol = 3)

##利用scCATCH进行细胞类型注释
sce.big <- readRDS("../sce.big_DE.rds")
all_clu_markers <- sce.big@misc$allmarkers_res.0.3 
clu_markers_10 <- all_clu_markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
library(scCATCH)
clu_ann <- scCATCH(clu_markers_10,
                   species = "Human",
                   cancer = NULL,
                   tissue = "Blood")
sce.big@misc$cellcytype_res_0.3_top10 <- clu_ann
celltype <- data.frame(ClusterID=clu_ann$cluster, celltype=clu_ann$cell_type, stringsAsFactors = F)
write.table(clu_ann,"celltype_top10_res_0.3.txt",sep="\t",row.names=F,quote=F)

#分群重新命名
clu_ann <- read.table("celltype_final.txt",header = T,sep= "\t",stringsAsFactors =F)
celltype <- data.frame(ClusterID=clu_ann$cluster, celltype=clu_ann$cell_type, stringsAsFactors = F)
sce.big@meta.data$my_celltype_res0.3="unknown"
for(i in 1:nrow(celltype)){
  sce.big@meta.data[which(sce.big@meta.data$integrated_snn_res.0.3 == celltype$ClusterID[i]),'my_celltype_res0.3'] <- celltype$celltype[i]}

Idents(object = sce.big) <- "integrated_snn_res.0.3"
pdf("04Dimplot_tsne_umap_pca_celltype_res0.3.pdf")
DimPlot(sce.big, reduction = "tsne", label = TRUE, pt.size = 0.5,group.by  = "my_celltype_res0.3" )
DimPlot(sce.big, reduction = "umap", label = TRUE, pt.size = 0.5,group.by  = "my_celltype_res0.3")
DimPlot(sce.big, reduction = "pca", label = TRUE, pt.size = 0.5,group.by  = "my_celltype_res0.3")
dev.off()

########################################
##添加分组信息
table(sce.big@meta.data$orig.ident)
sce.big@meta.data$group="S"
sce.big@meta.data[grep("YKR3259",sce.big$orig.ident),]$group="HP"
sce.big@meta.data[grep("YKR3260",sce.big$orig.ident),]$group="HP"

saveRDS(sce.big,file = "sce.big_celltype.rds")
Idents(object = sce.big) <- "my_celltype_res0.3"
##每组中cluster细胞组成
pdf("05per_group_percent_of_cluster_celltype0.3.pdf",width = 9)
test<- as.data.frame(table(sce.big@meta.data$my_celltype_res0.3,sce.big@meta.data$group))
ggplot(test,aes(x= Var2,y=Freq,fill= Var1))+
  geom_col( position="fill")+
  labs(x='',y='Fraction of sample per cluster (%)')+
  theme(panel.background=element_rect(fill='transparent',color ="gray"),
        axis.text.x = element_text(hjust = 0.5, vjust =0.5,angle = 30,color = "black",size=9))+
  scale_x_discrete(limits=c("S","HP"))
dev.off()

##cluster之间的相关性
pdf("05per_cluster_corr_celltype0.3.pdf",width = 9)
tt <- log1p(AverageExpression(sce.big, verbose = FALSE)$RNA)
cormat <-round(cor(tt,method = "spearman"),2)
corrplot::corrplot(corr=cormat,method="color",order="AOE",addCoef.col="grey")
dev.off()

##细胞分群比例统计
pdf("05Percentage_of_cellnumbers_per_cluster_celltype0.3.pdf",width = 9)
library(ggplot2)
t <- as.data.frame(table(sce.big@meta.data$my_celltype_res0.3))
dt<- data.frame(A=t$Freq,B=as.character(t$Var1))
dt = dt[order(dt$A, decreasing = TRUE),]   
myLabel = as.vector(dt$B)   
myLabel = paste(myLabel, "(", round(dt$A / sum(dt$A) * 100, 2), "%)        ", sep = "")   
p = ggplot(dt, aes(x = "", y = A, fill = B)) + 
  geom_bar(stat = "identity", width = 1) +    
  coord_polar(theta = "y") + 
  labs(x = "", y = "", title = "") + 
  theme(axis.ticks = element_blank()) + 
  theme(legend.title = element_blank(), legend.position = "right") + 
  scale_fill_discrete(breaks = dt$B, labels = myLabel) + 
  theme(axis.text.x = element_blank())   
p
dev.off()



