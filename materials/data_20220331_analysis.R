rm(list = ls())
setwd("/home/useryk1/project/scRNA/huadongxu/data_20220331_analysis")
library(Seurat)
library(plyr)
library(reshape2)
library(dplyr)
library(tidyverse)
#library(AnnotationHub)
#library(patchwork)
#读取原始数据
# 分别读取每个10x样本的结果文件夹
samples=list.files("/home/useryk1/project/scRNA/huadongxu/data/",pattern = "^YKR")
samples
# pro=samples[1]
sceList = lapply(samples,function(pro){ 
  raw_sce <- CreateSeuratObject(Read10X( file.path('/home/useryk1/project/scRNA/huadongxu/data/',pro,'filtered_feature_bc_matrix') ),
                                project = pro, min.cells = 3, min.features = 200)
  raw_sce
})
sceList
sce.big <- merge(sceList[[1]], 
                 y = c(sceList[[2]] ,sceList[[3]], 
                       sceList[[4]]))
sce.big
#22784 features across 24136 samples within 1 assay
table(sce.big$orig.ident)
save(sce.big,file = 'sce.big.merge.Rdata')

## 合并成为一个R对象文件
load(file = 'sce.big.merge.Rdata')
#计算线粒体的含量
sce.big[["percent.mt"]] <- PercentageFeatureSet(object = sce.big, pattern = "^MT-")
#############################################################################
###QC之前
sce.big[["ident"]]="All_samples"
#利用小提琴图可视化QC的结果
Idents(object = sce.big) <- "ident"
pdf("00QC_before_metrics_violin_plot_all_samples.pdf")

VlnPlot(sce.big, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

VlnPlot(sce.big, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size=0, ncol = 3)

dev.off()

Idents(object = sce.big) <- "orig.ident"

pdf("00QC_before_metrics_violin_plot_bysamples.pdf",width = 20)

VlnPlot(sce.big, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

VlnPlot(sce.big, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size=0, ncol = 3)

dev.off()
###########################################################################

#过滤新立体含量小于25%，检测基因个数大于200小于4000的
sce.big <- subset(sce.big, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 25 &nCount_RNA >1000 & nCount_RNA < 20000)
sce.big
#22784 features across 20676 samples within 1 assay 
#############################################################################

###QC之后
#利用小提琴图可视化QC的结果
Idents(object = sce.big) <- "ident"
pdf("00QC_after_metrics_violin_plot_all_samples.pdf")

VlnPlot(sce.big, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

VlnPlot(sce.big, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size=0, ncol = 3)

dev.off()

Idents(object = sce.big) <- "orig.ident"

pdf("00QC_after_metrics_violin_plot_bysamples.pdf",width = 20)

VlnPlot(sce.big, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

VlnPlot(sce.big, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size=0, ncol = 3)

dev.off()
#################################################################################

#整合数据（CCA）
samples.list=SplitObject(sce.big, split.by = "orig.ident")
sceList=list()
for (i in 1:length(samples.list)) {
  sceList[[i]] <- NormalizeData(samples.list[[i]], verbose = FALSE)
  sceList[[i]] <- FindVariableFeatures(samples.list[[i]], selection.method = "vst", 
                                       nfeatures = 2000, verbose = FALSE)
}
sceList
sce.anchors <- FindIntegrationAnchors(object.list = sceList, dims = 1:30)
sce.integrated <- IntegrateData(anchorset = sce.anchors, dims = 1:30)
save(sce.integrated, file = 'sce.integrated.Rdata' )

sce.big<- sce.integrated

### 对RNA矩阵进行Normalize和Scale，以便进行差异表达和可视化
DefaultAssay(sce.big) <- "RNA"
mito.genes <- grep(pattern = "^MT-", x = rownames(x = sce.big), value = TRUE)
sce.big <- subset(sce.big,features=(rownames(x = sce.big)[!(rownames(x = sce.big) %in% mito.genes)]))
sce.big
#19730 features across 41215 samples within 2 assays 
sce.big <- NormalizeData(sce.big, normalization.method = "LogNormalize", scale.factor = 1e4)
sce.big <- FindVariableFeatures(sce.big, selection.method = "vst", nfeatures = 2000,verbose = FALSE)
sce.big <- ScaleData(sce.big, vars.to.regress = c("nCount_RNA", "percent.mt"),model.use = "linear")

### DefaultAssay设置为“integrated”矩阵并进行下游分析
DefaultAssay(sce.big) <- "integrated"
sce.big <- ScaleData(object = sce.big, vars.to.regress = c("nCount_RNA", "percent.mt"), model.use = "linear")
length(VariableFeatures(sce.big))
##特征提取：PCA线性降维
sce.big <- RunPCA(object = sce.big, pc.genes = VariableFeatures(sce.big), do.print = TRUE, pcs.print = 1:5, genes.print = 5)
sce.big <- ProjectDim(object = sce.big)

pdf("00dimheatmap_pca_10.pdf")
DimHeatmap(sce.big, dims = 1:10, cells = 500, balanced = TRUE)
dev.off()
###Check the PCs number
###鉴定数据集的可用维度
#JackStraw()函数中, 使用基于零分布的置换检验方法
#对PCA分析结果可以进行一系列的可视化： VizDimReduction, DimPlot, DimHeatmap
pdf("01jackStraw_dims_30s.pdf")
sce.big <- JackStraw(object = sce.big, num.replicate = 100,dims = 50)
sce.big <- ScoreJackStraw(object = sce.big, dims = 1:30)
plot1<-JackStrawPlot(sce.big, dims = 1:30)
plot1
dev.off()
### 通过ElbowPlot选择PC维度进行下游降维和聚类
pdf("01ElbowPlot_dims_30s.pdf")
plot2<-ElbowPlot(sce.big,ndims = 50)
plot2
dev.off()

pdf("01VizDimLoadings_dims_30s.pdf",width=15,height = 40)

VizDimLoadings(sce.big, dims = 1:30, reduction = "pca")

dev.off()

pdf("01DimHeatmap_dims_30s.pdf",height = 20)

DimHeatmap(sce.big, dims = 1:30, cells = 500, balanced = TRUE)

dev.off()

####Computing nearest neighbor graph, Computing SNN

sce.big <- FindNeighbors(sce.big, dims = 1:30)

sce.big <- FindClusters(sce.big, resolution = c(0.3,0.4,0.5,0.8,1.2))

head(sce.big@meta.data)

table(sce.big$integrated_snn_res.0.3)
table(sce.big$integrated_snn_res.0.4)
table(sce.big$integrated_snn_res.0.5)
table(sce.big$integrated_snn_res.0.8)
table(sce.big$integrated_snn_res.1.2)

###Run TSNE and UMAP

sce.big <- RunTSNE(object = sce.big,dims = 1:30, verbose = FALSE)

sce.big <- RunUMAP(object = sce.big,dims = 1:30, verbose = FALSE)

###res.0.3
Idents(object = sce.big) <- "integrated_snn_res.0.3"

pdf("02Dimplot_tsne_umap_pca_res0.3.pdf")

DimPlot(object = sce.big, reduction = "tsne")

DimPlot(object = sce.big, reduction = "umap")

DimPlot(object = sce.big, reduction = "pca")

dev.off()


###res.0.5

Idents(object = sce.big) <- "integrated_snn_res.0.5"

pdf("02Dimplot_tsne_umap_pca_res0.5.pdf")

DimPlot(object = sce.big, reduction = "tsne")

DimPlot(object = sce.big, reduction = "umap")

DimPlot(object = sce.big, reduction = "pca")

dev.off()


###res.0.8

Idents(object = sce.big) <- "integrated_snn_res.0.8"

pdf("02Dimplot_tsne_umap_pca_res0.8.pdf")

DimPlot(object = sce.big, reduction = "tsne")

DimPlot(object = sce.big, reduction = "umap")

DimPlot(object = sce.big, reduction = "pca")

dev.off()

###res.1.2

Idents(object = sce.big) <- "integrated_snn_res.1.2"

pdf("02Dimplot_tsne_umap_pca_res1.2.pdf")

DimPlot(object = sce.big, reduction = "tsne")

DimPlot(object = sce.big, reduction = "umap")

DimPlot(object = sce.big, reduction = "pca")

dev.off()

Idents(object = sce.big) <- "orig.ident"

pdf("02Dimplot_tsne_umap_pca_single_samples.pdf")

DimPlot(object = sce.big, reduction = "tsne")

DimPlot(object = sce.big, reduction = "umap")

DimPlot(object = sce.big, reduction = "pca")

dev.off()

head(sce.big@meta.data)

### 使用RNA矩阵进行marker gene的鉴定，而不是批次矫正后的integrated矩阵，
###此处需要注意，除Seurat之外，其他多种差异表达分析的算法也推荐使用原始表达值，而不是批次矫正后的数值；
###此外，RNA矩阵用作文章中作图数据来源。
###据此来看，integrated矩阵仅仅是在去除批次效应、降维和聚类过程中用到，差异表达与数据可视化使用的都是RNA矩阵。
DefaultAssay(sce.big) <- "RNA"

####DE analysis for resolution 0.3
Idents(object = sce.big) <- "integrated_snn_res.0.3"

sce_allmarkers_res.0.3 <- FindAllMarkers(object = sce.big, only.pos = TRUE, assay = 'RNA',min.pct = 0.1, logfc.threshold = 0.25, test.use = 'MAST')

sce.big@misc$allmarkers_res.0.3 <- sce_allmarkers_res.0.3

saveRDS(sce.big,file = "sce.big_DE.rds")

sce_allmarkers_res.0.3_arrange <-  sce.big@misc$allmarkers_res.0.3 %>% dplyr::arrange(cluster, p_val_adj)

write.table(data.frame(sce_allmarkers_res.0.3_arrange),"AllMarkerGenes_res.0.3_min.pct.0.1_logfc0.25.txt",sep="\t",row.names=F,quote=F)

sce_allmarkers_res.0.3_top10 <- sce_allmarkers_res.0.3_arrange %>% group_by(cluster) %>% top_n(10, avg_logFC)

write.table(data.frame(sce_allmarkers_res.0.3_top10),"Top10MarkerGenes_res.0.3_min.pct.0.1_logfc0.25.txt",sep="\t",row.names=F,quote=F)


sce.big <- readRDS("sce.big_DE.rds")
##细胞分群比例统计
pdf("05Percentage_of_cellnumbers_per_cluster.pdf",width = 9)
library(ggplot2)
t <- as.data.frame(table(sce.big@meta.data$integrated_snn_res.0.3))
dt<- data.frame(A=t$Freq,B=as.character(as.numeric(t$Var1)))
dt = dt[order(dt$A, decreasing = TRUE),]   
myLabel = as.vector(dt$B)   
myLabel = paste("cluster",myLabel, "(", round(dt$A / sum(dt$A) * 100, 2), "%)        ", sep = "")   
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

##每个cluster中细胞组成
pdf("05per_cluster_percent_of_sample.pdf",width = 9)
test<- as.data.frame(table(sce.big@meta.data$integrated_snn_res.0.3,sce.big@meta.data$orig.ident))
test$cluster <- paste( "cluster", as.numeric(test$Var1), sep = "") 
ggplot(test,aes(x= Var2,y=Freq,fill= cluster))+
  geom_col( position="fill")+
  labs(x='',y='Fraction of sample per cluster (%)')+
  theme(panel.background=element_rect(fill='transparent',color ="gray"),
        axis.text.x = element_text(hjust = 0.5, vjust =0.5,angle = 30,color = "black",size=9))
dev.off()

##cluster之间的相关性
pdf("05per_cluster_corr.pdf",width = 9)
tt <- log1p(AverageExpression(sce.big, verbose = FALSE)$RNA)
cormat <-round(cor(tt,method = "spearman"),2)
corrplot::corrplot(corr=cormat,method="color",order="AOE",addCoef.col="grey")
dev.off()

##Marker基因热图展示
clu_markers <- read.table("AllMarkerGenes_res.0.3_min.pct.0.1_logfc0.25.txt",sep="\t",header = TRUE)
clu_markers_top10 <- read.table("Top10MarkerGenes_res.0.3_min.pct.0.1_logfc0.25.txt",sep="\t",header = TRUE,stringsAsFactors = FALSE)
#tmp <- sample(colnames(sce.big),5000) %>% sort()
#sub_sce.big <-sce.big[,tmp] 
subobj <- subset(sce.big, downsample = 100)
DoHeatmap(subobj, features = as.character(head(unique(clu_markers_top10$gene[clu_markers_top10$gene %in% VariableFeatures(sce.big)]), 50)),disp.min=-2.5, disp.max=2.5,label = F,slot = "scale.data")
ggsave("04top10_doheatmap_resolution0.3.pdf",width = 9, height = 9)

sce_allmarkers_res.0.3_top2 <- sce.big@misc$allmarkers_res.0.3 %>% group_by(cluster) %>% top_n(2, avg_logFC)
DotPlot(subobj, features =as.character(unique(sce_allmarkers_res.0.3_top2$gene)),cols = c("lightgrey", "blue"))+ theme(axis.text.y = element_text(size = 7),axis.text.x = element_text(size = 7,angle = 90))
ggsave("04top2_dotheatmap_resolution0.3.pdf",width = 9, height = 5)

##marker基因展示
dir.create("marker_gene")
setwd("marker_gene")
for (j in 1:length(clu_markers_top10$gene)) {
  gene_name=clu_markers_top10$gene[j]
  VlnPlot(sce.big,gene_name, pt.size=0)
  ggsave(filename=paste0('VlnPlot_',gene_name,"_markers.pdf"))
  FeaturePlot(object = sce.big,gene_name)
  ggsave(filename=paste0('FeaturePlot_',gene_name,"_markers.pdf"))
}

