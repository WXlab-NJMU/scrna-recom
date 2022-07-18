args<-commandArgs(T)
if (length(args) != 3) {
  cat("[usage:] Rscript script.R <genelist> <samplename> <outdir>\n")
  quit("no")
}

genelist <- args[1]
samplename <- args[2]
outdir <- args[3]

if ( !file.exists(outdir)) {
  dir.create(outdir)
}

setwd(outdir)

gene <- read.table(genelist,stringsAsFactors = F)

#注释
library(clusterProfiler)
library(org.Hs.eg.db)
df <- bitr(gene$V1, fromType = "SYMBOL", toType = c( "ENTREZID" ), OrgDb = org.Hs.eg.db )
target_gene_id <- df$ENTREZID
display_number = c(10, 10,10 )
## GO enrichment with clusterProfiler
ego_MF <- enrichGO(OrgDb="org.Hs.eg.db",
                   gene = target_gene_id,
                   pvalueCutoff = 0.4,
                   ont = "MF",
                   readable=TRUE)
ego_result_MF <- as.data.frame(ego_MF@result)[1:display_number[1], ]
write.table(as.data.frame(ego_MF@result),file=paste(outdir,paste(samplename,"GO_MF.xls",sep="_"),sep="/"),sep = "\t",quote = F)

ego_CC <- enrichGO(OrgDb="org.Hs.eg.db",
                   gene = target_gene_id,
                   pvalueCutoff = 0.05,
                   ont = "CC",
                   readable=TRUE)
ego_result_CC <- na.omit(as.data.frame(ego_CC@result)[1:display_number[2], ])
write.table(as.data.frame( ego_CC@result), file=paste(outdir,paste(samplename,"GO_CC.xls",sep="_"),sep="/"),sep = "\t",quote = F)

ego_BP <- enrichGO(OrgDb="org.Hs.eg.db",
                   gene = target_gene_id,
                   pvalueCutoff = 0.05,
                   ont = "BP",
                   readable=TRUE)
ego_result_BP <- na.omit(as.data.frame(ego_BP@result)[1:display_number[3], ])
write.table(as.data.frame( ego_BP@result), file=paste(outdir,paste(samplename,"GO_BP.xls",sep="_"),sep="/"),sep = "\t",quote = F)

options(stringsAsFactors = FALSE)
go_enrich_df <- data.frame(ID=c(ego_result_BP$ID, ego_result_CC$ID, ego_result_MF$ID),
                           Description=c(ego_result_BP$Description, ego_result_CC$Description, ego_result_MF$Description),
                           GeneNumber=c(ego_result_BP$Count, ego_result_CC$Count, ego_result_MF$Count),
                           type=factor(c(rep("biological process", display_number[3]), rep("cellular component", display_number[2]),
                                         rep("molecular function", display_number[1])), levels=c("molecular function", "cellular component", "biological process")))


## numbers as data on x axis
go_enrich_df$number <- factor(rev(1:nrow(go_enrich_df)))
## shorten the names of GO terms
shorten_names <- function(x, n_word=4, n_char=20){
  if (length(strsplit(x, " ")[[1]]) > n_word || (nchar(x) > 20))
  {
    if (nchar(x) > 20) x <- substr(x, 1, 20)
    x <- paste(paste(strsplit(x, " ")[[1]][1:min(length(strsplit(x," ")[[1]]), n_word)],
                     collapse=" "), "...", sep="")
    return(x)
  } 
  else
  {
    return(x)
  }
}

labels=lapply(go_enrich_df$Description,shorten_names)
names(labels) = rev(1:nrow(go_enrich_df))

library(ggplot2)
pdf(file=paste(outdir,paste(samplename,"go_enrichment_of_targets.pdf",sep="_"),sep="/"))
ggplot(data=go_enrich_df, aes(x=number, y=GeneNumber, fill=type)) +
  geom_bar(stat="identity", width=0.8) +
  coord_flip()+
  theme_bw() +
  scale_x_discrete(labels=labels) +
  xlab("GO term") +
  theme(axis.text.x=element_text(angle = 45,hjust=1, vjust=1,face = "bold", color="gray50")) +
  labs(title = "The Most Enriched GO Terms")+theme(plot.title = element_text(hjust = 0.5))

dev.off()

png(file=paste(outdir,paste(samplename,"go_enrichment_of_targets.png",sep="_"),sep="/"))
ggplot(data=go_enrich_df, aes(x=number, y=GeneNumber, fill=type)) +
  geom_bar(stat="identity", width=0.8) + 
  coord_flip()+
  theme_bw() + 
  scale_x_discrete(labels=labels) +
  xlab("GO term") + 
  theme(axis.text.x=element_text(angle = 45,hjust=1, vjust=1,face = "bold", color="gray50")) +
  labs(title = "The Most Enriched GO Terms")+theme(plot.title = element_text(hjust = 0.5))

dev.off()

#KEGG分析
gene<-na.omit(unique(df$ENTREZID))
kegg <- enrichKEGG(gene   = gene,   #KEGG分析
                   organism  = 'hsa',                 #物种
                   keyType = "kegg", 
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.2,
                   minGSSize = 1,
                   use_internal_data =FALSE)
write.table(as.data.frame(kegg@result), file=paste(outdir,paste(samplename,"kegg.txt",sep="_"),sep="/"),sep = "\t",quote = F)
#head(kegg)
#barplot(kegg, drop=TRUE, showCategory=12,title = "EnrichmentKEGG") #条形图
#dotplot(kegg,showCategory=10,title="Enrichmentkegg_dot","string")    #气泡图



#富集不显著柱状图###
nrow_kegg<-nrow(as.data.frame(kegg@result))
if(nrow_kegg>10)
{
	pathway<<-as.data.frame(kegg@result)[1:10,]
}else{

	pathway<<-as.data.frame(kegg@result)[1:nrow_kegg,]
}
pdf(file=paste(outdir,paste(samplename,"kegg_bar.pdf",sep="_"),sep="/"))
pathbar = ggplot(pathway,aes(x=Description,y=Count))
pathbar + geom_bar(stat="identity",aes(fill=pvalue),width=0.5) + coord_flip()+scale_fill_gradient( low = "blue", high = "red" )
dev.off()

png(file=paste(outdir,paste(samplename,"kegg_bar.png",sep="_"),sep="/"))
pathbar = ggplot(pathway,aes(x=Description,y=Count))
pathbar + geom_bar(stat="identity",aes(fill=pvalue),width=0.5) + coord_flip()+scale_fill_gradient( low = "blue", high = "red" ) 
dev.off()

#富集不显著气泡图#
pdf(file=paste(outdir,paste(samplename,"kegg_dot.pdf",sep="_"),sep="/"))
p = ggplot(pathway,aes(GeneRatio,Description ))
p=p + geom_point()
p=p + geom_point(aes(size=Count))
pbubble = p+ geom_point(aes(size=Count,color=p.adjust))
pr = pbubble+scale_color_gradient(low="blue",high = "red")
pr = pr+labs(color=expression(p.adjust),size="Count",x="GeneRatio",y="Description",title="")
pr + theme_bw()
dev.off()

png(file=paste(outdir,paste(samplename,"kegg_dot.png",sep="_"),sep="/"))
p = ggplot(pathway,aes(GeneRatio,Description ))
p=p + geom_point()
p=p + geom_point(aes(size=Count))
pbubble = p+ geom_point(aes(size=Count,color=p.adjust))
pr = pbubble+scale_color_gradient(low="blue",high = "red")
pr = pr+labs(color=expression(p.adjust),size="Count",x="GeneRatio",y="Description",title="")
pr + theme_bw()
dev.off()

