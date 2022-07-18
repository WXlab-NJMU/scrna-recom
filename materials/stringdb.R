##protein_stringdb##
setwd("/home/useryk1/project/scRNA/huadongxu/data_20220331_analysis/ppi")
library(STRINGdb)
string_db <- STRINGdb$new(version='10',species =9606) ###小鼠10090#构建对象## 大鼠10116#人类是9606
data <- read.table("total.gene",header = T,sep="\t")
mapped <- string_db$map(data,"gene",removeUnmappedRows=TRUE)#用map函数将基因id转化string id#
hits <- mapped$STRING_id[1:100]# 选择前100个
#mapped1<-mapped[1:100,]
#up_loc<-which(mapped1$logFC>0)
#down_loc<-which(mapped1$logFC<0)
#up.len<-length(up_loc)
#down.len<-length(down_loc)
#up_data<-mapped1[up_loc,]
#down_data<-mapped1[down_loc,]
#up.color<-colorRampPalette(c("red","white"))(up.len)
#down.color<-colorRampPalette(c("green","white"))(down.len)
#example_mapped_pval05<-rbind(data.frame(up_data,color=up.color),data.frame(down_data,color=down.color))
#payload_id<-string_db$post_payload(example_mapped_pval05$STRING_id,
#                                   colors = example_mapped_pval05$color)
#pdf(file="mRNA_string_PPI.pdf")
#string_db$plot_network(hits,payload_id = payload_id)
#dev.off()
#提取interaction，cystocape绘制关系图#
info <- string_db$get_interactions(hits)
write.table(mapped,file="id_mapp.txt",sep="\t",quote = F,row.names = F)
write.table(info,file="string_inter.txt",sep="\t",quote = F,row.names = F)
