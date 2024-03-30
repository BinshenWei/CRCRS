setwd("/Users/weibinshen/Desktop/未命名文件夹/富集分析")
TCGA_COAD_TPM<-read.table('exprSet_Count.txt', header = TRUE, sep = "\t", check.names = FALSE)

# 读取 Excel 文件
library(readxl)
library(limma)
mycounts <- as.data.frame(avereps(TCGA_COAD_TPM[,-1],ID = TCGA_COAD_TPM$gene_name))
mycounts[1:5,1:5]


group <- read_excel("group.xlsx")
group$CRCRS <- as.factor(group$CRCRS)
group$CRCRS <- relevel(group$CRCRS, ref = "Low")
library(DESeq2)
dds <-DESeqDataSetFromMatrix(countData=round(mycounts), 
                             colData=group, 
                             design=~CRCRS)
dds <- DESeq(dds)
resultsNames(dds)

res <- results(dds, name="CRCRS_High_vs_Low")
#res <- res[res$padj < 0.05 & abs(res$log2FoldChange) >=1, ]
#res
#plotMA(res)
rownamesres<-res@rownames
res <- cbind(res, rownamesres)
res<- as.data.frame(res)
library(dplyr)
library(knitr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
#排序
data_sort <- res %>%
  arrange(desc(log2FoldChange))
head(data_sort)
tail(data_sort)
gene_list <- data_sort$log2FoldChange
geneID<-data_sort$rownamesres
Gene.id <- bitr(
  geneID, fromType = "SYMBOL",
  toType = c("ENTREZID"),
  OrgDb = org.Hs.eg.db
) 
gene_list <- as.numeric(gene_list)

names(gene_list) <- Gene.id$SYMBOL
head(gene_list)

#GO
res <- gseGO(
  gene_list,    
  ont = "ALL",    
  OrgDb = org.Hs.eg.db,    
  keyType = "SYMBOL", 
  minGSSize = 100,
  maxGSSize = 200,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",    
)

head(res,2)
dotplot(res,
        split=".sign") + facet_grid(.~.sign)


a=res@result
write.xlsx(a,'GO.xlsx')
#KRGG
gene_list <- data_sort$log2FoldChange

gene_list_KEGG <- as.numeric(gene_list)

names(gene_list_KEGG) <- Gene.id$ENTREZID
head(gene_list_KEGG)

res_KEGG <- gseKEGG(
  gene_list_KEGG,    
  organism = "hsa", 
  minGSSize = 40,
  maxGSSize = 200,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)
dotplot(
  res_KEGG,showCategory=15,
  split=".sign") + facet_grid(.~.sign)
b=res_KEGG
write.xlsx(b,'kEGG.xlsx')
