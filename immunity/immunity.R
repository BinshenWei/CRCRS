setwd("/Users/weibinshen/Desktop/未命名文件夹/免疫")
source('CIBERSORT.R')
TCGA_TME.results<-CIBERSORT('LM22.txt','exprSet_Count.txt', perm = 1000, QN = F) #perm置换次数=1000，QN分位数归一化=TRU
phenotype = read.xlsx("group.xlsx")
TCGA_TME.results<-read.table("CIBERSORT-Results.txt",header = TRUE, sep = "\t", check.names = FALSE)
cibersort_data <- as.data.frame(TCGA_TME.results[,1:24])
group_list <- phenotype$CRCRS
table(group_list)
cibersort_data$group <- group_list
#cibersort_data$sample <- row.names(cibersort_data)
library(reshape2)

# 2.2 融合数据
cibersort_data_New = melt(cibersort_data)
## Using group, sample as id variables

colnames(cibersort_data_New)=c("Sample","Group","Celltype","Composition")  #设置行名
head(cibersort_data_New)

if(T){
  mytheme <- theme(plot.title = element_text(size = 12,color="black",hjust = 0.5),
                   axis.title = element_text(size = 12,color ="black"), 
                   axis.text = element_text(size= 12,color = "black"),
                   panel.grid.minor.y = element_blank(),
                   panel.grid.minor.x = element_blank(),
                   axis.text.x = element_text(angle = 45, hjust = 1 ),
                   panel.grid=element_blank(),
                   legend.position = "top",
                   legend.text = element_text(size= 12),
                   legend.title= element_text(size= 12)
  ) }
library(ggplot2)
library(ggpubr)
box_TME <- ggplot(cibersort_data_New, aes(x = Celltype, y = Composition))+ 
  labs(y="Cell composition",x= NULL,title = "TME Cell composition")+  
  geom_boxplot(aes(fill = Group),position=position_dodge(0.5),width=0.5,outlier.alpha = 0)+ 
  scale_fill_manual(values = c("#1CB4B8", "#EB7369"))+
  theme_classic() + mytheme + 
  stat_compare_means(aes(group =  Group),
                     label = "p.signif",
                     method = "wilcox.test",
                     hide.ns = T)



box_TME
box_TME$theme
