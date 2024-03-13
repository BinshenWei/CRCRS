library(maftools)


setwd("/Users/weibinshen/Desktop/未命名文件夹/TMB/gdc_download_20240113_091547.950190")
files <- list.files(pattern = '*.gz',recursive = TRUE)
all_mut <- data.frame()
for (file in files) {
  mut <- read.delim(file,skip = 7, header = T, fill = TRUE,sep = "\t")
  all_mut <- rbind(all_mut,mut)
}

all_mut <- read.maf(all_mut)
##计算TMB
tmb <- tmb(maf =all_mut ,
           captureSize = 50,
           logScale = T)
head(tmb)
write.xlsx(tmb,'TMB.xlsx')

#计算mutant-allele tumor heterogeneity
barcode <- unique(all_mut@data$Tumor_Sample_Barcode)
head(barcode)
MATH <- data.frame()
for (i in barcode){
  out.math = inferHeterogeneity(maf =all_mut , tsb = i)
  Tumor_Sample_Barcode=unique(out.math$clusterData$Tumor_Sample_Barcode)
  m = unique(out.math$clusterData$MATH)
  out = data.frame(Tumor_Sample_Barcode, m)
  MATH = rbind(MATH, out)
}
head(MATH)
write.xlsx(MATH,'MATH.xlsx')