setwd("/Users/weibinshen/Desktop/未命名文件夹/体细胞突变")
library(openxlsx)
#library(questionr)
library(epiDisplay)

input1 = "all_mut_01.csv"
groupf = "group.xlsx"

data = read.csv(input1,check.names=F)
group.data = read.xlsx(groupf)
high.sample = group.data[group.data[,2]=="High",1]
low.sample = group.data[group.data[,2]=="Low",1]
high.sample.pos = match(high.sample,colnames(data))
low.sample.pos = match(low.sample,colnames(data))

high.mut = apply(data[,high.sample.pos],1,function(x)sum(x!=0))
high.nonmut = apply(data[,high.sample.pos],1,function(x)sum(x==0))
low.mut = apply(data[,low.sample.pos],1,function(x)sum(x!=0))
low.nonmut = apply(data[,low.sample.pos],1,function(x)sum(x==0))

newdata = data.frame(high.mut=high.mut,high.nonmut=high.nonmut,low.mut=low.mut,low.nonmut=low.nonmut)
rownames(newdata) = data[,1]
gene_or_p = t(apply(newdata,1,function(x){
  cc = cci(x[1],x[2],x[3],x[4])
  or = cc$or
  ct = chisq.test(cc$table,correct=F)
  pvalue = ct$p.value
  return(c(or,pvalue))}))
colnames(gene_or_p) = c("or","pvalue")
gene_or_p = apply(gene_or_p,2,as.numeric)
## filter
gene_or_p[is.na(gene_or_p)]=1
outdata = data.frame(gene = data[,1],or=gene_or_p[,1],pvalue=gene_or_p[,2])
#write.table(outdata,"gene_or_pvalue.xls",sep="\t",col.names=T,row.names=F,quote=F)
################## draw  ####
outdata = outdata[order(outdata$pvalue),]
outdata = outdata[outdata$or>1,]
outdata = outdata[outdata$pvalue<0.05,]
write.xlsx(outdata,'体细胞突变.xlsx')
nr = 20
df = head(outdata,n=nr)
df =read.xlsx('体细胞突变.xlsx')
genelist = df$gene
mutdata = data[match(genelist,data[,1]),]
low_mutdata = mutdata[,match(low.sample,colnames(mutdata))]
high_mutdata = mutdata[,match(high.sample,colnames(mutdata))]
draw.data = cbind(low_mutdata,high_mutdata)
rownames(draw.data) = genelist
nc = ncol(draw.data)

pdf("mutation.pdf",width=20)

layout(mat=matrix(c(1,2),nc=2),width=c(7,3))
par(mar=c(2,6,1,0),xpd=T)
plot(1,type="n",xlim=c(0,nc),ylim=c(-nr,3),xaxt="n",yaxt="n",xaxs="i",yaxs="i",axes=F,xlab="",ylab="")
rect(0,-nr,nc,0,col="grey90",border="grey90") # bg 
for(i in 1:nr){
  muti = as.numeric(draw.data[i,])
  coli = ifelse(muti==1,"black",NA)
  borderi=coli
  rect(0:(nc-1),-i,1:nc,-i+1,col=coli,border=borderi)
  text(-2,-i+0.5,rownames(draw.data)[i],cex=1.3,adj=1,xpd=T)
  #text(nc+0.5,-i+0.5,paste("Pvalue=",round(df[i,3],4)," OR=",round(df[i,2],2)),cex=1,adj=0)
}
rect(0,0.5,length(low.sample),2.5,col="#4575B4",border="#4575B4")  # low
rect(length(low.sample),0.5,nc,2.5,col="#D73027",border="#D73027") # high
text(-2,1.5,"Low",adj=1,cex=1.5) # low
#text(nc,1.5,"High",adj=0,cex=1.5)

### right
par(mar=c(2,0.1,1,1))
plot(1,type="n",xlim=c(0,10),ylim=c(-nr,3),xaxt="n",yaxt="n",xaxs="i",yaxs="i",axes=F,xlab="",ylab="")
for(i in 1:nr){
  text(0.1,-i+0.5,paste("Pvalue=",round(df[i,3],4)," OR=",round(df[i,2],2)),cex=1,adj=0)
}
text(0.1,1.5,"High",adj=0,cex=1.5)
legend("topright",legend=c("NO","YES"),title="Mutation",
       title.cex=1.5,bty="n",pch=15,col=c("grey90","black"))
dev.off()