library(survival)
library(plyr)
##制作table1
library(tableone)
##制作森林图
library(forestplot)
setwd("/Users/weibinshen/Desktop/未命名文件夹/亚组")
#4.查看数据前6行

PLCO<-read.xlsx('TCGA的CRCRS分为高低.xlsx')
head(PLCO)
#5.查看数据数据性质
str(PLCO)
PLCO$Status<-factor(PLCO$Status)
summary(PLCO$Status)
categorical_vars <- c("Gender", "M","T","N",'Age','CRCRS')
PLCO[, categorical_vars] <- lapply(PLCO[, categorical_vars], as.factor)
#1-建一个名为"HR","5%CI","95%CI","P"的4列表：a表。
#用于放批量获取的结果。
a<-c("HR","5%CI","95%CI","P")
##2建两个for循环，第一个是批量Cox亚组分析，
##第二个是将结果批量提取放入a表
for(i in 5:dim(PLCO)[2]){ 
  cox1<-coxph(Surv(OS,Status==1)~CRCRS,subset=(PLCO[,i]=='1'),data=PLCO)
  aa<-summary(cox1)
  for(j in 1:dim(aa$coefficients)[1]){
    a<-rbind(a,c(round(aa$coefficients[j,2],2),
                 round(aa$conf.int[j,3],2),
                 round(aa$conf.int[j,4],2),
                 round(aa$coefficients[j,5],3)))}}

#3-a表加入第一列名字，
#即变量T/ER/PR/KI67/LVI/N/Grade中亚变量名为1的实际含义
rownames(a)<-c("Characteristics",
               "Female",
               "Age<=70",
               'I-II',
               "M0",
               "N0",
               "T1");a <- data.frame(a)
a$HR.CI95<-paste0(a$X1," (",a$X2,"-",a$X3,")");a


#1- 建立b表
b<-c("HR","5%CI","95%CI","P")

#2- 所有变量为2的亚组分析for循环
for(i in 5:dim(PLCO)[2]){ 
  cox1<-coxph(Surv(OS,Status==1)~CRCRS,subset=(PLCO[,i]=='2'),data=PLCO)
  bb<-summary(cox1)
  for(j in 1:dim(bb$coefficients)[1]){
    b<-rbind(b,c(round(bb$coefficients[j,2],2),
                 round(bb$conf.int[j,3],2),
                 round(bb$conf.int[j,4],2),
                 round(bb$coefficients[j,5],3)))}}

#3- b表第一列名字，即变量T/ER/PR/KI67/LVI/N/Grade中亚变量名为2的实际含义
rownames(b)<-c("Characteristics",
               "Male",
               "Age>70",
               "III-IV",
               "M1",
               "N1",
               "T2")
b <- data.frame(b);b

#4- 加入HR(95% CI)格式
b$HR.CI95<-paste0(b$X1," (",b$X2,"-",b$X3,")");b


#1- 建立c表
c<-c("HR","5%CI","95%CI","P")

#2- 亚变量为3的亚组分析
for(i in 9:dim(PLCO)[2]){ 
  cox2<-coxph(Surv(OS,Status==1)~CRCRS,subset=(PLCO[,i]=='3'),data=PLCO)
  cc<-summary(cox2)
  for(j in 1:dim(cc$coefficients)[1]){
    c<-rbind(c,c(round(cc$coefficients[j,2],2),
                 round(cc$conf.int[j,3],2),
                 round(cc$conf.int[j,4],2),
                 round(cc$coefficients[j,5],3)))}}

#3- c表的第一列名字，即变量N/Grade中亚变量名为3的实际含义
rownames(c)<-c("Characteristics",
               "N2",
               "T3")
c <- data.frame(c);c
#4- 加入HR（95% CI）列
c$HR.CI95<-paste0(c$X1," (",c$X2,"-",c$X3,")");c


d<-c("HR","5%CI","95%CI","P")

#2- 亚变量为3的亚组分析
for(i in 10:dim(PLCO)[2]){ 
  cox2<-coxph(Surv(OS,Status==1)~CRCRS,subset=(PLCO[,i]=='4'),data=PLCO)
  dd<-summary(cox2)
  for(j in 1:dim(dd$coefficients)[1]){
    d<-rbind(d,c(round(dd$coefficients[j,2],2),
                 round(dd$conf.int[j,3],2),
                 round(dd$conf.int[j,4],2),
                 round(dd$coefficients[j,5],3)))}}

#3- c表的第一列名字，即变量N/Grade中亚变量名为3的实际含义
rownames(d)<-c("Characteristics",
               "T4")
d <- data.frame(d);d
#4- 加入HR（95% CI）列
d$HR.CI95<-paste0(d$X1," (",d$X2,"-",d$X3,")");d


#1- RT单因素分析
cox3<-coxph(Surv(OS,Status==1)~CRCRS,data=PLCO)
EE<-summary(cox3)

#2- 提取信息
EE$coefficients 
EE$conf.int    
HR<- round(EE$coefficients[,2],2)     
PValue <- round(EE$coefficients[,5],3) 
CI5<- round(EE$conf.int[,3],2)
CI95<- round(EE$conf.int[,4],2)
HR.CI95<-paste0(HR," (",CI5,"-",CI95,")")

#3- E表是最后建的，就是把信息提取出来组合为d表
e<- data.frame("X1" = HR,
               "X2" = CI5,
               "X3" = CI95,
               "X4" = PValue,
               "HR.CI95" = HR.CI95)
#4- 行名，命名为所有人群
rownames(e)<-"All parents"           



#1- 查看数据变量名
names(PLCO)
#"status" "time""RT" "AGE""T" "ER" "PR" "KI67""LVI" "N" "Grade" 
#2- 指定需进行基线统计的变量，纳入的变量为上面进行的亚组分析的变量
myVars <- c( "Gender" , "Age"    , "Stage"  , "M"  ,    "N"   ,    "T"      )
#指定基线表中哪些变量是分类变量
catVars<-c("Gender" , "Age"    , "Stage"  , "M"  ,    "N"   ,    "T"      )
#3- 构建table1函数，RT为分类标准
table<- print(CreateTableOne(vars=myVars,
                             data=PLCO,
                             strata="CRCRS",
                             factorVars=catVars),
              catDigits = 2,contDigits = 2,showAllLevels=TRUE)     
#4- 去掉第4.5列（上图p和test列），新表命名为table1（17行3列）
table1<-table[,c(-4,-5)]
#1- 按照table1的格式将abcd表合成一个表
res1<-rbind(e[1,],#abcd表合为res1表
            a[2,],b[2,],#T的亚变量
            a[3,],b[3,],#ER
            a[4,],b[4,],#PR
            a[5,],b[5,],#KI67
            a[6,],b[6,],c[2,],#LVI
            a[7,],b[7,],c[3,],d[2,])#N
#2- res1表列名进入表格内，并命名为Characteristics
res1<-tibble::rownames_to_column(res1, var = "Characteristics")
res2<-cbind(table1,res1)
res3<-res2[,-1]
#5- 将合成的res3进行列排序，亚变量名放首位
res4<-data.frame(res3[c(3,2,1,4:6,8,7)])

N<-rbind(res4[1, ], #提取res4表第一行
         c("Gender",rep(NA, 7)),#插入一行第一个单元格为T stage，后7列为NA
         res4[c(2:3), ],#提取res4表第2，3行
         c("Age",rep(NA, 7)),#插入一行第一个单元格为ER stage，后7列为NA
         res4[c(4:5),],#以此类推
         c("Stage",rep(NA, 7)),
         res4[c(6:7),],
         c("Distant metastasis",rep(NA, 7)), 
         res4[c(8:9),],
         c("Lymph node metastasis",rep(NA, 7)),
         res4[c(10:12),],
         c("Depth of invasion",rep(NA, 7)), 
         res4[c(13:16),])
#3- p=0的变为<0.001
#2-插入行名
for(i in 2:8) {N[, i] = as.character(N[, i])}#先让变量性质变为character类型
result <- rbind( c("Characteristics","High CRCRS","Low CRCRS",NA,NA,NA,"HR (95%CI)","P Value"),
                 N[c(1:22),])
#3- p=0的变为<0.001
result$X4[result$X4==0]<-"<0.001"
#4- N分期放在T分期下面
result<-data.frame(result,row.names=NULL);result
result <- rbind(result [c(1:11,19:23,15:18,12:14),])



fig<-forestplot(result[,c(1,2,3,7,8)], #12378列显示为原数字格式
                mean=result[,4],   
                lower=result[,5],  
                upper=result[,6], 
                zero=1,            
                boxsize=0.4,      
                graph.pos= 4,#图放在第四列
                hrzl_lines=list("1" = gpar(lty=1,lwd=2),
                                "2" = gpar(lty=2),
                                "24"= gpar(lwd=2,lty=1)),
                graphwidth = unit(.25,"npc"),
                xlab="HR<1,意味着放疗能改善局部区域复发风险",
                xticks=c(0,1,8) ,
                #----------------#字体
                is.summary=c(T,
                             F,
                             T,F,F,
                             T,F,F,
                             T,F,F,
                             T,F,F,F,F,
                             T,F,F,F,
                             T,F,F
                ),#T=粗体
                txt_gp=fpTxtGp(label=gpar(cex=1),
                               ticks=gpar(cex=1.1), 
                               xlab=gpar(cex=1), 
                               title=gpar(cex=2)),
                #----------------#线条粗细（x轴、置信区间）
                lwd.zero=2,
                lwd.ci=2,
                lwd.xaxis=1, 
                lty.ci=1,
                ci.vertices =T,
                ci.vertices.height=0.2, 
                clip=c(0,3),
                #----------------#行间距、字间距/box形状                 
                ineheight=unit(8, 'mm'), 
                line.margin=unit(8, 'mm'),
                colgap=unit(6, 'mm'),
                col=fpColors(zero = "black",
                             box = 'black', 
                             lines = 'black'),
                fn.ci_norm="fpDrawCircleCI", 
                title="亚组分析森林图")
fig




PLCO<-read.xlsx('PLCO.xlsx')
head(PLCO)
#5.查看数据数据性质
str(PLCO)
PLCO$Status<-factor(PLCO$Status)
summary(PLCO$Status)
categorical_vars <- c("Gender", "M","T","N",'Age','CRCRS')
PLCO[, categorical_vars] <- lapply(PLCO[, categorical_vars], as.factor)
#1-建一个名为"HR","5%CI","95%CI","P"的4列表：a表。
#用于放批量获取的结果。
a<-c("HR","5%CI","95%CI","P")
##2建两个for循环，第一个是批量Cox亚组分析，
##第二个是将结果批量提取放入a表
for(i in 5:dim(PLCO)[2]){ 
  cox1<-coxph(Surv(OS,Status==1)~CRCRS,subset=(PLCO[,i]=='1'),data=PLCO)
  aa<-summary(cox1)
  for(j in 1:dim(aa$coefficients)[1]){
    a<-rbind(a,c(round(aa$coefficients[j,2],2),
                 round(aa$conf.int[j,3],2),
                 round(aa$conf.int[j,4],2),
                 round(aa$coefficients[j,5],3)))}}

#3-a表加入第一列名字，
#即变量T/ER/PR/KI67/LVI/N/Grade中亚变量名为1的实际含义
rownames(a)<-c("Characteristics",
               "Female",
               "Age<=70",
               'I-II',
               "M0",
               "N0",
               "T1");a <- data.frame(a)
a$HR.CI95<-paste0(a$X1," (",a$X2,"-",a$X3,")");a


#1- 建立b表
b<-c("HR","5%CI","95%CI","P")

#2- 所有变量为2的亚组分析for循环
for(i in 5:dim(PLCO)[2]){ 
  cox1<-coxph(Surv(OS,Status==1)~CRCRS,subset=(PLCO[,i]=='2'),data=PLCO)
  bb<-summary(cox1)
  for(j in 1:dim(bb$coefficients)[1]){
    b<-rbind(b,c(round(bb$coefficients[j,2],2),
                 round(bb$conf.int[j,3],2),
                 round(bb$conf.int[j,4],2),
                 round(bb$coefficients[j,5],3)))}}

#3- b表第一列名字，即变量T/ER/PR/KI67/LVI/N/Grade中亚变量名为2的实际含义
rownames(b)<-c("Characteristics",
               "Male",
               "Age>70",
               "III-IV",
               "M1",
               "N1",
               "T2")
b <- data.frame(b);b

#4- 加入HR(95% CI)格式
b$HR.CI95<-paste0(b$X1," (",b$X2,"-",b$X3,")");b


#1- 建立c表
c<-c("HR","5%CI","95%CI","P")

#2- 亚变量为3的亚组分析
for(i in 9:dim(PLCO)[2]){ 
  cox2<-coxph(Surv(OS,Status==1)~CRCRS,subset=(PLCO[,i]=='3'),data=PLCO)
  cc<-summary(cox2)
  for(j in 1:dim(cc$coefficients)[1]){
    c<-rbind(c,c(round(cc$coefficients[j,2],2),
                 round(cc$conf.int[j,3],2),
                 round(cc$conf.int[j,4],2),
                 round(cc$coefficients[j,5],3)))}}

#3- c表的第一列名字，即变量N/Grade中亚变量名为3的实际含义
rownames(c)<-c("Characteristics",
               "N2",
               "T3")
c <- data.frame(c);c
#4- 加入HR（95% CI）列
c$HR.CI95<-paste0(c$X1," (",c$X2,"-",c$X3,")");c


d<-c("HR","5%CI","95%CI","P")

#2- 亚变量为3的亚组分析
for(i in 10:dim(PLCO)[2]){ 
  cox2<-coxph(Surv(OS,Status==1)~CRCRS,subset=(PLCO[,i]=='4'),data=PLCO)
  dd<-summary(cox2)
  for(j in 1:dim(dd$coefficients)[1]){
    d<-rbind(d,c(round(dd$coefficients[j,2],2),
                 round(dd$conf.int[j,3],2),
                 round(dd$conf.int[j,4],2),
                 round(dd$coefficients[j,5],3)))}}

#3- c表的第一列名字，即变量N/Grade中亚变量名为3的实际含义
rownames(d)<-c("Characteristics",
               "T4")
d <- data.frame(d);d
#4- 加入HR（95% CI）列
d$HR.CI95<-paste0(d$X1," (",d$X2,"-",d$X3,")");d


#1- RT单因素分析
cox3<-coxph(Surv(OS,Status==1)~CRCRS,data=PLCO)
EE<-summary(cox3)

#2- 提取信息
EE$coefficients 
EE$conf.int    
HR<- round(EE$coefficients[,2],2)     
PValue <- round(EE$coefficients[,5],3) 
CI5<- round(EE$conf.int[,3],2)
CI95<- round(EE$conf.int[,4],2)
HR.CI95<-paste0(HR," (",CI5,"-",CI95,")")

#3- E表是最后建的，就是把信息提取出来组合为d表
e<- data.frame("X1" = HR,
               "X2" = CI5,
               "X3" = CI95,
               "X4" = PValue,
               "HR.CI95" = HR.CI95)
#4- 行名，命名为所有人群
rownames(e)<-"All parents"           



#1- 查看数据变量名
names(PLCO)
#"status" "time""RT" "AGE""T" "ER" "PR" "KI67""LVI" "N" "Grade" 
#2- 指定需进行基线统计的变量，纳入的变量为上面进行的亚组分析的变量
myVars <- c( "Gender" , "Age"    , "Stage"  , "M"  ,    "N"   ,    "T"      )
#指定基线表中哪些变量是分类变量
catVars<-c("Gender" , "Age"    , "Stage"  , "M"  ,    "N"   ,    "T"      )
#3- 构建table1函数，RT为分类标准
table<- print(CreateTableOne(vars=myVars,
                             data=PLCO,
                             strata="CRCRS",
                             factorVars=catVars),
              catDigits = 2,contDigits = 2,showAllLevels=TRUE)     
#4- 去掉第4.5列（上图p和test列），新表命名为table1（17行3列）
table1<-table[,c(-4,-5)]
#1- 按照table1的格式将abcd表合成一个表
res1<-rbind(e[1,],#abcd表合为res1表
            a[2,],b[2,],#T的亚变量
            a[3,],b[3,],#ER
            a[4,],b[4,],#PR
            a[5,],b[5,],#KI67
            a[6,],b[6,],c[2,],#LVI
            a[7,],b[7,],c[3,],d[2,])#N
#2- res1表列名进入表格内，并命名为Characteristics
res1<-tibble::rownames_to_column(res1, var = "Characteristics")
res2<-cbind(table1,res1)
res3<-res2[,-1]
#5- 将合成的res3进行列排序，亚变量名放首位
res4<-data.frame(res3[c(3,2,1,4:6,8,7)])

N<-rbind(res4[1, ], #提取res4表第一行
         c("Gender",rep(NA, 7)),#插入一行第一个单元格为T stage，后7列为NA
         res4[c(2:3), ],#提取res4表第2，3行
         c("Age",rep(NA, 7)),#插入一行第一个单元格为ER stage，后7列为NA
         res4[c(4:5),],#以此类推
         c("Stage",rep(NA, 7)),
         res4[c(6:7),],
         c("Distant metastasis",rep(NA, 7)), 
         res4[c(8:9),],
         c("Lymph node metastasis",rep(NA, 7)),
         res4[c(10:12),],
         c("Depth of invasion",rep(NA, 7)), 
         res4[c(13:16),])
#3- p=0的变为<0.001
#2-插入行名
for(i in 2:8) {N[, i] = as.character(N[, i])}#先让变量性质变为character类型
result <- rbind( c("Characteristics","High CRCRS","Low CRCRS",NA,NA,NA,"HR (95%CI)","P Value"),
                 N[c(1:22),])
#3- p=0的变为<0.001
result$X4[result$X4==0]<-"<0.001"
#4- N分期放在T分期下面
result<-data.frame(result,row.names=NULL);result
result <- rbind(result [c(1:11,19:23,15:18,12:14),])



fig<-forestplot(result[,c(1,2,3,7,8)], #12378列显示为原数字格式
                mean=result[,4],   
                lower=result[,5],  
                upper=result[,6], 
                zero=1,            
                boxsize=0.4,      
                graph.pos= 4,#图放在第四列
                hrzl_lines=list("1" = gpar(lty=1,lwd=2),
                                "2" = gpar(lty=2),
                                "24"= gpar(lwd=2,lty=1)),
                graphwidth = unit(.25,"npc"),
                xlab="HR<1,意味着放疗能改善局部区域复发风险",
                xticks=c(0,1,7) ,
                #----------------#字体
                is.summary=c(T,
                             F,
                             T,F,F,
                             T,F,F,
                             T,F,F,
                             T,F,F,F,F,
                             T,F,F,F,
                             T,F,F
                ),#T=粗体
                txt_gp=fpTxtGp(label=gpar(cex=1),
                               ticks=gpar(cex=1.1), 
                               xlab=gpar(cex=1), 
                               title=gpar(cex=2)),
                #----------------#线条粗细（x轴、置信区间）
                lwd.zero=2,
                lwd.ci=2,
                lwd.xaxis=1, 
                lty.ci=1,
                ci.vertices =T,
                ci.vertices.height=0.2, 
                clip=c(0,3),
                #----------------#行间距、字间距/box形状                 
                ineheight=unit(8, 'mm'), 
                line.margin=unit(8, 'mm'),
                colgap=unit(6, 'mm'),
                col=fpColors(zero = "black",
                             box = 'black', 
                             lines = 'black'),
                fn.ci_norm="fpDrawCircleCI", 
                title="亚组分析森林图")
fig