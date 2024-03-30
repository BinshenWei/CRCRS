setwd("/Users/weibinshen/Desktop/未命名文件夹/nomagram")
PLCO<-read.xlsx('PLCO.xlsx')
TCGA<-read.xlsx('TCGA.xlsx')
categorical_vars <- c("Gender", "M","T","N",'Age')
PLCO[, categorical_vars] <- lapply(PLCO[, categorical_vars], as.factor)
TCGA[, categorical_vars] <- lapply(TCGA[, categorical_vars], as.factor)
#train nomogram
library(survival)
library(rms)
library(timeROC)
library(survival)
library(caret)


dd <- datadist(PLCO)

options(datadist = "dd")

model <- cph(Surv(OS, Status) ~ CRCRS+T+N+M+Age,data = PLCO ,x = TRUE, y = T,surv = T)
surv_OS <- Survival(model)
nom_OS <- nomogram(model, fun = list(function(x) surv_OS(730, x),
                                     function(x) surv_OS(1095, x),
                                     function(x) surv_OS(1825, x)
),
fun.at = c(seq(.1, .9, by = .1), .95),
funlabel = c("2-year probability of OS","3-year probability of OS","5-year probability of OS"),lp=0)
plot(nom_OS)

test.ph <- cox.zph(model)
test.ph
ggplot <- ggcoxzph(test.ph, linear.predictions = TRUE, pval = TRUE) +
  labs(x = "Time", y = "Scaled Schoenfeld Residuals")
ggplot
train_nomogram <- predict(model)
nomogram_PLCO <- timeROC(T = PLCO$OS,
                                        delta = PLCO$Status,
                                        marker = train_nomogram,
                                        cause = 1,
                                        times = c(365*2,365*3,365 * 5),
                                        iid = TRUE)
nomogram_PLCO$AUC
nomogram_PLCO$iid

#Age
model_Age <- cph(Surv(OS, Status) ~ Age,data = PLCO ,x = TRUE, y = T,surv = T)
train_Age <- predict(model_Age)
nomogram_Age <- timeROC(T = PLCO$OS,
                         delta = PLCO$Status,
                         marker = train_Age,
                         cause = 1,
                         times = c(365*2,365*3,365 * 5),
                         iid = TRUE)
nomogram_Age$AUC

#T
model_T <- cph(Surv(OS, Status) ~ T,data = PLCO ,x = TRUE, y = T,surv = T)
train_T <- predict(model_T)
nomogram_T <- timeROC(T = PLCO$OS,
                         delta = PLCO$Status,
                         marker = train_T,
                         cause = 1,
                         times = c(365*2,365*3,365 * 5),
                         iid = TRUE)
nomogram_T$AUC
#N
model_N <- cph(Surv(OS, Status) ~ N,data = PLCO ,x = TRUE, y = T,surv = T)
train_N <- predict(model_N)
nomogram_N <- timeROC(T = PLCO$OS,
                      delta = PLCO$Status,
                      marker = train_N,
                      cause = 1,
                      times = c(365*2,365*3,365 * 5),
                      iid = TRUE)
nomogram_N$AUC
#M
model_M <- cph(Surv(OS, Status) ~ M,data = PLCO ,x = TRUE, y = T,surv = T)
train_M <- predict(model_M)
nomogram_M <- timeROC(T = PLCO$OS,
                      delta = PLCO$Status,
                      marker = train_M,
                      cause = 1,
                      times = c(365*2,365*3,365 * 5),
                      iid = TRUE)
nomogram_M$AUC
#TNM
model_TNM <- cph(Surv(OS, Status) ~ T+N+M,data = PLCO ,x = TRUE, y = T,surv = T)
train_TNM <- predict(model_TNM)
nomogram_TNM <- timeROC(T = PLCO$OS,
                      delta = PLCO$Status,
                      marker = train_TNM,
                      cause = 1,
                      times = c(365*2,365*3,365 * 5),
                      iid = TRUE)
nomogram_TNM$AUC
#CRCRS
model_CRCRS<-cph(Surv(OS, Status) ~ CRCRS,data = PLCO ,x = TRUE, y = T,surv = T)
train_CRCRS <- predict(model_CRCRS)
nomogram_CRCRS <- timeROC(T = PLCO$OS,
                        delta = PLCO$Status,
                        marker = train_CRCRS,
                        cause = 1,
                        times = c(365*2,365*3,365 * 5),
                        iid = TRUE)
nomogram_CRCRS$AUC
compare(nomogram_PLCO,nomogram_TNM,adjusted = TRUE)


#TCGA
#nomogram
nomogram_TCGA <- predict(model, newdata = TCGA)
nomogram_TCGA_AUC <- timeROC(T = TCGA$OS,
                                        delta = TCGA$Status,
                                        marker = nomogram_TCGA,
                                        cause = 1,
                                        times =c(365*2,365*3,365 * 5),
                                        iid = TRUE)
nomogram_TCGA_AUC$AUC
nomogram_TCGA_AUC$inference
#Age

nomogram_TCGA_Age <- predict(model_Age, newdata = TCGA)
test_Age <- timeROC(T = TCGA$OS,
                             delta = TCGA$Status,
                             marker =nomogram_TCGA_Age ,
                             cause = 1,
                             times = c(365*2,365*3,365 * 5),
                             iid = TRUE)
test_Age$AUC


#T

nomogram_TCGA_T <- predict(model_T, newdata = TCGA)
test_T <- timeROC(T = TCGA$OS,
                    delta = TCGA$Status,
                    marker =nomogram_TCGA_T ,
                    cause = 1,
                    times = c(365*2,365*3,365 * 5),
                    iid = TRUE)
test_T$AUC

#N
nomogram_TCGA_N <- predict(model_N, newdata = TCGA)
test_N <- timeROC(T = TCGA$OS,
                  delta = TCGA$Status,
                  marker =nomogram_TCGA_N ,
                  cause = 1,
                  times = c(365*2,365*3,365 * 5),
                  iid = TRUE)
test_N$AUC

#M

nomogram_TCGA_M <- predict(model_M, newdata = TCGA)
test_M <- timeROC(T = TCGA$OS,
                  delta = TCGA$Status,
                  marker =nomogram_TCGA_M ,
                  cause = 1,
                  times = c(365*2,365*3,365 * 5),
                  iid = TRUE)
test_M$AUC
#TNM
nomogram_TCGA_TNM <- predict(model_TNM, newdata = TCGA)
test_TNM <- timeROC(T = TCGA$OS,
                  delta = TCGA$Status,
                  marker =nomogram_TCGA_TNM ,
                  cause = 1,
                  times = c(365*2,365*3,365 * 5),
                  iid = TRUE)
test_TNM$AUC
#CRCRS
nomogram_TCGA_CRCRS <- predict(model_CRCRS, newdata = TCGA)
test_CRCRS <- timeROC(T = TCGA$OS,
                    delta = TCGA$Status,
                    marker =nomogram_TCGA_CRCRS ,
                    cause = 1,
                    times = c(365*2,365*3,365 * 5),
                    iid = TRUE)
test_CRCRS$AUC
compare(nomogram_TCGA_AUC,test_TNM,adjusted = TRUE)


#train
par(mar = c(4, 4, 0.5, 0.5))
plot(nomogram_PLCO, time = 365 * 3, col = '#b12b23', lwd = 1, title = 'PLCO')
plot(nomogram_TNM, time = 365 * 3, add = TRUE, col = '#223e9c', lwd = 1)
plot(nomogram_CRCRS, time = 365 * 3, add = TRUE, col = '#f6c619', lwd = 1)
plot(nomogram_Age, time = 365 * 3, add = TRUE,col = '#e1abbc',lwd = 1) 
plot(nomogram_T, time = 365 * 3, add = TRUE, col = '#0eb0c8', lwd = 1)
plot(nomogram_N, time = 365 * 3, add = TRUE, col = '#ea894e', lwd = 1)
plot(nomogram_M, time = 365 * 3, add = TRUE, col = '#0f6657', lwd = 1)
#plot(roc_T, time = 365*5, add = TRUE, col = '#44c1f0', lwd = 1)
#plot(roc_N, time = 365*5, add = TRUE, col = '#e4ce00', lwd = 1)
#plot(roc_M, time = 365*5, add = TRUE, col = '#f18800', lwd = 1)


fig <- legend('bottomright',
              legend = c('Nomogram : 0.84','TNM : 0.80','CRCRS : 0.64','Age : 0.60','T : 0.69','N : 0.71','M : 0.69'),
              col = c('#b12b23','#223e9c','#f6c619','#e1abbc', '#0eb0c8', '#ea894e','#0f6657'),
              lty = 1, lwd = 1, 
              bty='n', adj = c(0, 0.5)
)
#test
par(mar = c(4, 4, 0.5, 0.5))
plot(nomogram_TCGA_AUC, time = 365 * 3, col = '#b12b23', lwd = 1, title = 'TCGA')
plot(test_TNM, time = 365 * 3, add = TRUE, col = '#223e9c', lwd = 1)
plot(test_CRCRS, time = 365 * 3, add = TRUE, col = '#f6c619', lwd = 1)
plot(test_Age, time = 365 * 3, add = TRUE,col = '#e1abbc',lwd = 1)
plot(test_T, time = 365 * 3, add = TRUE, col = '#0eb0c8', lwd = 1)
plot(test_N, time = 365 * 3, add = TRUE, col = '#ea894e', lwd = 1)
plot(test_M, time = 365 * 3, add = TRUE, col = '#0f6657', lwd = 1)
#plot(roc_T, time = 365*5, add = TRUE, col = '#44c1f0', lwd = 1)
#plot(roc_N, time = 365*5, add = TRUE, col = '#e4ce00', lwd = 1)
#plot(roc_M, time = 365*5, add = TRUE, col = '#f18800', lwd = 1)


fig <- legend('bottomright',
              legend = c('Nomogram : 0.77','TNM : 0.71','CRCRS : 0.65','Age : 0.55', 'T : 0.62','N : 0.68','M : 0.61'),
              col = c('#b12b23','#223e9c','#f6c619','#e1abbc', '#0eb0c8', '#ea894e','#0f6657'),
              lty = 1, lwd = 1, 
              bty='n', adj = c(0, 0.5)
)


#DCA曲线
options(unzip ='internal')
devtools::install_github('yikeshu0611/ggDCA')
library(ggDCA)


d_train <- dca(model,model_TNM,model_CRCRS,model_Age,model_T,model_N,model_M,times=c(365*5))
ggplot(d_train)
d_test <- dca(model,model_TNM,model_CRCRS,model_Age,model_T,model_N,model_M,times=c(365*5),new.data=TCGA)
ggplot(d_test)

#RMS
library(rms)
#PLCO
rms_PLCO_nomogram<-calibrate(model,
                   cmethod='KM',
                   method = 'boot',
                   u=1825,
                   m=200,
                   B=500)

rms_PLCO_Age<-calibrate(model_Age,
                             cmethod='KM',
                             method = 'boot',
                             u=1825,
                             m=111,
                             B=1000)

rms_PLCO_T<-calibrate(model_T,
                        cmethod='KM',
                        method = 'boot',
                        u=1825,
                        m=111,
                        B=1000)

rms_PLCO_N<-calibrate(model_N,
                      cmethod='KM',
                      method = 'boot',
                      u=1825,
                      m=111,
                      B=1000)

rms_PLCO_M<-calibrate(model_M,
                      cmethod='KM',
                      method = 'boot',
                      u=1825,
                      m=111,
                      B=1000)
rms_PLCO_TNM<-calibrate(model_TNM,
                      cmethod='KM',
                      method = 'boot',
                      u=1825,
                      m=111,
                      B=1000)
rms_PLCO_CRCRS<-calibrate(model_CRCRS,
                      cmethod='KM',
                      method = 'boot',
                      u=1825,
                      m=111,
                      B=1000)

plot_error<-function(x,y,sd,len=1,col='black'){
  len<-len*0.05
  arrows(x0=x,y0=y,x1=x,y1=y-sd*y,col=col,angle=90,length = len)
  arrows(x0=x,y0=y,x1=x,y1=y+sd*y,col=col,angle=90,length = len)
}
plot(x=1,type='n',
     xlim=c(0,1),
     ylim=c(0,1),
     xaxs='i',
     yaxs='i',
     xlab='Predicted Probability',
     ylab='Observed Probability',
     legend=FALSE,
     subtitles=FALSE,
     cex=1.5,
     cex.axis=1.5,
     cex.lab=1.5)
x1<-rms_PLCO_nomogram[,c('mean.predicted')]
y1<-rms_PLCO_nomogram[,c('KM')]
sd1<-rms_PLCO_nomogram[,c('std.err')]
points(x1,y1,
       type='o',
       pch=0.5,
       col='#b12b23',
       lty=1,
       lwd=2)
plot_error(x1,y1,
           sd=sd1,col='#b12b23')

x2<-rms_PLCO_Age[,c('mean.predicted')]
y2<-rms_PLCO_Age[,c('KM')]
sd2<-rms_PLCO_Age[,c('std.err')]
points(x2,y2,
       type='o',
       pch=0.5,
       col='#e1abbc',
       lty=1,
       lwd=2)
plot_error(x2,y2,
           sd=sd2,col='#e1abbc')

x3<-rms_PLCO_T[,c('mean.predicted')]
y3<-rms_PLCO_T[,c('KM')]
sd3<-rms_PLCO_T[,c('std.err')]
points(x3,y3,
       type='o',
       pch=0.5,
       col='#0eb0c8',
       lty=1,
       lwd=2)
plot_error(x3,y3,
           sd=sd3,col='#0eb0c8')

x4<-rms_PLCO_N[,c('mean.predicted')]
y4<-rms_PLCO_N[,c('KM')]
sd4<-rms_PLCO_N[,c('std.err')]
points(x4,y4,
       type='o',
       pch=0.5,
       col='#ea894e',
       lty=1,
       lwd=2)
plot_error(x4,y4,
           sd=sd4,col='#ea894e')

x5<-rms_PLCO_M[,c('mean.predicted')]
y5<-rms_PLCO_M[,c('KM')]
sd5<-rms_PLCO_M[,c('std.err')]
points(x5,y5,
       type='o',
       pch=0.5,
       col='#0f6657',
       lty=1,
       lwd=2)
plot_error(x5,y5,
           sd=sd5,col='#0f6657')

x6<-rms_PLCO_TNM[,c('mean.predicted')]
y6<-rms_PLCO_TNM[,c('KM')]
sd6<-rms_PLCO_TNM[,c('std.err')]
points(x6,y6,
       type='o',
       pch=0.5,
       col='#223e9c',
       lty=1,
       lwd=2)
plot_error(x6,y6,
           sd=sd6,col='#223e9c')

x7<-rms_PLCO_CRCRS[,c('mean.predicted')]
y7<-rms_PLCO_CRCRS[,c('KM')]
sd7<-rms_PLCO_CRCRS[,c('std.err')]
points(x7,y7,
       type='o',
       pch=0.5,
       col='#f6c619',
       lty=1,
       lwd=2)
plot_error(x7,y7,
           sd=sd7,col='#f6c619')
abline(0,1,lty=3,lwd=2,col='black')


#TCGA
#nomgram
test_nomogram <- c((summary(survfit(model, newdata=TCGA), times=1825)$surv))
cuts_nomogram <- unique(quantile(c(0, 1, test_nomogram), seq(0, 1, length = 5), na.rm = TRUE))
cuts_nomogram
suppressMessages(library(rms))

km_surv_nomogram <- groupkm(test_nomogram,
                     Srv = Surv(TCGA$OS,TCGA$Status),
                     u = 1825,
                     cuts = cuts_nomogram)
km_surv_nomogram
#Age
test_Age <- c((summary(survfit(model_Age, newdata=TCGA), times=1825)$surv))
cuts_Age <- unique(quantile(c(0, 1, test_Age), seq(0, 1, length = 10), na.rm = TRUE))
cuts_Age
suppressMessages(library(rms))

km_surv_Age <- groupkm(test_Age,
                            Srv = Surv(TCGA$OS,TCGA$Status),
                            u = 1825,
                            cuts = cuts_Age)
km_surv_Age



#T
test_T <- c((summary(survfit(model_T, newdata=TCGA), times=1825)$surv))
cuts_T <- unique(quantile(c(0, 1, test_T), seq(0, 1, length = 200), na.rm = TRUE))
cuts_T
suppressMessages(library(rms))

km_surv_t <- groupkm(test_T,
                     Srv = Surv(TCGA$OS,TCGA$Status),
                     u = 1825,
                     cuts = cuts_T)
km_surv_t
#N
test_N <- c((summary(survfit(model_N, newdata=TCGA), times=1825)$surv))
cuts_N <- unique(quantile(c(0, 1, test_N), seq(0, 1, length = 200), na.rm = TRUE))
cuts_N
suppressMessages(library(rms))

km_surv_N <- groupkm(test_N,
                     Srv = Surv(TCGA$OS,TCGA$Status),
                     u = 1825,
                     cuts = cuts_N)
km_surv_N
#M
test_M <- c((summary(survfit(model_M, newdata=TCGA), times=1825)$surv))
cuts_M <- unique(quantile(c(0, 1, test_M), seq(0, 1, length = 200), na.rm = TRUE))
cuts_M
suppressMessages(library(rms))

km_surv_M <- groupkm(test_M,
                     Srv = Surv(TCGA$OS,TCGA$Status),
                     u = 1825,
                     cuts = cuts_M)
km_surv_M
#TNM

test_TNM <- c((summary(survfit(model_TNM, newdata=TCGA), times=1825)$surv))
cuts_TNM <- unique(quantile(c(0, 1, test_TNM), seq(0, 1, length = 5), na.rm = TRUE))
cuts_TNM
suppressMessages(library(rms))

km_surv_TNM <- groupkm(test_TNM,
                     Srv = Surv(TCGA$OS,TCGA$Status),
                     u = 1825,
                     cuts = cuts_TNM)
km_surv_TNM

#CRCRS

test_CRCRS <- c((summary(survfit(model_CRCRS, newdata=TCGA), times=1825)$surv))
cuts_CRCRS <- unique(quantile(c(0, 1, test_CRCRS), seq(0, 1, length = 7), na.rm = TRUE))
cuts_CRCRS
suppressMessages(library(rms))

km_surv_CRCRS <- groupkm(test_CRCRS,
                       Srv = Surv(TCGA$OS,TCGA$Status),
                       u = 1825,
                       cuts = cuts_CRCRS)
km_surv_CRCRS


plot_TCGA_error<-function(x,y,sd,len=1,col='black'){
  len<-len*0.05
  arrows(x0=x,y0=y,x1=x,y1=y-sd*y,col=col,angle=90,length = len)
  arrows(x0=x,y0=y,x1=x,y1=y+sd*y,col=col,angle=90,length = len)
}

plot(x=1,type='n',
     xlim=c(0,1),
     ylim=c(0,1),
     xaxs='i',
     yaxs='i',
     xlab='Predicted Probability',
     ylab='Observed Probability',
     legend=FALSE,
     subtitles=FALSE,
     cex=1.5,
     cex.axis=1.5,
     cex.lab=1.5)
x1<-km_surv_nomogram[,c('x')]
y1<-km_surv_nomogram[,c('KM')]
sd1<-km_surv_nomogram[,c('std.err')]
points(x1,y1,
       type='o',
       pch=0.5,
       col='#b12b23',
       lty=1,
       lwd=2)
plot_TCGA_error(x1,y1,
           sd=sd1,col='#b12b23')
x2<-km_surv_Age[,c('x')]
y2<-km_surv_Age[,c('KM')]
sd2<-km_surv_Age[,c('std.err')]
points(x2,y2,
       type='o',
       pch=0.5,
       col='#e1abbc',
       lty=1,
       lwd=2)
plot_TCGA_error(x2,y2,
           sd=sd2,col='#e1abbc')

x3<-km_surv_t[,c('x')]
y3<-km_surv_t[,c('KM')]
sd3<-km_surv_t[,c('std.err')]
points(x3,y3,
       type='o',
       pch=0.5,
       col='#0eb0c8',
       lty=1,
       lwd=2)
plot_TCGA_error(x3,y3,
           sd=sd3,col='#0eb0c8')

x4<-km_surv_N[,c('x')]
y4<-km_surv_N[,c('KM')]
sd4<-km_surv_N[,c('std.err')]
points(x4,y4,
       type='o',
       pch=0.5,
       col='#ea894e',
       lty=1,
       lwd=2)
plot_TCGA_error(x4,y4,
           sd=sd4,col='#ea894e')

x5<-km_surv_M[,c('x')]
y5<-km_surv_M[,c('KM')]
sd5<-km_surv_M[,c('std.err')]
points(x5,y5,
       type='o',
       pch=0.5,
       col='#0f6657',
       lty=1,
       lwd=2)
plot_TCGA_error(x5,y5,
           sd=sd5,col='#0f6657')

x6<-km_surv_TNM[,c('x')]
y6<-km_surv_TNM[,c('KM')]
sd6<-km_surv_TNM[,c('std.err')]
points(x6,y6,
       type='o',
       pch=0.5,
       col='#223e9c',
       lty=1,
       lwd=2)
plot_TCGA_error(x6,y6,
           sd=sd6,col='#223e9c')

x7<-km_surv_CRCRS[,c('x')]
y7<-km_surv_CRCRS[,c('KM')]
sd7<-km_surv_CRCRS[,c('std.err')]
points(x7,y7,
       type='o',
       pch=0.5,
       col='#f6c619',
       lty=1,
       lwd=2)
plot_TCGA_error(x7,y7,
           sd=sd7,col='#f6c619')
abline(0,1,lty=3,lwd=2,col='black')




# 模型预测的概率



library(survival)
set.seed(123)
trian_PLCO <- c((summary(survfit(model, newdata=PLCO), times=1825)$surv))
head(train_PLCO)
## [1] 0.6958190 0.9347173 0.8681793 0.7350013 0.8694693 0.8213824



library(rms)
#PLCO
rms_PLCO_nomogram_2<-calibrate(model,
                             cmethod='KM',
                             method = 'boot',
                             u=730,
                             m=111,
                             B=500)

rms_PLCO_nomogram_3<-calibrate(model,
                             cmethod='KM',
                             method = 'boot',
                             u=1095,
                             m=111,
                             B=500)
rms_PLCO_nomogram_5<-calibrate(model,
                               cmethod='KM',
                               method = 'boot',
                               u=1825,
                               m=111,
                               B=500)

plot_error<-function(x,y,sd,len=1,col='black'){
  len<-len*0.05
  arrows(x0=x,y0=y,x1=x,y1=y-sd*y,col=col,angle=90,length = len)
  arrows(x0=x,y0=y,x1=x,y1=y+sd*y,col=col,angle=90,length = len)
}
plot(x=1,type='n',
     xlim=c(0,1),
     ylim=c(0,1),
     xaxs='i',
     yaxs='i',
     xlab='Predicted Probability',
     ylab='Observed Probability',
     legend=FALSE,
     subtitles=FALSE,
     cex=1.5,
     cex.axis=1.5,
     cex.lab=1.5)
x1<-rms_PLCO_nomogram_2[,c('mean.predicted')]
y1<-rms_PLCO_nomogram_2[,c('KM')]
sd1<-rms_PLCO_nomogram_2[,c('std.err')]
points(x1,y1,
       type='o',
       pch=0.5,
       col='#b12b23',
       lty=1,
       lwd=2)
plot_error(x1,y1,
           sd=sd1,col='#b12b23')

x2<-rms_PLCO_nomogram_3[,c('mean.predicted')]
y2<-rms_PLCO_nomogram_3[,c('KM')]
sd2<-rms_PLCO_nomogram_3[,c('std.err')]
points(x2,y2,
       type='o',
       pch=0.5,
       col='#e1abbc',
       lty=1,
       lwd=2)
plot_error(x2,y2,
           sd=sd2,col='#e1abbc')

x3<-rms_PLCO_nomogram_5[,c('mean.predicted')]
y3<-rms_PLCO_nomogram_5[,c('KM')]
sd3<-rms_PLCO_nomogram_5[,c('std.err')]
points(x3,y3,
       type='o',
       pch=0.5,
       col='#0eb0c8',
       lty=1,
       lwd=2)
plot_error(x3,y3,
           sd=sd3,col='#0eb0c8')
abline(0,1,lty=3,lwd=2,col='black')



#TCGA
#nomgram
test_nomogram_2 <- c((summary(survfit(model, newdata=TCGA), times=730)$surv))
cuts_nomogram_2 <- unique(quantile(c(0, 1, test_nomogram_2), seq(0, 1, length = 8), na.rm = TRUE))
cuts_nomogram_2
suppressMessages(library(rms))

km_surv_nomogram_2 <- groupkm(test_nomogram_2,
                            Srv = Surv(TCGA$OS,TCGA$Status),
                            u = 730,
                            cuts = cuts_nomogram_2)
km_surv_nomogram_2

test_nomogram_3 <- c((summary(survfit(model, newdata=TCGA), times=1095)$surv))
cuts_nomogram_3 <- unique(quantile(c(0, 1, test_nomogram_3), seq(0, 1, length = 5), na.rm = TRUE))
cuts_nomogram_3
suppressMessages(library(rms))

km_surv_nomogram_3 <- groupkm(test_nomogram_3,
                              Srv = Surv(TCGA$OS,TCGA$Status),
                              u = 1095,
                              cuts = cuts_nomogram_3)
km_surv_nomogram_3

test_nomogram_5 <- c((summary(survfit(model, newdata=TCGA), times=1825)$surv))
cuts_nomogram_5 <- unique(quantile(c(0, 1, test_nomogram_5), seq(0, 1, length = 5), na.rm = TRUE))
cuts_nomogram_5
suppressMessages(library(rms))

km_surv_nomogram_5 <- groupkm(test_nomogram_5,
                              Srv = Surv(TCGA$OS,TCGA$Status),
                              u = 1825,
                              cuts = cuts_nomogram_5)
km_surv_nomogram_5


plot_TCGA_error<-function(x,y,sd,len=1,col='black'){
  len<-len*0.05
  arrows(x0=x,y0=y,x1=x,y1=y-sd*y,col=col,angle=90,length = len)
  arrows(x0=x,y0=y,x1=x,y1=y+sd*y,col=col,angle=90,length = len)
}

plot(x=1,type='n',
     xlim=c(0,1),
     ylim=c(0,1),
     xaxs='i',
     yaxs='i',
     xlab='Predicted Probability',
     ylab='Observed Probability',
     legend=FALSE,
     subtitles=FALSE,
     cex=1.5,
     cex.axis=1.5,
     cex.lab=1.5)
x1<-km_surv_nomogram_2[,c('x')]
y1<-km_surv_nomogram_2[,c('KM')]
sd1<-km_surv_nomogram_2[,c('std.err')]
points(x1,y1,
       type='o',
       pch=0.5,
       col='#b12b23',
       lty=1,
       lwd=2)
plot_TCGA_error(x1,y1,
                sd=sd1,col='#b12b23')
x2<-km_surv_nomogram_3[,c('x')]
y2<-km_surv_nomogram_3[,c('KM')]
sd2<-km_surv_nomogram_3[,c('std.err')]
points(x2,y2,
       type='o',
       pch=0.5,
       col='#e1abbc',
       lty=1,
       lwd=2)
plot_TCGA_error(x2,y2,
                sd=sd2,col='#e1abbc')

x3<-km_surv_nomogram_5[,c('x')]
y3<-km_surv_nomogram_5[,c('KM')]
sd3<-km_surv_nomogram_5[,c('std.err')]
points(x3,y3,
       type='o',
       pch=0.5,
       col='#0eb0c8',
       lty=1,
       lwd=2)
plot_TCGA_error(x3,y3,
                sd=sd3,col='#0eb0c8')
abline(0,1,lty=3,lwd=2,col='black')

#PLCO C-index
options(scipen = 999)
#nomogram
train_nomogram<-predict(model,PLCO)
cindex_train_GIRS_TNM=rcorr.cens(x=-train_nomogram,Surv(PLCO$OS,PLCO$Status))
cindex_train_GIRS_TNM
SD<-rcorr.cens(train_nomogram,Surv(PLCO$OS,PLCO$Status))['S.D.']
cindex_train_GIRS_TNM-1.96*SD/2
cindex_train_GIRS_TNM+1.96*SD/2

#Age
train_Age<-predict(model_Age,PLCO)
cindex_train_Age<-rcorr.cens(x=-train_Age,Surv(PLCO$OS,PLCO$Status))
cindex_train_Age
SD<-rcorr.cens(train_Age,Surv(PLCO$OS,PLCO$Status))['S.D.']
cindex_train_Age-1.96*SD/2
cindex_train_Age+1.96*SD/2
#T
train_T<-predict(model_T,PLCO)
cindex_train_T<-rcorr.cens(x=-train_T,Surv(PLCO$OS,PLCO$Status))
cindex_train_T
SD<-rcorr.cens(train_T,Surv(PLCO$OS,PLCO$Status))['S.D.']
cindex_train_T-1.96*SD/2
cindex_train_T+1.96*SD/2
#N
train_N<-predict(model_N,PLCO)
cindex_train_N<-rcorr.cens(x=-train_N,Surv(PLCO$OS,PLCO$Status))
cindex_train_N
SD<-rcorr.cens(train_N,Surv(PLCO$OS,PLCO$Status))['S.D.']
cindex_train_N-1.96*SD/2
cindex_train_N+1.96*SD/2
#M
train_M<-predict(model_M,PLCO)
cindex_train_M<-rcorr.cens(x=-train_M,Surv(PLCO$OS,PLCO$Status))
cindex_train_M
SD<-rcorr.cens(train_M,Surv(PLCO$OS,PLCO$Status))['S.D.']
cindex_train_M-1.96*SD/2
cindex_train_M+1.96*SD/2
#TNM
train_TNM<-predict(model_TNM,PLCO)
cindex_train_TNM<-rcorr.cens(x=-train_TNM,Surv(PLCO$OS,PLCO$Status))
cindex_train_TNM
SD<-rcorr.cens(train_TNM,Surv(PLCO$OS,PLCO$Status))['S.D.']
cindex_train_TNM-1.96*SD/2
cindex_train_TNM+1.96*SD/2
#CRCRRS
train_CRCRS<-predict(model_CRCRS,PLCO)
cindex_train_CRCRS<-rcorr.cens(x=-train_CRCRS,Surv(PLCO$OS,PLCO$Status))
cindex_train_CRCRS
SD<-rcorr.cens(train_CRCRS,Surv(PLCO$OS,PLCO$Status))['S.D.']
cindex_train_CRCRS-1.96*SD/2
cindex_train_CRCRS+1.96*SD/2


#TCGA C-index
#nomogram
test_nomogram<-predict(model,TCGA)
cindex_test_GIRS_TNM=rcorr.cens(x=-test_nomogram,Surv(TCGA$OS,TCGA$Status))
cindex_test_GIRS_TNM
SD<-rcorr.cens(test_nomogram,Surv(TCGA$OS,TCGA$Status))['S.D.']
cindex_test_GIRS_TNM-1.96*SD/2
cindex_test_GIRS_TNM+1.96*SD/2
#Age
test_Age<-predict(model_Age,TCGA)
cindex_test_GIRS_Age=rcorr.cens(x=-test_Age,Surv(TCGA$OS,TCGA$Status))
cindex_test_GIRS_Age
SD<-rcorr.cens(test_Age,Surv(TCGA$OS,TCGA$Status))['S.D.']
cindex_test_GIRS_Age-1.96*SD/2
cindex_test_GIRS_Age+1.96*SD/2
#T
test_T<-predict(model_T,TCGA)
cindex_test_GIRS_T=rcorr.cens(x=-test_T,Surv(TCGA$OS,TCGA$Status))
cindex_test_GIRS_T
SD<-rcorr.cens(test_T,Surv(TCGA$OS,TCGA$Status))['S.D.']
cindex_test_GIRS_T-1.96*SD/2
cindex_test_GIRS_T+1.96*SD/2
#N
test_N<-predict(model_N,TCGA)
cindex_test_GIRS_N=rcorr.cens(x=-test_N,Surv(TCGA$OS,TCGA$Status))
cindex_test_GIRS_N
SD<-rcorr.cens(test_N,Surv(TCGA$OS,TCGA$Status))['S.D.']
cindex_test_GIRS_N-1.96*SD/2
cindex_test_GIRS_N+1.96*SD/2
#M
test_M<-predict(model_M,TCGA)
cindex_test_GIRS_M=rcorr.cens(x=-test_M,Surv(TCGA$OS,TCGA$Status))
cindex_test_GIRS_M
SD<-rcorr.cens(test_M,Surv(TCGA$OS,TCGA$Status))['S.D.']
cindex_test_GIRS_M-1.96*SD/2
cindex_test_GIRS_M+1.96*SD/2
#TNM
test_TNM<-predict(model_TNM,TCGA)
cindex_test_GIRS_TNM=rcorr.cens(x=-test_TNM,Surv(TCGA$OS,TCGA$Status))
cindex_test_GIRS_TNM
SD<-rcorr.cens(test_TNM,Surv(TCGA$OS,TCGA$Status))['S.D.']
cindex_test_GIRS_TNM-1.96*SD/2
cindex_test_GIRS_TNM+1.96*SD/2
#CRCRS
test_CRCRS<-predict(model_CRCRS,TCGA)
cindex_test_GIRS_CRCRS=rcorr.cens(x=-test_CRCRS,Surv(TCGA$OS,TCGA$Status))
cindex_test_GIRS_CRCRS
SD<-rcorr.cens(test_CRCRS,Surv(TCGA$OS,TCGA$Status))['S.D.']
cindex_test_GIRS_CRCRS-1.96*SD/2
cindex_test_GIRS_CRCRS+1.96*SD/2

#训练集C-index对比
library(compareC)
compareC(timeX=PLCO$OS,PLCO$Status,
         scoreY = -train_nomogram,scoreZ = -train_TNM)
compareC(timeX=PLCO$OS,PLCO$Status,
         scoreY = -train_nomogram,scoreZ = -train_T)
compareC(timeX=PLCO$OS,PLCO$Status,
         scoreY = -train_nomogram,scoreZ = -train_N)
compareC(timeX=PLCO$OS,PLCO$Status,
         scoreY = -train_nomogram,scoreZ = -train_M)
compareC(timeX=PLCO$OS,PLCO$Status,
         scoreY = -train_nomogram,scoreZ = -train_Age)
compareC(timeX=PLCO$OS,PLCO$Status,
         scoreY = -train_nomogram,scoreZ = -train_CRCRS)
#测试集C-index对比
compareC(timeX=TCGA$OS,TCGA$Status,
         scoreY = -test_nomogram,scoreZ = -test_TNM)
compareC(timeX=TCGA$OS,TCGA$Status,
         scoreY = -test_nomogram,scoreZ = -test_T)
compareC(timeX=TCGA$OS,TCGA$Status,
         scoreY = -test_nomogram,scoreZ = -test_N)
compareC(timeX=TCGA$OS,TCGA$Status,
         scoreY = -test_nomogram,scoreZ = -test_M)
compareC(timeX=TCGA$OS,TCGA$Status,
         scoreY = -test_nomogram,scoreZ = -test_Age)
compareC(timeX=TCGA$OS,TCGA$Status,
         scoreY = -test_nomogram,scoreZ = -test_CRCRS)


Cindex<-read.xlsx('cindex.xlsx')
head((Cindex))
## 加载相关的包 
library(ggplot2) 
library(ggthemr) 
library(tidyverse) 
library(dplyr) 
library(ggpubr) 
library(devEMF) 
library(ggsignif) 
names(Cindex) <- c("group","value","SD") 
compared <- list(c("Nomogram","TNM"),c("Nomogram",'CRCRS'),c("Nomogram","T"),c("Nomogram","N"),c("Nomogram","M"),c("Nomogram","Age")) 
Cindex$group<-factor(Cindex$group,levels = c('Nomogram','TNM','CRCRS','T','N',"M","Age"))
ggplot(Cindex, aes(x=group, y=value, fill = group))+ 
  geom_bar(stat = "identity", position = "dodge",color = "black")+ 
  scale_y_continuous(expand = c(0,0),limits = c(0,1))+ ##y从0开始 
  #
  scale_fill_manual(values = c("#ECA8A9","#74AED4", "#D3E2B7","#CFAFD4",'#F7C97E','#eff4fb','#fcf1f0'))+ 
  theme_bw()+theme(legend.position = "none")+ 
  ##  
  geom_errorbar(aes(ymax = value+SD, ymin = value-SD),  
                position = position_dodge(0.01),width = 0.2)+
  ## 添加显著标记 
  geom_signif(annotations = c("***","***","***","***","***","***"),  
              y_position = c(0.75,0.79,0.83,0.87,0.91,0.95),xmin = c(1,1,1,1,1,1), xmax = c(2,3,4,5,6,7))+ 
  labs(y = "C-index", x = "") +
  theme(text = element_text(family = "Arial")) 


Cindex_TCGA<-read.xlsx('TCGA-index.xlsx')

names(Cindex_TCGA) <- c("group","value","SD") 
compared <- list(c("Nomogram","TNM"),c("Nomogram",'CRCRS'),c("Nomogram","T"),c("Nomogram","N"),c("Nomogram","M"),c("Nomogram","Age")) 
Cindex_TCGA$group<-factor(Cindex_TCGA$group,levels = c('Nomogram','TNM','CRCRS','T','N',"M","Age"))
ggplot(Cindex_TCGA, aes(x=group, y=value, fill = group))+ 
  geom_bar(stat = "identity", position = "dodge",color = "black")+ 
  scale_y_continuous(expand = c(0,0),limits = c(0,1.2))+ ##y从0开始 
  #
  scale_fill_manual(values = c("#ECA8A9","#74AED4", "#D3E2B7","#CFAFD4",'#F7C97E','#eff4fb','#fcf1f0'))+ 
  theme_bw()+theme(legend.position = "none")+ 
  ##  
  geom_errorbar(aes(ymax = value+SD, ymin = value-SD),  
                position = position_dodge(0.01),width = 0.2)+
  ## 
  geom_signif(annotations = c("***","***","***","***","***","***"),  
              y_position = c(0.83,0.88,0.93,0.98,1.03,1.08),xmin = c(1,1,1,1,1,1), xmax = c(2,3,4,5,6,7))+ 
  labs(y = "C-index", x = "") +
  theme(text = element_text(family = "Arial")) 






