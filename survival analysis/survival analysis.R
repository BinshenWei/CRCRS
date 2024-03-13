library(survival)
library(plyr)
library(openxlsx)
# Set working directory
setwd("/Users/weibinshen/Desktop/未命名文件夹/PLCO预测")
TCGA<-read.xlsx('TCGA.xlsx')
PLCO<-read.xlsx('PLCO.xlsx')
categorical_vars <- c("Gender", "M","T","N",'Age')
TCGA[, categorical_vars] <- lapply(TCGA[, categorical_vars], as.factor)
PLCO[, categorical_vars] <- lapply(PLCO[, categorical_vars], as.factor)
#Val[, categorical_vars] <- lapply(Val[, categorical_vars], as.factor)
#Codes for cut-point of GIRS determination
library(survminer)
library(survival)
cut_point_train <- surv_cutpoint(data =PLCO , time = 'OS', event = 'Status', variables = 'CRCRS',minprop = 0.000001, progressbar = TRUE)
res.cat <- surv_categorize(cut_point_train)
cut_point_train$cutpoint
plot(cut_point_train, palette ='npg', legend = "",bins =15, yscale='sqrt')
cat_point_train <- surv_categorize(cut_point_train)
fit <- survfit(Surv(OS, Status) ~CRCRS, data =cat_point_train )
OS_trainset <- ggsurvplot(fit, palette = "jco",
                          risk.table = TRUE, pval = TRUE, conf.int = TRUE,
                          risk.table.col = "strata",
                          linetype = "strata",
                          legend.title = "Risk levels", legend.labs = c("High risk", "Low risk"),
                          pval.size = 6, fontsize = 4.5, font.legend = 14,
                          surv.median.line = "hv",
                          legend = c(0.8, 0.25), pval.coord = c(0.15, 0.08),
                          risk.table.y.text.col = TRUE,
                          risk.table.y.text = FALSE,
                          tables.height = 0.2,
                          tables.theme = theme_survminer(),
                          xlab = "Time in days", risk.table.height = 0.25,
                          ylab = "Overall survival Probability",
                          xlim = c(0, 2920),
                          break.time.by =365,# Change the xlim to restrict the time range to 5 years (1825 days)
                          data =cat_point_train )

OS_trainset

TCGA$CRCRS = ifelse(TCGA$CRCRS>cut_point_train$cutpoint$cutpoint,'high','low')
fit <- survfit(Surv(OS, Status) ~CRCRS, data =TCGA )
OS_val <- ggsurvplot(fit,palette = "jco",
                     risk.table =TRUE,pval =TRUE,conf.int=T,
                     risk.table.col = "strata", # Change risk table color by groups
                     linetype = "strata", # Change line type by groups
                     legend.title = "Risk levels",legend.labs=c("High risk","Low risk"),
                     pval.size=6,fontsize=4.5,font.legend=14,
                     surv.median.line = "hv",
                     legend=c(0.85, 0.8),pval.coord=c(0.15, 0.08),
                     risk.table.y.text.col = T,
                     risk.table.y.text = FALSE,
                     tables.height = 0.2,
                     tables.theme = theme_survminer(),#theme_cleantable(),
                     xlab ="Time in days",risk.table.height=.25,
                     ylab="Overall survival Probality",
                     xlim = c(0, 2920),
                     break.time.by = 365,
                     data = TCGA)
OS_val





TCGA<-read.xlsx('TCGA.xlsx')
PLCO<-read.xlsx('PLCO.xlsx')
categorical_vars <- c("Gender", "M","T","N",'Age')
TCGA[, categorical_vars] <- lapply(TCGA[, categorical_vars], as.factor)
PLCO[, categorical_vars] <- lapply(PLCO[, categorical_vars], as.factor)
#Val[, categorical_vars] <- lapply(Val[, categorical_vars], as.factor)

# Perform one-way Cox regression analysis on PLCO data frame
train_surv <- Surv(time = PLCO$OS, event = PLCO$Status)
vars <- colnames(PLCO)[c(4:9)]
UniCox <- function(x) {
  FML <- as.formula(paste0('train_surv ~', x))
  fit <- coxph(FML, data = Val)
  fitSum <- summary(fit)
  HR <- round(fitSum$coefficients[,2], 4)
  PValue <- round(fitSum$coefficients[,5], 5)
  CI <- paste0(round(fitSum$conf.int[,3:4], 5))
  Unicox <- data.frame('Characteristics' = x,
                       'Hazard Ratio' = HR,
                       'CI95' = CI,
                       'P Value' = PValue)
  return(Unicox)
}

# Perform one-way Cox regression analysis and output results
results_train <- lapply(vars, UniCox)
results_train 
results_train_df <- do.call(rbind, results_train)
print(results_train_df)

# Perform a multifactor Cox regression analysis on the train data frame.
str(Val)
options(scipen = 999)
multCox <- coxph(Surv(OS, Status) ~T+N + M+CRCRS+Age+Gender,data = Val)
summary(multCox)
multCoxSum <- summary(multCox)
multCoxHR <- signif(as.matrix(multCoxSum$coefficients)[,2],4)
multCoxPValue <- signif(as.matrix(multCoxSum$coefficients)[,5],2)
multCoxCI <- signif(as.matrix(multCoxSum$conf.int[,3:4],2))
multi_res <- data.frame(
  'Hazard Ratio' = multCoxHR,
  'CI95.lower' = multCoxCI[,1],
  'CI95.upper' = multCoxCI[,2],
  'P Value' = multCoxPValue,
  stringsAsFactors = FALSE
)

print(multi_res)



# Performing one-way Cox regression analysis on the TCGA data frame
test_surv <- Surv(time = TCGA$OS, event = TCGA$Status)
vars <- colnames(TCGA)[c(4:9)]
UniCox <- function(x) {
  FML <- as.formula(paste0('test_surv ~', x))
  fit <- coxph(FML, data = TCGA)
  fitSum <- summary(fit)
  HR <- round(fitSum$coefficients[,2], 4)
  PValue <- round(fitSum$coefficients[,5], 5)
  CI <- paste0(round(fitSum$conf.int[,3:4], 5))
  Unicox <- data.frame('Characteristics' = x,
                       'Hazard Ratio' = HR,
                       'CI95' = CI,
                       'P Value' = PValue)
  return(Unicox)
}
results_test <- lapply(vars, UniCox)
results_test 
results_test_df <- do.call(rbind, results_test)
print(results_test_df)

# Perform multifactor Cox regression analysis on TCGA data frame
str(TCGA)
options(scipen = 999)
multCox <- coxph(Surv(OS, Status) ~T+N + M+CRCRS+Age,data = TCGA)
summary(multCox)
multCoxSum <- summary(multCox)
multCoxHR <- signif(as.matrix(multCoxSum$coefficients)[,2],4)
multCoxPValue <- signif(as.matrix(multCoxSum$coefficients)[,5],2)
multCoxCI <- signif(as.matrix(multCoxSum$conf.int[,3:4],2))
multi_res <- data.frame(
  'Hazard Ratio' = multCoxHR,
  'CI95.lower' = multCoxCI[,1],
  'CI95.upper' = multCoxCI[,2],
  'P Value' = multCoxPValue,
  stringsAsFactors = FALSE
)

print(multi_res)


#risk forest diagram
library(survminer)
library(survival)
TCGA<-read.xlsx('TCGA.xlsx')
PLCO<-read.xlsx('PLCO.xlsx')
categorical_vars <- c("Gender", "M","T","N",'Age')
TCGA[, categorical_vars] <- lapply(TCGA[, categorical_vars], as.factor)
PLCO[, categorical_vars] <- lapply(PLCO[, categorical_vars], as.factor)
multCox <- coxph(Surv(OS, Status) ~ T + N + M + CRCRS+Age+Gender , data = PLCO)

pdf("risk forest diagram.pdf",12,8)
ggforest(multCox, data = PLCO, fontsize = 1)
dev.off()
#TCGA
multTCGA_Cox <- coxph(Surv(OS, Status) ~ T + N + M + CRCRS+Age+Gender , data = TCGA)
ggforest(multTCGA_Cox, data = TCGA, fontsize = 1)


