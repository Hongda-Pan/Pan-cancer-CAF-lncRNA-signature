setwd("D:\\BaiduNetdiskWorkspace\\05CAF\\9.8cor.GSVA.bestseparation")
library(survival)
library(survminer)
#---------------------------------------------ACC------------------------------------
svdata <- read.table("3score.time.ACC.txt",sep="\t",header = T,row.names = 1)
dim(svdata)
pct = 0.1
expr <- svdata[,3:ncol(svdata)]
#expr <- log2(expr + 1)
colnames(expr) <- make.names(colnames(expr))
svdata <- cbind.data.frame(svdata[,1:2],expr)
res.cut <- surv_cutpoint(svdata, time = "futime", 
                         event = "fustat", 
                         variables = names(svdata)[3:ncol(svdata)], 
                         minprop = 0.3) #Ĭ??????sample???ܵ???30%
res.cat <- surv_categorize(res.cut)
res.cat$futime=res.cat$futime/30   #####ʱ?�?????????
my.surv <- Surv(res.cat$futime, res.cat$fustat)
pl <- pl2 <- list() # ??ʼ?????????ߺ??¼??????????ߵ??б?
ptable <- data.frame(matrix(NA,ncol(svdata)-2,5))
colnames(ptable) <- c("ID", "pvalue", "HR", "CIlow", "CIup")
n = 0
for (i in colnames(res.cat)[3:ncol(svdata)]) {
  group <- res.cat[,i] 
  survival_dat <- data.frame(group = group)
  fit <- survfit(my.surv ~ group)
  group <- factor(group, levels = c("low", "high"))
  data.survdiff <- survdiff(my.surv ~ group)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  CancerType="ACC"
  n = n + 1
  ptable$ID[n] <- i
  ptable$pvalue[n] <- p.val
  ptable$HR[n] <- HR
  ptable$CIlow[n] <- low95
  ptable$CIup[n] <- up95
  ptable$CancerType[n] <- CancerType
    HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
  CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
  svsort <- svdata[order(svdata[,i]),]
  pl[[i]]<-ggsurvplot(fit, data = survival_dat ,
                      #ggtheme = theme_bw(), #??Ҫ??????????????
                      conf.int = T, #???????????????????????????Ͱ?F?ĳ?T
                      #conf.int.style = "step",#?????????????ͣ????ɸ?Ϊribbon
                      censor = F, #????ʾ?۲?ֵ???ڵ?λ??
                      palette = c("#e64b35","#1b90be"), #?ߵ???ɫ??Ӧ?ߡ???
                      xlab="Time(Months)",
                      legend = c(0.8, 0.9),
                      legend.title = i,#??????д??ͼ????Ŀ??λ??
                      font.legend = 11,#ͼ??????????С
                      #font.title = 12,font.x = 10,font.y = 10,#??????????????С
                      #??ͼ???ϱ????ߵͷֽ????ı???��????????sample??��
                      legend.labs=c(paste0(">",round(svsort[fit$n[2],i],2),"(",fit$n[1],")"),
                                    paste0("<",round(svsort[fit$n[2],i],2),"(",fit$n[2],")")),
                      #????pvalue??HR??95% CI
                      #̫С??p value??Ϊp < 0.001
                      pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                                 paste("p = ",round(p.val,3), sep = "")),
                                   HR, CI, sep = "\n"))
ggsave(paste0(i,".ACC.pdf"),width = 4,height = 4)}
length(pl)
F02=t(ptable)#write.csv(ptable, "pvalue_output_BLCA.csv", quote = F, row.names = F)


#---------------------------------------------BLCA------------------------------------
svdata <- read.table("3score.time.BLCA.txt",sep="\t",header = T,row.names = 1)
dim(svdata)
pct = 0.1
expr <- svdata[,3:ncol(svdata)]
#expr <- log2(expr + 1)
colnames(expr) <- make.names(colnames(expr))
svdata <- cbind.data.frame(svdata[,1:2],expr)
res.cut <- surv_cutpoint(svdata, time = "futime", 
                         event = "fustat", 
                         variables = names(svdata)[3:ncol(svdata)], 
                         minprop = 0.3) #Ĭ??????sample???ܵ???30%
res.cat <- surv_categorize(res.cut)
res.cat$futime=res.cat$futime/30   #####ʱ?�?????????
my.surv <- Surv(res.cat$futime, res.cat$fustat)
pl <- pl2 <- list() # ??ʼ?????????ߺ??¼??????????ߵ??б?
ptable <- data.frame(matrix(NA,ncol(svdata)-2,5))
colnames(ptable) <- c("ID", "pvalue", "HR", "CIlow", "CIup")
n = 0
for (i in colnames(res.cat)[3:ncol(svdata)]) {
  group <- res.cat[,i] 
  survival_dat <- data.frame(group = group)
  fit <- survfit(my.surv ~ group)
  group <- factor(group, levels = c("low", "high"))
  data.survdiff <- survdiff(my.surv ~ group)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  CancerType="BLCA"
  n = n + 1
  ptable$ID[n] <- i
  ptable$pvalue[n] <- p.val
  ptable$HR[n] <- HR
  ptable$CIlow[n] <- low95
  ptable$CIup[n] <- up95
  ptable$CancerType[n] <- CancerType
  HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
  CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
  svsort <- svdata[order(svdata[,i]),]
  pl[[i]]<-ggsurvplot(fit, data = survival_dat ,
                      #ggtheme = theme_bw(), #??Ҫ??????????????
                      conf.int = T, #???????????????????????????Ͱ?F?ĳ?T
                      #conf.int.style = "step",#?????????????ͣ????ɸ?Ϊribbon
                      censor = F, #????ʾ?۲?ֵ???ڵ?λ??
                      palette = c("#e64b35","#1b90be"), #?ߵ???ɫ??Ӧ?ߡ???
                      xlab="Time(Months)",
                      legend = c(0.8, 0.9),
                      legend.title = i,#??????д??ͼ????Ŀ??λ??
                      font.legend = 11,#ͼ??????????С
                      #font.title = 12,font.x = 10,font.y = 10,#??????????????С
                      #??ͼ???ϱ????ߵͷֽ????ı???��????????sample??��
                      legend.labs=c(paste0(">",round(svsort[fit$n[2],i],2),"(",fit$n[1],")"),
                                    paste0("<",round(svsort[fit$n[2],i],2),"(",fit$n[2],")")),
                      #????pvalue??HR??95% CI
                      #̫С??p value??Ϊp < 0.001
                      pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                                 paste("p = ",round(p.val,3), sep = "")),
                                   HR, CI, sep = "\n"))
  ggsave(paste0(i,".BLCA.pdf"),width = 4,height = 4)}
length(pl)
F03=t(ptable)#write.csv(ptable, "pvalue_output_BLCA.csv", quote = F, row.names = F)


#---------------------------------------------BRCA------------------------------------
svdata <- read.table("3score.time.BRCA.txt",sep="\t",header = T,row.names = 1)
dim(svdata)
pct = 0.1
expr <- svdata[,3:ncol(svdata)]
#expr <- log2(expr + 1)
colnames(expr) <- make.names(colnames(expr))
svdata <- cbind.data.frame(svdata[,1:2],expr)
res.cut <- surv_cutpoint(svdata, time = "futime", 
                         event = "fustat", 
                         variables = names(svdata)[3:ncol(svdata)], 
                         minprop = 0.3) #Ĭ??????sample???ܵ???30%
res.cat <- surv_categorize(res.cut)
res.cat$futime=res.cat$futime/30   #####ʱ?�?????????
my.surv <- Surv(res.cat$futime, res.cat$fustat)
pl <- pl2 <- list() # ??ʼ?????????ߺ??¼??????????ߵ??б?
ptable <- data.frame(matrix(NA,ncol(svdata)-2,5))
colnames(ptable) <- c("ID", "pvalue", "HR", "CIlow", "CIup")
n = 0
for (i in colnames(res.cat)[3:ncol(svdata)]) {
  group <- res.cat[,i] 
  survival_dat <- data.frame(group = group)
  fit <- survfit(my.surv ~ group)
  group <- factor(group, levels = c("low", "high"))
  data.survdiff <- survdiff(my.surv ~ group)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  CancerType="BRCA"
  n = n + 1
  ptable$ID[n] <- i
  ptable$pvalue[n] <- p.val
  ptable$HR[n] <- HR
  ptable$CIlow[n] <- low95
  ptable$CIup[n] <- up95
  ptable$CancerType[n] <- CancerType
  HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
  CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
  svsort <- svdata[order(svdata[,i]),]
  pl[[i]]<-ggsurvplot(fit, data = survival_dat ,
                      #ggtheme = theme_bw(), #??Ҫ??????????????
                      conf.int = T, #???????????????????????????Ͱ?F?ĳ?T
                      #conf.int.style = "step",#?????????????ͣ????ɸ?Ϊribbon
                      censor = F, #????ʾ?۲?ֵ???ڵ?λ??
                      palette = c("#e64b35","#1b90be"), #?ߵ???ɫ??Ӧ?ߡ???
                      xlab="Time(Months)",
                      legend = c(0.8, 0.9),
                      legend.title = i,#??????д??ͼ????Ŀ??λ??
                      font.legend = 11,#ͼ??????????С
                      #font.title = 12,font.x = 10,font.y = 10,#??????????????С
                      #??ͼ???ϱ????ߵͷֽ????ı???��????????sample??��
                      legend.labs=c(paste0(">",round(svsort[fit$n[2],i],2),"(",fit$n[1],")"),
                                    paste0("<",round(svsort[fit$n[2],i],2),"(",fit$n[2],")")),
                      #????pvalue??HR??95% CI
                      #̫С??p value??Ϊp < 0.001
                      pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                                 paste("p = ",round(p.val,3), sep = "")),
                                   HR, CI, sep = "\n"))
  ggsave(paste0(i,".BRCA.pdf"),width = 4,height = 4)}
length(pl)
F04=t(ptable)#write.csv(ptable, "pvalue_output_BLCA.csv", quote = F, row.names = F)


#---------------------------------------------CESC------------------------------------
svdata <- read.table("3score.time.CESC.txt",sep="\t",header = T,row.names = 1)
dim(svdata)
pct = 0.1
expr <- svdata[,3:ncol(svdata)]
#expr <- log2(expr + 1)
colnames(expr) <- make.names(colnames(expr))
svdata <- cbind.data.frame(svdata[,1:2],expr)
res.cut <- surv_cutpoint(svdata, time = "futime", 
                         event = "fustat", 
                         variables = names(svdata)[3:ncol(svdata)], 
                         minprop = 0.3) #Ĭ??????sample???ܵ???30%
res.cat <- surv_categorize(res.cut)
res.cat$futime=res.cat$futime/30   #####ʱ?�?????????
my.surv <- Surv(res.cat$futime, res.cat$fustat)
pl <- pl2 <- list() # ??ʼ?????????ߺ??¼??????????ߵ??б?
ptable <- data.frame(matrix(NA,ncol(svdata)-2,5))
colnames(ptable) <- c("ID", "pvalue", "HR", "CIlow", "CIup")
n = 0
for (i in colnames(res.cat)[3:ncol(svdata)]) {
  group <- res.cat[,i] 
  survival_dat <- data.frame(group = group)
  fit <- survfit(my.surv ~ group)
  group <- factor(group, levels = c("low", "high"))
  data.survdiff <- survdiff(my.surv ~ group)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  CancerType="CESC"
  n = n + 1
  ptable$ID[n] <- i
  ptable$pvalue[n] <- p.val
  ptable$HR[n] <- HR
  ptable$CIlow[n] <- low95
  ptable$CIup[n] <- up95
  ptable$CancerType[n] <- CancerType
  HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
  CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
  svsort <- svdata[order(svdata[,i]),]
  pl[[i]]<-ggsurvplot(fit, data = survival_dat ,
                      #ggtheme = theme_bw(), #??Ҫ??????????????
                      conf.int = T, #???????????????????????????Ͱ?F?ĳ?T
                      #conf.int.style = "step",#?????????????ͣ????ɸ?Ϊribbon
                      censor = F, #????ʾ?۲?ֵ???ڵ?λ??
                      palette = c("#e64b35","#1b90be"), #?ߵ???ɫ??Ӧ?ߡ???
                      xlab="Time(Months)",
                      legend = c(0.8, 0.9),
                      legend.title = i,#??????д??ͼ????Ŀ??λ??
                      font.legend = 11,#ͼ??????????С
                      #font.title = 12,font.x = 10,font.y = 10,#??????????????С
                      #??ͼ???ϱ????ߵͷֽ????ı???��????????sample??��
                      legend.labs=c(paste0(">",round(svsort[fit$n[2],i],2),"(",fit$n[1],")"),
                                    paste0("<",round(svsort[fit$n[2],i],2),"(",fit$n[2],")")),
                      #????pvalue??HR??95% CI
                      #̫С??p value??Ϊp < 0.001
                      pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                                 paste("p = ",round(p.val,3), sep = "")),
                                   HR, CI, sep = "\n"))
  ggsave(paste0(i,".CESC.pdf"),width = 4,height = 4)}
length(pl)
F05=t(ptable)#write.csv(ptable, "pvalue_output_BLCA.csv", quote = F, row.names = F)


#---------------------------------------------CHOL------------------------------------
svdata <- read.table("3score.time.CHOL.txt",sep="\t",header = T,row.names = 1)
dim(svdata)
pct = 0.1
expr <- svdata[,3:ncol(svdata)]
#expr <- log2(expr + 1)
colnames(expr) <- make.names(colnames(expr))
svdata <- cbind.data.frame(svdata[,1:2],expr)
res.cut <- surv_cutpoint(svdata, time = "futime", 
                         event = "fustat", 
                         variables = names(svdata)[3:ncol(svdata)], 
                         minprop = 0.3) #Ĭ??????sample???ܵ???30%
res.cat <- surv_categorize(res.cut)
res.cat$futime=res.cat$futime/30   #####ʱ?�?????????
my.surv <- Surv(res.cat$futime, res.cat$fustat)
pl <- pl2 <- list() # ??ʼ?????????ߺ??¼??????????ߵ??б?
ptable <- data.frame(matrix(NA,ncol(svdata)-2,5))
colnames(ptable) <- c("ID", "pvalue", "HR", "CIlow", "CIup")
n = 0
for (i in colnames(res.cat)[3:ncol(svdata)]) {
  group <- res.cat[,i] 
  survival_dat <- data.frame(group = group)
  fit <- survfit(my.surv ~ group)
  group <- factor(group, levels = c("low", "high"))
  data.survdiff <- survdiff(my.surv ~ group)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  CancerType="CHOL"
  n = n + 1
  ptable$ID[n] <- i
  ptable$pvalue[n] <- p.val
  ptable$HR[n] <- HR
  ptable$CIlow[n] <- low95
  ptable$CIup[n] <- up95
  ptable$CancerType[n] <- CancerType
  HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
  CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
  svsort <- svdata[order(svdata[,i]),]
  pl[[i]]<-ggsurvplot(fit, data = survival_dat ,
                      #ggtheme = theme_bw(), #??Ҫ??????????????
                      conf.int = T, #???????????????????????????Ͱ?F?ĳ?T
                      #conf.int.style = "step",#?????????????ͣ????ɸ?Ϊribbon
                      censor = F, #????ʾ?۲?ֵ???ڵ?λ??
                      palette = c("#e64b35","#1b90be"), #?ߵ???ɫ??Ӧ?ߡ???
                      xlab="Time(Months)",
                      legend = c(0.8, 0.9),
                      legend.title = i,#??????д??ͼ????Ŀ??λ??
                      font.legend = 11,#ͼ??????????С
                      #font.title = 12,font.x = 10,font.y = 10,#??????????????С
                      #??ͼ???ϱ????ߵͷֽ????ı???��????????sample??��
                      legend.labs=c(paste0(">",round(svsort[fit$n[2],i],2),"(",fit$n[1],")"),
                                    paste0("<",round(svsort[fit$n[2],i],2),"(",fit$n[2],")")),
                      #????pvalue??HR??95% CI
                      #̫С??p value??Ϊp < 0.001
                      pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                                 paste("p = ",round(p.val,3), sep = "")),
                                   HR, CI, sep = "\n"))
  ggsave(paste0(i,".CHOL.pdf"),width = 4,height = 4)}
length(pl)
F06=t(ptable)#write.csv(ptable, "pvalue_output_BLCA.csv", quote = F, row.names = F)


#---------------------------------------------COAD------------------------------------
svdata <- read.table("3score.time.COAD.txt",sep="\t",header = T,row.names = 1)
dim(svdata)
pct = 0.1
expr <- svdata[,3:ncol(svdata)]
#expr <- log2(expr + 1)
colnames(expr) <- make.names(colnames(expr))
svdata <- cbind.data.frame(svdata[,1:2],expr)
res.cut <- surv_cutpoint(svdata, time = "futime", 
                         event = "fustat", 
                         variables = names(svdata)[3:ncol(svdata)], 
                         minprop = 0.3) #Ĭ??????sample???ܵ???30%
res.cat <- surv_categorize(res.cut)
res.cat$futime=res.cat$futime/30   #####ʱ?�?????????
my.surv <- Surv(res.cat$futime, res.cat$fustat)
pl <- pl2 <- list() # ??ʼ?????????ߺ??¼??????????ߵ??б?
ptable <- data.frame(matrix(NA,ncol(svdata)-2,5))
colnames(ptable) <- c("ID", "pvalue", "HR", "CIlow", "CIup")
n = 0
for (i in colnames(res.cat)[3:ncol(svdata)]) {
  group <- res.cat[,i] 
  survival_dat <- data.frame(group = group)
  fit <- survfit(my.surv ~ group)
  group <- factor(group, levels = c("low", "high"))
  data.survdiff <- survdiff(my.surv ~ group)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  CancerType="COAD"
  n = n + 1
  ptable$ID[n] <- i
  ptable$pvalue[n] <- p.val
  ptable$HR[n] <- HR
  ptable$CIlow[n] <- low95
  ptable$CIup[n] <- up95
  ptable$CancerType[n] <- CancerType
  HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
  CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
  svsort <- svdata[order(svdata[,i]),]
  pl[[i]]<-ggsurvplot(fit, data = survival_dat ,
                      #ggtheme = theme_bw(), #??Ҫ??????????????
                      conf.int = T, #???????????????????????????Ͱ?F?ĳ?T
                      #conf.int.style = "step",#?????????????ͣ????ɸ?Ϊribbon
                      censor = F, #????ʾ?۲?ֵ???ڵ?λ??
                      palette = c("#e64b35","#1b90be"), #?ߵ???ɫ??Ӧ?ߡ???
                      xlab="Time(Months)",
                      legend = c(0.8, 0.9),
                      legend.title = i,#??????д??ͼ????Ŀ??λ??
                      font.legend = 11,#ͼ??????????С
                      #font.title = 12,font.x = 10,font.y = 10,#??????????????С
                      #??ͼ???ϱ????ߵͷֽ????ı???��????????sample??��
                      legend.labs=c(paste0(">",round(svsort[fit$n[2],i],2),"(",fit$n[1],")"),
                                    paste0("<",round(svsort[fit$n[2],i],2),"(",fit$n[2],")")),
                      #????pvalue??HR??95% CI
                      #̫С??p value??Ϊp < 0.001
                      pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                                 paste("p = ",round(p.val,3), sep = "")),
                                   HR, CI, sep = "\n"))
  ggsave(paste0(i,".COAD.pdf"),width = 4,height = 4)}
length(pl)
F07=t(ptable)#write.csv(ptable, "pvalue_output_BLCA.csv", quote = F, row.names = F)


#---------------------------------------------DLBC------------------------------------
svdata <- read.table("3score.time.DLBC.txt",sep="\t",header = T,row.names = 1)
dim(svdata)
pct = 0.1
expr <- svdata[,3:ncol(svdata)]
#expr <- log2(expr + 1)
colnames(expr) <- make.names(colnames(expr))
svdata <- cbind.data.frame(svdata[,1:2],expr)
res.cut <- surv_cutpoint(svdata, time = "futime", 
                         event = "fustat", 
                         variables = names(svdata)[3:ncol(svdata)], 
                         minprop = 0.3) #Ĭ??????sample???ܵ???30%
res.cat <- surv_categorize(res.cut)
res.cat$futime=res.cat$futime/30   #####ʱ?�?????????
my.surv <- Surv(res.cat$futime, res.cat$fustat)
pl <- pl2 <- list() # ??ʼ?????????ߺ??¼??????????ߵ??б?
ptable <- data.frame(matrix(NA,ncol(svdata)-2,5))
colnames(ptable) <- c("ID", "pvalue", "HR", "CIlow", "CIup")
n = 0
for (i in colnames(res.cat)[3:ncol(svdata)]) {
  group <- res.cat[,i] 
  survival_dat <- data.frame(group = group)
  fit <- survfit(my.surv ~ group)
  group <- factor(group, levels = c("low", "high"))
  data.survdiff <- survdiff(my.surv ~ group)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  CancerType="DLBC"
  n = n + 1
  ptable$ID[n] <- i
  ptable$pvalue[n] <- p.val
  ptable$HR[n] <- HR
  ptable$CIlow[n] <- low95
  ptable$CIup[n] <- up95
  ptable$CancerType[n] <- CancerType
  HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
  CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
  svsort <- svdata[order(svdata[,i]),]
  pl[[i]]<-ggsurvplot(fit, data = survival_dat ,
                      #ggtheme = theme_bw(), #??Ҫ??????????????
                      conf.int = T, #???????????????????????????Ͱ?F?ĳ?T
                      #conf.int.style = "step",#?????????????ͣ????ɸ?Ϊribbon
                      censor = F, #????ʾ?۲?ֵ???ڵ?λ??
                      palette = c("#e64b35","#1b90be"), #?ߵ???ɫ??Ӧ?ߡ???
                      xlab="Time(Months)",
                      legend = c(0.8, 0.9),
                      legend.title = i,#??????д??ͼ????Ŀ??λ??
                      font.legend = 11,#ͼ??????????С
                      #font.title = 12,font.x = 10,font.y = 10,#??????????????С
                      #??ͼ???ϱ????ߵͷֽ????ı???��????????sample??��
                      legend.labs=c(paste0(">",round(svsort[fit$n[2],i],2),"(",fit$n[1],")"),
                                    paste0("<",round(svsort[fit$n[2],i],2),"(",fit$n[2],")")),
                      #????pvalue??HR??95% CI
                      #̫С??p value??Ϊp < 0.001
                      pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                                 paste("p = ",round(p.val,3), sep = "")),
                                   HR, CI, sep = "\n"))
  ggsave(paste0(i,".DLBC.pdf"),width = 4,height = 4)}
length(pl)
F08=t(ptable)#write.csv(ptable, "pvalue_output_BLCA.csv", quote = F, row.names = F)


#---------------------------------------------ESCA------------------------------------
svdata <- read.table("3score.time.ESCA.txt",sep="\t",header = T,row.names = 1)
dim(svdata)
pct = 0.1
expr <- svdata[,3:ncol(svdata)]
#expr <- log2(expr + 1)
colnames(expr) <- make.names(colnames(expr))
svdata <- cbind.data.frame(svdata[,1:2],expr)
res.cut <- surv_cutpoint(svdata, time = "futime", 
                         event = "fustat", 
                         variables = names(svdata)[3:ncol(svdata)], 
                         minprop = 0.3) #Ĭ??????sample???ܵ???30%
res.cat <- surv_categorize(res.cut)
res.cat$futime=res.cat$futime/30   #####ʱ?�?????????
my.surv <- Surv(res.cat$futime, res.cat$fustat)
pl <- pl2 <- list() # ??ʼ?????????ߺ??¼??????????ߵ??б?
ptable <- data.frame(matrix(NA,ncol(svdata)-2,5))
colnames(ptable) <- c("ID", "pvalue", "HR", "CIlow", "CIup")
n = 0
for (i in colnames(res.cat)[3:ncol(svdata)]) {
  group <- res.cat[,i] 
  survival_dat <- data.frame(group = group)
  fit <- survfit(my.surv ~ group)
  group <- factor(group, levels = c("low", "high"))
  data.survdiff <- survdiff(my.surv ~ group)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  CancerType="ESCA"
  n = n + 1
  ptable$ID[n] <- i
  ptable$pvalue[n] <- p.val
  ptable$HR[n] <- HR
  ptable$CIlow[n] <- low95
  ptable$CIup[n] <- up95
  ptable$CancerType[n] <- CancerType
  HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
  CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
  svsort <- svdata[order(svdata[,i]),]
  pl[[i]]<-ggsurvplot(fit, data = survival_dat ,
                      #ggtheme = theme_bw(), #??Ҫ??????????????
                      conf.int = T, #???????????????????????????Ͱ?F?ĳ?T
                      #conf.int.style = "step",#?????????????ͣ????ɸ?Ϊribbon
                      censor = F, #????ʾ?۲?ֵ???ڵ?λ??
                      palette = c("#e64b35","#1b90be"), #?ߵ???ɫ??Ӧ?ߡ???
                      xlab="Time(Months)",
                      legend = c(0.8, 0.9),
                      legend.title = i,#??????д??ͼ????Ŀ??λ??
                      font.legend = 11,#ͼ??????????С
                      #font.title = 12,font.x = 10,font.y = 10,#??????????????С
                      #??ͼ???ϱ????ߵͷֽ????ı???��????????sample??��
                      legend.labs=c(paste0(">",round(svsort[fit$n[2],i],2),"(",fit$n[1],")"),
                                    paste0("<",round(svsort[fit$n[2],i],2),"(",fit$n[2],")")),
                      #????pvalue??HR??95% CI
                      #̫С??p value??Ϊp < 0.001
                      pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                                 paste("p = ",round(p.val,3), sep = "")),
                                   HR, CI, sep = "\n"))
  ggsave(paste0(i,".ESCA.pdf"),width = 4,height = 4)}
length(pl)
F09=t(ptable)#write.csv(ptable, "pvalue_output_BLCA.csv", quote = F, row.names = F)


#---------------------------------------------GBM------------------------------------
svdata <- read.table("3score.time.GBM.txt",sep="\t",header = T,row.names = 1)
dim(svdata)
pct = 0.1
expr <- svdata[,3:ncol(svdata)]
#expr <- log2(expr + 1)
colnames(expr) <- make.names(colnames(expr))
svdata <- cbind.data.frame(svdata[,1:2],expr)
res.cut <- surv_cutpoint(svdata, time = "futime", 
                         event = "fustat", 
                         variables = names(svdata)[3:ncol(svdata)], 
                         minprop = 0.3) #Ĭ??????sample???ܵ???30%
res.cat <- surv_categorize(res.cut)
res.cat$futime=res.cat$futime/30   #####ʱ?�?????????
my.surv <- Surv(res.cat$futime, res.cat$fustat)
pl <- pl2 <- list() # ??ʼ?????????ߺ??¼??????????ߵ??б?
ptable <- data.frame(matrix(NA,ncol(svdata)-2,5))
colnames(ptable) <- c("ID", "pvalue", "HR", "CIlow", "CIup")
n = 0
for (i in colnames(res.cat)[3:ncol(svdata)]) {
  group <- res.cat[,i] 
  survival_dat <- data.frame(group = group)
  fit <- survfit(my.surv ~ group)
  group <- factor(group, levels = c("low", "high"))
  data.survdiff <- survdiff(my.surv ~ group)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  CancerType="GBM"
  n = n + 1
  ptable$ID[n] <- i
  ptable$pvalue[n] <- p.val
  ptable$HR[n] <- HR
  ptable$CIlow[n] <- low95
  ptable$CIup[n] <- up95
  ptable$CancerType[n] <- CancerType
  HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
  CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
  svsort <- svdata[order(svdata[,i]),]
  pl[[i]]<-ggsurvplot(fit, data = survival_dat ,
                      #ggtheme = theme_bw(), #??Ҫ??????????????
                      conf.int = T, #???????????????????????????Ͱ?F?ĳ?T
                      #conf.int.style = "step",#?????????????ͣ????ɸ?Ϊribbon
                      censor = F, #????ʾ?۲?ֵ???ڵ?λ??
                      palette = c("#e64b35","#1b90be"), #?ߵ???ɫ??Ӧ?ߡ???
                      xlab="Time(Months)",
                      legend = c(0.8, 0.9),
                      legend.title = i,#??????д??ͼ????Ŀ??λ??
                      font.legend = 11,#ͼ??????????С
                      #font.title = 12,font.x = 10,font.y = 10,#??????????????С
                      #??ͼ???ϱ????ߵͷֽ????ı???��????????sample??��
                      legend.labs=c(paste0(">",round(svsort[fit$n[2],i],2),"(",fit$n[1],")"),
                                    paste0("<",round(svsort[fit$n[2],i],2),"(",fit$n[2],")")),
                      #????pvalue??HR??95% CI
                      #̫С??p value??Ϊp < 0.001
                      pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                                 paste("p = ",round(p.val,3), sep = "")),
                                   HR, CI, sep = "\n"))
  ggsave(paste0(i,".GBM.pdf"),width = 4,height = 4)}
length(pl)
F10=t(ptable)#write.csv(ptable, "pvalue_output_BLCA.csv", quote = F, row.names = F)


#---------------------------------------------HNSC------------------------------------
svdata <- read.table("3score.time.HNSC.txt",sep="\t",header = T,row.names = 1)
dim(svdata)
pct = 0.1
expr <- svdata[,3:ncol(svdata)]
#expr <- log2(expr + 1)
colnames(expr) <- make.names(colnames(expr))
svdata <- cbind.data.frame(svdata[,1:2],expr)
res.cut <- surv_cutpoint(svdata, time = "futime", 
                         event = "fustat", 
                         variables = names(svdata)[3:ncol(svdata)], 
                         minprop = 0.3) #Ĭ??????sample???ܵ???30%
res.cat <- surv_categorize(res.cut)
res.cat$futime=res.cat$futime/30   #####ʱ?�?????????
my.surv <- Surv(res.cat$futime, res.cat$fustat)
pl <- pl2 <- list() # ??ʼ?????????ߺ??¼??????????ߵ??б?
ptable <- data.frame(matrix(NA,ncol(svdata)-2,5))
colnames(ptable) <- c("ID", "pvalue", "HR", "CIlow", "CIup")
n = 0
for (i in colnames(res.cat)[3:ncol(svdata)]) {
  group <- res.cat[,i] 
  survival_dat <- data.frame(group = group)
  fit <- survfit(my.surv ~ group)
  group <- factor(group, levels = c("low", "high"))
  data.survdiff <- survdiff(my.surv ~ group)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  CancerType="HNSC"
  n = n + 1
  ptable$ID[n] <- i
  ptable$pvalue[n] <- p.val
  ptable$HR[n] <- HR
  ptable$CIlow[n] <- low95
  ptable$CIup[n] <- up95
  ptable$CancerType[n] <- CancerType
  HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
  CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
  svsort <- svdata[order(svdata[,i]),]
  pl[[i]]<-ggsurvplot(fit, data = survival_dat ,
                      #ggtheme = theme_bw(), #??Ҫ??????????????
                      conf.int = T, #???????????????????????????Ͱ?F?ĳ?T
                      #conf.int.style = "step",#?????????????ͣ????ɸ?Ϊribbon
                      censor = F, #????ʾ?۲?ֵ???ڵ?λ??
                      palette = c("#e64b35","#1b90be"), #?ߵ???ɫ??Ӧ?ߡ???
                      xlab="Time(Months)",
                      legend = c(0.8, 0.9),
                      legend.title = i,#??????д??ͼ????Ŀ??λ??
                      font.legend = 11,#ͼ??????????С
                      #font.title = 12,font.x = 10,font.y = 10,#??????????????С
                      #??ͼ???ϱ????ߵͷֽ????ı???��????????sample??��
                      legend.labs=c(paste0(">",round(svsort[fit$n[2],i],2),"(",fit$n[1],")"),
                                    paste0("<",round(svsort[fit$n[2],i],2),"(",fit$n[2],")")),
                      #????pvalue??HR??95% CI
                      #̫С??p value??Ϊp < 0.001
                      pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                                 paste("p = ",round(p.val,3), sep = "")),
                                   HR, CI, sep = "\n"))
  ggsave(paste0(i,".HNSC.pdf"),width = 4,height = 4)}
length(pl)
F11=t(ptable)#write.csv(ptable, "pvalue_output_BLCA.csv", quote = F, row.names = F)


#---------------------------------------------KICH------------------------------------
svdata <- read.table("3score.time.KICH.txt",sep="\t",header = T,row.names = 1)
dim(svdata)
pct = 0.1
expr <- svdata[,3:ncol(svdata)]
#expr <- log2(expr + 1)
colnames(expr) <- make.names(colnames(expr))
svdata <- cbind.data.frame(svdata[,1:2],expr)
res.cut <- surv_cutpoint(svdata, time = "futime", 
                         event = "fustat", 
                         variables = names(svdata)[3:ncol(svdata)], 
                         minprop = 0.3) #Ĭ??????sample???ܵ???30%
res.cat <- surv_categorize(res.cut)
res.cat$futime=res.cat$futime/30   #####ʱ?�?????????
my.surv <- Surv(res.cat$futime, res.cat$fustat)
pl <- pl2 <- list() # ??ʼ?????????ߺ??¼??????????ߵ??б?
ptable <- data.frame(matrix(NA,ncol(svdata)-2,5))
colnames(ptable) <- c("ID", "pvalue", "HR", "CIlow", "CIup")
n = 0
for (i in colnames(res.cat)[3:ncol(svdata)]) {
  group <- res.cat[,i] 
  survival_dat <- data.frame(group = group)
  fit <- survfit(my.surv ~ group)
  group <- factor(group, levels = c("low", "high"))
  data.survdiff <- survdiff(my.surv ~ group)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  CancerType="KICH"
  n = n + 1
  ptable$ID[n] <- i
  ptable$pvalue[n] <- p.val
  ptable$HR[n] <- HR
  ptable$CIlow[n] <- low95
  ptable$CIup[n] <- up95
  ptable$CancerType[n] <- CancerType
  HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
  CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
  svsort <- svdata[order(svdata[,i]),]
  pl[[i]]<-ggsurvplot(fit, data = survival_dat ,
                      #ggtheme = theme_bw(), #??Ҫ??????????????
                      conf.int = T, #???????????????????????????Ͱ?F?ĳ?T
                      #conf.int.style = "step",#?????????????ͣ????ɸ?Ϊribbon
                      censor = F, #????ʾ?۲?ֵ???ڵ?λ??
                      palette = c("#e64b35","#1b90be"), #?ߵ???ɫ??Ӧ?ߡ???
                      xlab="Time(Months)",
                      legend = c(0.8, 0.9),
                      legend.title = i,#??????д??ͼ????Ŀ??λ??
                      font.legend = 11,#ͼ??????????С
                      #font.title = 12,font.x = 10,font.y = 10,#??????????????С
                      #??ͼ???ϱ????ߵͷֽ????ı???��????????sample??��
                      legend.labs=c(paste0(">",round(svsort[fit$n[2],i],2),"(",fit$n[1],")"),
                                    paste0("<",round(svsort[fit$n[2],i],2),"(",fit$n[2],")")),
                      #????pvalue??HR??95% CI
                      #̫С??p value??Ϊp < 0.001
                      pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                                 paste("p = ",round(p.val,3), sep = "")),
                                   HR, CI, sep = "\n"))
  ggsave(paste0(i,".KICH.pdf"),width = 4,height = 4)}
length(pl)
F12=t(ptable)#write.csv(ptable, "pvalue_output_BLCA.csv", quote = F, row.names = F)


#---------------------------------------------KIRC------------------------------------
svdata <- read.table("3score.time.KIRC.txt",sep="\t",header = T,row.names = 1)
dim(svdata)
pct = 0.1
expr <- svdata[,3:ncol(svdata)]
#expr <- log2(expr + 1)
colnames(expr) <- make.names(colnames(expr))
svdata <- cbind.data.frame(svdata[,1:2],expr)
res.cut <- surv_cutpoint(svdata, time = "futime", 
                         event = "fustat", 
                         variables = names(svdata)[3:ncol(svdata)], 
                         minprop = 0.3) #Ĭ??????sample???ܵ???30%
res.cat <- surv_categorize(res.cut)
res.cat$futime=res.cat$futime/30   #####ʱ?�?????????
my.surv <- Surv(res.cat$futime, res.cat$fustat)
pl <- pl2 <- list() # ??ʼ?????????ߺ??¼??????????ߵ??б?
ptable <- data.frame(matrix(NA,ncol(svdata)-2,5))
colnames(ptable) <- c("ID", "pvalue", "HR", "CIlow", "CIup")
n = 0
for (i in colnames(res.cat)[3:ncol(svdata)]) {
  group <- res.cat[,i] 
  survival_dat <- data.frame(group = group)
  fit <- survfit(my.surv ~ group)
  group <- factor(group, levels = c("low", "high"))
  data.survdiff <- survdiff(my.surv ~ group)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  CancerType="KIRC"
  n = n + 1
  ptable$ID[n] <- i
  ptable$pvalue[n] <- p.val
  ptable$HR[n] <- HR
  ptable$CIlow[n] <- low95
  ptable$CIup[n] <- up95
  ptable$CancerType[n] <- CancerType
  HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
  CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
  svsort <- svdata[order(svdata[,i]),]
  pl[[i]]<-ggsurvplot(fit, data = survival_dat ,
                      #ggtheme = theme_bw(), #??Ҫ??????????????
                      conf.int = T, #???????????????????????????Ͱ?F?ĳ?T
                      #conf.int.style = "step",#?????????????ͣ????ɸ?Ϊribbon
                      censor = F, #????ʾ?۲?ֵ???ڵ?λ??
                      palette = c("#e64b35","#1b90be"), #?ߵ???ɫ??Ӧ?ߡ???
                      xlab="Time(Months)",
                      legend = c(0.8, 0.9),
                      legend.title = i,#??????д??ͼ????Ŀ??λ??
                      font.legend = 11,#ͼ??????????С
                      #font.title = 12,font.x = 10,font.y = 10,#??????????????С
                      #??ͼ???ϱ????ߵͷֽ????ı???��????????sample??��
                      legend.labs=c(paste0(">",round(svsort[fit$n[2],i],2),"(",fit$n[1],")"),
                                    paste0("<",round(svsort[fit$n[2],i],2),"(",fit$n[2],")")),
                      #????pvalue??HR??95% CI
                      #̫С??p value??Ϊp < 0.001
                      pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                                 paste("p = ",round(p.val,3), sep = "")),
                                   HR, CI, sep = "\n"))
  ggsave(paste0(i,".KIRC.pdf"),width = 4,height = 4)}
length(pl)
F13=t(ptable)#write.csv(ptable, "pvalue_output_BLCA.csv", quote = F, row.names = F)


#---------------------------------------------KIRP------------------------------------
svdata <- read.table("3score.time.KIRP.txt",sep="\t",header = T,row.names = 1)
dim(svdata)
pct = 0.1
expr <- svdata[,3:ncol(svdata)]
#expr <- log2(expr + 1)
colnames(expr) <- make.names(colnames(expr))
svdata <- cbind.data.frame(svdata[,1:2],expr)
res.cut <- surv_cutpoint(svdata, time = "futime", 
                         event = "fustat", 
                         variables = names(svdata)[3:ncol(svdata)], 
                         minprop = 0.3) #Ĭ??????sample???ܵ???30%
res.cat <- surv_categorize(res.cut)
res.cat$futime=res.cat$futime/30   #####ʱ?�?????????
my.surv <- Surv(res.cat$futime, res.cat$fustat)
pl <- pl2 <- list() # ??ʼ?????????ߺ??¼??????????ߵ??б?
ptable <- data.frame(matrix(NA,ncol(svdata)-2,5))
colnames(ptable) <- c("ID", "pvalue", "HR", "CIlow", "CIup")
n = 0
for (i in colnames(res.cat)[3:ncol(svdata)]) {
  group <- res.cat[,i] 
  survival_dat <- data.frame(group = group)
  fit <- survfit(my.surv ~ group)
  group <- factor(group, levels = c("low", "high"))
  data.survdiff <- survdiff(my.surv ~ group)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  CancerType="KIRP"
  n = n + 1
  ptable$ID[n] <- i
  ptable$pvalue[n] <- p.val
  ptable$HR[n] <- HR
  ptable$CIlow[n] <- low95
  ptable$CIup[n] <- up95
  ptable$CancerType[n] <- CancerType
  HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
  CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
  svsort <- svdata[order(svdata[,i]),]
  pl[[i]]<-ggsurvplot(fit, data = survival_dat ,
                      #ggtheme = theme_bw(), #??Ҫ??????????????
                      conf.int = T, #???????????????????????????Ͱ?F?ĳ?T
                      #conf.int.style = "step",#?????????????ͣ????ɸ?Ϊribbon
                      censor = F, #????ʾ?۲?ֵ???ڵ?λ??
                      palette = c("#e64b35","#1b90be"), #?ߵ???ɫ??Ӧ?ߡ???
                      xlab="Time(Months)",
                      legend = c(0.8, 0.9),
                      legend.title = i,#??????д??ͼ????Ŀ??λ??
                      font.legend = 11,#ͼ??????????С
                      #font.title = 12,font.x = 10,font.y = 10,#??????????????С
                      #??ͼ???ϱ????ߵͷֽ????ı???��????????sample??��
                      legend.labs=c(paste0(">",round(svsort[fit$n[2],i],2),"(",fit$n[1],")"),
                                    paste0("<",round(svsort[fit$n[2],i],2),"(",fit$n[2],")")),
                      #????pvalue??HR??95% CI
                      #̫С??p value??Ϊp < 0.001
                      pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                                 paste("p = ",round(p.val,3), sep = "")),
                                   HR, CI, sep = "\n"))
  ggsave(paste0(i,".KIRP.pdf"),width = 4,height = 4)}
length(pl)
F14=t(ptable)#write.csv(ptable, "pvalue_output_BLCA.csv", quote = F, row.names = F)


#---------------------------------------------LAML------------------------------------
svdata <- read.table("3score.time.LAML.txt",sep="\t",header = T,row.names = 1)
dim(svdata)
pct = 0.1
expr <- svdata[,3:ncol(svdata)]
#expr <- log2(expr + 1)
colnames(expr) <- make.names(colnames(expr))
svdata <- cbind.data.frame(svdata[,1:2],expr)
res.cut <- surv_cutpoint(svdata, time = "futime", 
                         event = "fustat", 
                         variables = names(svdata)[3:ncol(svdata)], 
                         minprop = 0.3) #Ĭ??????sample???ܵ???30%
res.cat <- surv_categorize(res.cut)
res.cat$futime=res.cat$futime/30   #####ʱ?�?????????
my.surv <- Surv(res.cat$futime, res.cat$fustat)
pl <- pl2 <- list() # ??ʼ?????????ߺ??¼??????????ߵ??б?
ptable <- data.frame(matrix(NA,ncol(svdata)-2,5))
colnames(ptable) <- c("ID", "pvalue", "HR", "CIlow", "CIup")
n = 0
for (i in colnames(res.cat)[3:ncol(svdata)]) {
  group <- res.cat[,i] 
  survival_dat <- data.frame(group = group)
  fit <- survfit(my.surv ~ group)
  group <- factor(group, levels = c("low", "high"))
  data.survdiff <- survdiff(my.surv ~ group)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  CancerType="LAML"
  n = n + 1
  ptable$ID[n] <- i
  ptable$pvalue[n] <- p.val
  ptable$HR[n] <- HR
  ptable$CIlow[n] <- low95
  ptable$CIup[n] <- up95
  ptable$CancerType[n] <- CancerType
  HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
  CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
  svsort <- svdata[order(svdata[,i]),]
  pl[[i]]<-ggsurvplot(fit, data = survival_dat ,
                      #ggtheme = theme_bw(), #??Ҫ??????????????
                      conf.int = T, #???????????????????????????Ͱ?F?ĳ?T
                      #conf.int.style = "step",#?????????????ͣ????ɸ?Ϊribbon
                      censor = F, #????ʾ?۲?ֵ???ڵ?λ??
                      palette = c("#e64b35","#1b90be"), #?ߵ???ɫ??Ӧ?ߡ???
                      xlab="Time(Months)",
                      legend = c(0.8, 0.9),
                      legend.title = i,#??????д??ͼ????Ŀ??λ??
                      font.legend = 11,#ͼ??????????С
                      #font.title = 12,font.x = 10,font.y = 10,#??????????????С
                      #??ͼ???ϱ????ߵͷֽ????ı???��????????sample??��
                      legend.labs=c(paste0(">",round(svsort[fit$n[2],i],2),"(",fit$n[1],")"),
                                    paste0("<",round(svsort[fit$n[2],i],2),"(",fit$n[2],")")),
                      #????pvalue??HR??95% CI
                      #̫С??p value??Ϊp < 0.001
                      pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                                 paste("p = ",round(p.val,3), sep = "")),
                                   HR, CI, sep = "\n"))
  ggsave(paste0(i,".LAML.pdf"),width = 4,height = 4)}
length(pl)
F15=t(ptable)#write.csv(ptable, "pvalue_output_BLCA.csv", quote = F, row.names = F)


#---------------------------------------------LGG------------------------------------
svdata <- read.table("3score.time.LGG.txt",sep="\t",header = T,row.names = 1)
dim(svdata)
pct = 0.1
expr <- svdata[,3:ncol(svdata)]
#expr <- log2(expr + 1)
colnames(expr) <- make.names(colnames(expr))
svdata <- cbind.data.frame(svdata[,1:2],expr)
res.cut <- surv_cutpoint(svdata, time = "futime", 
                         event = "fustat", 
                         variables = names(svdata)[3:ncol(svdata)], 
                         minprop = 0.3) #Ĭ??????sample???ܵ???30%
res.cat <- surv_categorize(res.cut)
res.cat$futime=res.cat$futime/30   #####ʱ?�?????????
my.surv <- Surv(res.cat$futime, res.cat$fustat)
pl <- pl2 <- list() # ??ʼ?????????ߺ??¼??????????ߵ??б?
ptable <- data.frame(matrix(NA,ncol(svdata)-2,5))
colnames(ptable) <- c("ID", "pvalue", "HR", "CIlow", "CIup")
n = 0
for (i in colnames(res.cat)[3:ncol(svdata)]) {
  group <- res.cat[,i] 
  survival_dat <- data.frame(group = group)
  fit <- survfit(my.surv ~ group)
  group <- factor(group, levels = c("low", "high"))
  data.survdiff <- survdiff(my.surv ~ group)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  CancerType="LGG"
  n = n + 1
  ptable$ID[n] <- i
  ptable$pvalue[n] <- p.val
  ptable$HR[n] <- HR
  ptable$CIlow[n] <- low95
  ptable$CIup[n] <- up95
  ptable$CancerType[n] <- CancerType
  HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
  CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
  svsort <- svdata[order(svdata[,i]),]
  pl[[i]]<-ggsurvplot(fit, data = survival_dat ,
                      #ggtheme = theme_bw(), #??Ҫ??????????????
                      conf.int = T, #???????????????????????????Ͱ?F?ĳ?T
                      #conf.int.style = "step",#?????????????ͣ????ɸ?Ϊribbon
                      censor = F, #????ʾ?۲?ֵ???ڵ?λ??
                      palette = c("#e64b35","#1b90be"), #?ߵ???ɫ??Ӧ?ߡ???
                      xlab="Time(Months)",
                      legend = c(0.8, 0.9),
                      legend.title = i,#??????д??ͼ????Ŀ??λ??
                      font.legend = 11,#ͼ??????????С
                      #font.title = 12,font.x = 10,font.y = 10,#??????????????С
                      #??ͼ???ϱ????ߵͷֽ????ı???��????????sample??��
                      legend.labs=c(paste0(">",round(svsort[fit$n[2],i],2),"(",fit$n[1],")"),
                                    paste0("<",round(svsort[fit$n[2],i],2),"(",fit$n[2],")")),
                      #????pvalue??HR??95% CI
                      #̫С??p value??Ϊp < 0.001
                      pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                                 paste("p = ",round(p.val,3), sep = "")),
                                   HR, CI, sep = "\n"))
  ggsave(paste0(i,".LGG.pdf"),width = 4,height = 4)}
length(pl)
F16=t(ptable)#write.csv(ptable, "pvalue_output_BLCA.csv", quote = F, row.names = F)


#---------------------------------------------LIHC------------------------------------
svdata <- read.table("3score.time.LIHC.txt",sep="\t",header = T,row.names = 1)
dim(svdata)
pct = 0.1
expr <- svdata[,3:ncol(svdata)]
#expr <- log2(expr + 1)
colnames(expr) <- make.names(colnames(expr))
svdata <- cbind.data.frame(svdata[,1:2],expr)
res.cut <- surv_cutpoint(svdata, time = "futime", 
                         event = "fustat", 
                         variables = names(svdata)[3:ncol(svdata)], 
                         minprop = 0.3) #Ĭ??????sample???ܵ???30%
res.cat <- surv_categorize(res.cut)
res.cat$futime=res.cat$futime/30   #####ʱ?�?????????
my.surv <- Surv(res.cat$futime, res.cat$fustat)
pl <- pl2 <- list() # ??ʼ?????????ߺ??¼??????????ߵ??б?
ptable <- data.frame(matrix(NA,ncol(svdata)-2,5))
colnames(ptable) <- c("ID", "pvalue", "HR", "CIlow", "CIup")
n = 0
for (i in colnames(res.cat)[3:ncol(svdata)]) {
  group <- res.cat[,i] 
  survival_dat <- data.frame(group = group)
  fit <- survfit(my.surv ~ group)
  group <- factor(group, levels = c("low", "high"))
  data.survdiff <- survdiff(my.surv ~ group)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  CancerType="LIHC"
  n = n + 1
  ptable$ID[n] <- i
  ptable$pvalue[n] <- p.val
  ptable$HR[n] <- HR
  ptable$CIlow[n] <- low95
  ptable$CIup[n] <- up95
  ptable$CancerType[n] <- CancerType
  HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
  CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
  svsort <- svdata[order(svdata[,i]),]
  pl[[i]]<-ggsurvplot(fit, data = survival_dat ,
                      #ggtheme = theme_bw(), #??Ҫ??????????????
                      conf.int = T, #???????????????????????????Ͱ?F?ĳ?T
                      #conf.int.style = "step",#?????????????ͣ????ɸ?Ϊribbon
                      censor = F, #????ʾ?۲?ֵ???ڵ?λ??
                      palette = c("#e64b35","#1b90be"), #?ߵ???ɫ??Ӧ?ߡ???
                      xlab="Time(Months)",
                      legend = c(0.8, 0.9),
                      legend.title = i,#??????д??ͼ????Ŀ??λ??
                      font.legend = 11,#ͼ??????????С
                      #font.title = 12,font.x = 10,font.y = 10,#??????????????С
                      #??ͼ???ϱ????ߵͷֽ????ı???��????????sample??��
                      legend.labs=c(paste0(">",round(svsort[fit$n[2],i],2),"(",fit$n[1],")"),
                                    paste0("<",round(svsort[fit$n[2],i],2),"(",fit$n[2],")")),
                      #????pvalue??HR??95% CI
                      #̫С??p value??Ϊp < 0.001
                      pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                                 paste("p = ",round(p.val,3), sep = "")),
                                   HR, CI, sep = "\n"))
  ggsave(paste0(i,".LIHC.pdf"),width = 4,height = 4)}
length(pl)
F17=t(ptable)#write.csv(ptable, "pvalue_output_BLCA.csv", quote = F, row.names = F)


#---------------------------------------------LUAD------------------------------------
svdata <- read.table("3score.time.LUAD.txt",sep="\t",header = T,row.names = 1)
dim(svdata)
pct = 0.1
expr <- svdata[,3:ncol(svdata)]
#expr <- log2(expr + 1)
colnames(expr) <- make.names(colnames(expr))
svdata <- cbind.data.frame(svdata[,1:2],expr)
res.cut <- surv_cutpoint(svdata, time = "futime", 
                         event = "fustat", 
                         variables = names(svdata)[3:ncol(svdata)], 
                         minprop = 0.3) #Ĭ??????sample???ܵ???30%
res.cat <- surv_categorize(res.cut)
res.cat$futime=res.cat$futime/30   #####ʱ?�?????????
my.surv <- Surv(res.cat$futime, res.cat$fustat)
pl <- pl2 <- list() # ??ʼ?????????ߺ??¼??????????ߵ??б?
ptable <- data.frame(matrix(NA,ncol(svdata)-2,5))
colnames(ptable) <- c("ID", "pvalue", "HR", "CIlow", "CIup")
n = 0
for (i in colnames(res.cat)[3:ncol(svdata)]) {
  group <- res.cat[,i] 
  survival_dat <- data.frame(group = group)
  fit <- survfit(my.surv ~ group)
  group <- factor(group, levels = c("low", "high"))
  data.survdiff <- survdiff(my.surv ~ group)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  CancerType="LUAD"
  n = n + 1
  ptable$ID[n] <- i
  ptable$pvalue[n] <- p.val
  ptable$HR[n] <- HR
  ptable$CIlow[n] <- low95
  ptable$CIup[n] <- up95
  ptable$CancerType[n] <- CancerType
  HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
  CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
  svsort <- svdata[order(svdata[,i]),]
  pl[[i]]<-ggsurvplot(fit, data = survival_dat ,
                      #ggtheme = theme_bw(), #??Ҫ??????????????
                      conf.int = T, #???????????????????????????Ͱ?F?ĳ?T
                      #conf.int.style = "step",#?????????????ͣ????ɸ?Ϊribbon
                      censor = F, #????ʾ?۲?ֵ???ڵ?λ??
                      palette = c("#e64b35","#1b90be"), #?ߵ???ɫ??Ӧ?ߡ???
                      xlab="Time(Months)",
                      legend = c(0.8, 0.9),
                      legend.title = i,#??????д??ͼ????Ŀ??λ??
                      font.legend = 11,#ͼ??????????С
                      #font.title = 12,font.x = 10,font.y = 10,#??????????????С
                      #??ͼ???ϱ????ߵͷֽ????ı???��????????sample??��
                      legend.labs=c(paste0(">",round(svsort[fit$n[2],i],2),"(",fit$n[1],")"),
                                    paste0("<",round(svsort[fit$n[2],i],2),"(",fit$n[2],")")),
                      #????pvalue??HR??95% CI
                      #̫С??p value??Ϊp < 0.001
                      pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                                 paste("p = ",round(p.val,3), sep = "")),
                                   HR, CI, sep = "\n"))
  ggsave(paste0(i,".LUAD.pdf"),width = 4,height = 4)}
length(pl)
F18=t(ptable)#write.csv(ptable, "pvalue_output_BLCA.csv", quote = F, row.names = F)


#---------------------------------------------LUSC------------------------------------
svdata <- read.table("3score.time.LUSC.txt",sep="\t",header = T,row.names = 1)
dim(svdata)
pct = 0.1
expr <- svdata[,3:ncol(svdata)]
#expr <- log2(expr + 1)
colnames(expr) <- make.names(colnames(expr))
svdata <- cbind.data.frame(svdata[,1:2],expr)
res.cut <- surv_cutpoint(svdata, time = "futime", 
                         event = "fustat", 
                         variables = names(svdata)[3:ncol(svdata)], 
                         minprop = 0.3) #Ĭ??????sample???ܵ???30%
res.cat <- surv_categorize(res.cut)
res.cat$futime=res.cat$futime/30   #####ʱ?�?????????
my.surv <- Surv(res.cat$futime, res.cat$fustat)
pl <- pl2 <- list() # ??ʼ?????????ߺ??¼??????????ߵ??б?
ptable <- data.frame(matrix(NA,ncol(svdata)-2,5))
colnames(ptable) <- c("ID", "pvalue", "HR", "CIlow", "CIup")
n = 0
for (i in colnames(res.cat)[3:ncol(svdata)]) {
  group <- res.cat[,i] 
  survival_dat <- data.frame(group = group)
  fit <- survfit(my.surv ~ group)
  group <- factor(group, levels = c("low", "high"))
  data.survdiff <- survdiff(my.surv ~ group)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  CancerType="LUSC"
  n = n + 1
  ptable$ID[n] <- i
  ptable$pvalue[n] <- p.val
  ptable$HR[n] <- HR
  ptable$CIlow[n] <- low95
  ptable$CIup[n] <- up95
  ptable$CancerType[n] <- CancerType
  HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
  CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
  svsort <- svdata[order(svdata[,i]),]
  pl[[i]]<-ggsurvplot(fit, data = survival_dat ,
                      #ggtheme = theme_bw(), #??Ҫ??????????????
                      conf.int = T, #???????????????????????????Ͱ?F?ĳ?T
                      #conf.int.style = "step",#?????????????ͣ????ɸ?Ϊribbon
                      censor = F, #????ʾ?۲?ֵ???ڵ?λ??
                      palette = c("#e64b35","#1b90be"), #?ߵ???ɫ??Ӧ?ߡ???
                      xlab="Time(Months)",
                      legend = c(0.8, 0.9),
                      legend.title = i,#??????д??ͼ????Ŀ??λ??
                      font.legend = 11,#ͼ??????????С
                      #font.title = 12,font.x = 10,font.y = 10,#??????????????С
                      #??ͼ???ϱ????ߵͷֽ????ı???��????????sample??��
                      legend.labs=c(paste0(">",round(svsort[fit$n[2],i],2),"(",fit$n[1],")"),
                                    paste0("<",round(svsort[fit$n[2],i],2),"(",fit$n[2],")")),
                      #????pvalue??HR??95% CI
                      #̫С??p value??Ϊp < 0.001
                      pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                                 paste("p = ",round(p.val,3), sep = "")),
                                   HR, CI, sep = "\n"))
  ggsave(paste0(i,".LUSC.pdf"),width = 4,height = 4)}
length(pl)
F19=t(ptable)#write.csv(ptable, "pvalue_output_BLCA.csv", quote = F, row.names = F)


#---------------------------------------------MESO------------------------------------
svdata <- read.table("3score.time.MESO.txt",sep="\t",header = T,row.names = 1)
dim(svdata)
pct = 0.1
expr <- svdata[,3:ncol(svdata)]
#expr <- log2(expr + 1)
colnames(expr) <- make.names(colnames(expr))
svdata <- cbind.data.frame(svdata[,1:2],expr)
res.cut <- surv_cutpoint(svdata, time = "futime", 
                         event = "fustat", 
                         variables = names(svdata)[3:ncol(svdata)], 
                         minprop = 0.3) #Ĭ??????sample???ܵ???30%
res.cat <- surv_categorize(res.cut)
res.cat$futime=res.cat$futime/30   #####ʱ?�?????????
my.surv <- Surv(res.cat$futime, res.cat$fustat)
pl <- pl2 <- list() # ??ʼ?????????ߺ??¼??????????ߵ??б?
ptable <- data.frame(matrix(NA,ncol(svdata)-2,5))
colnames(ptable) <- c("ID", "pvalue", "HR", "CIlow", "CIup")
n = 0
for (i in colnames(res.cat)[3:ncol(svdata)]) {
  group <- res.cat[,i] 
  survival_dat <- data.frame(group = group)
  fit <- survfit(my.surv ~ group)
  group <- factor(group, levels = c("low", "high"))
  data.survdiff <- survdiff(my.surv ~ group)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  CancerType="MESO"
  n = n + 1
  ptable$ID[n] <- i
  ptable$pvalue[n] <- p.val
  ptable$HR[n] <- HR
  ptable$CIlow[n] <- low95
  ptable$CIup[n] <- up95
  ptable$CancerType[n] <- CancerType
  HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
  CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
  svsort <- svdata[order(svdata[,i]),]
  pl[[i]]<-ggsurvplot(fit, data = survival_dat ,
                      #ggtheme = theme_bw(), #??Ҫ??????????????
                      conf.int = T, #???????????????????????????Ͱ?F?ĳ?T
                      #conf.int.style = "step",#?????????????ͣ????ɸ?Ϊribbon
                      censor = F, #????ʾ?۲?ֵ???ڵ?λ??
                      palette = c("#e64b35","#1b90be"), #?ߵ???ɫ??Ӧ?ߡ???
                      xlab="Time(Months)",
                      legend = c(0.8, 0.9),
                      legend.title = i,#??????д??ͼ????Ŀ??λ??
                      font.legend = 11,#ͼ??????????С
                      #font.title = 12,font.x = 10,font.y = 10,#??????????????С
                      #??ͼ???ϱ????ߵͷֽ????ı???��????????sample??��
                      legend.labs=c(paste0(">",round(svsort[fit$n[2],i],2),"(",fit$n[1],")"),
                                    paste0("<",round(svsort[fit$n[2],i],2),"(",fit$n[2],")")),
                      #????pvalue??HR??95% CI
                      #̫С??p value??Ϊp < 0.001
                      pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                                 paste("p = ",round(p.val,3), sep = "")),
                                   HR, CI, sep = "\n"))
  ggsave(paste0(i,".MESO.pdf"),width = 4,height = 4)}
length(pl)
F20=t(ptable)#write.csv(ptable, "pvalue_output_BLCA.csv", quote = F, row.names = F)


#---------------------------------------------OV------------------------------------
svdata <- read.table("3score.time.OV.txt",sep="\t",header = T,row.names = 1)
dim(svdata)
pct = 0.1
expr <- svdata[,3:ncol(svdata)]
#expr <- log2(expr + 1)
colnames(expr) <- make.names(colnames(expr))
svdata <- cbind.data.frame(svdata[,1:2],expr)
res.cut <- surv_cutpoint(svdata, time = "futime", 
                         event = "fustat", 
                         variables = names(svdata)[3:ncol(svdata)], 
                         minprop = 0.3) #Ĭ??????sample???ܵ???30%
res.cat <- surv_categorize(res.cut)
res.cat$futime=res.cat$futime/30   #####ʱ?�?????????
my.surv <- Surv(res.cat$futime, res.cat$fustat)
pl <- pl2 <- list() # ??ʼ?????????ߺ??¼??????????ߵ??б?
ptable <- data.frame(matrix(NA,ncol(svdata)-2,5))
colnames(ptable) <- c("ID", "pvalue", "HR", "CIlow", "CIup")
n = 0
for (i in colnames(res.cat)[3:ncol(svdata)]) {
  group <- res.cat[,i] 
  survival_dat <- data.frame(group = group)
  fit <- survfit(my.surv ~ group)
  group <- factor(group, levels = c("low", "high"))
  data.survdiff <- survdiff(my.surv ~ group)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  CancerType="OV"
  n = n + 1
  ptable$ID[n] <- i
  ptable$pvalue[n] <- p.val
  ptable$HR[n] <- HR
  ptable$CIlow[n] <- low95
  ptable$CIup[n] <- up95
  ptable$CancerType[n] <- CancerType
  HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
  CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
  svsort <- svdata[order(svdata[,i]),]
  pl[[i]]<-ggsurvplot(fit, data = survival_dat ,
                      #ggtheme = theme_bw(), #??Ҫ??????????????
                      conf.int = T, #???????????????????????????Ͱ?F?ĳ?T
                      #conf.int.style = "step",#?????????????ͣ????ɸ?Ϊribbon
                      censor = F, #????ʾ?۲?ֵ???ڵ?λ??
                      palette = c("#e64b35","#1b90be"), #?ߵ???ɫ??Ӧ?ߡ???
                      xlab="Time(Months)",
                      legend = c(0.8, 0.9),
                      legend.title = i,#??????д??ͼ????Ŀ??λ??
                      font.legend = 11,#ͼ??????????С
                      #font.title = 12,font.x = 10,font.y = 10,#??????????????С
                      #??ͼ???ϱ????ߵͷֽ????ı???��????????sample??��
                      legend.labs=c(paste0(">",round(svsort[fit$n[2],i],2),"(",fit$n[1],")"),
                                    paste0("<",round(svsort[fit$n[2],i],2),"(",fit$n[2],")")),
                      #????pvalue??HR??95% CI
                      #̫С??p value??Ϊp < 0.001
                      pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                                 paste("p = ",round(p.val,3), sep = "")),
                                   HR, CI, sep = "\n"))
  ggsave(paste0(i,".OV.pdf"),width = 4,height = 4)}
length(pl)
F21=t(ptable)#write.csv(ptable, "pvalue_output_BLCA.csv", quote = F, row.names = F)


#---------------------------------------------PAAD------------------------------------
svdata <- read.table("3score.time.PAAD.txt",sep="\t",header = T,row.names = 1)
dim(svdata)
pct = 0.1
expr <- svdata[,3:ncol(svdata)]
#expr <- log2(expr + 1)
colnames(expr) <- make.names(colnames(expr))
svdata <- cbind.data.frame(svdata[,1:2],expr)
res.cut <- surv_cutpoint(svdata, time = "futime", 
                         event = "fustat", 
                         variables = names(svdata)[3:ncol(svdata)], 
                         minprop = 0.3) #Ĭ??????sample???ܵ???30%
res.cat <- surv_categorize(res.cut)
res.cat$futime=res.cat$futime/30   #####ʱ?�?????????
my.surv <- Surv(res.cat$futime, res.cat$fustat)
pl <- pl2 <- list() # ??ʼ?????????ߺ??¼??????????ߵ??б?
ptable <- data.frame(matrix(NA,ncol(svdata)-2,5))
colnames(ptable) <- c("ID", "pvalue", "HR", "CIlow", "CIup")
n = 0
for (i in colnames(res.cat)[3:ncol(svdata)]) {
  group <- res.cat[,i] 
  survival_dat <- data.frame(group = group)
  fit <- survfit(my.surv ~ group)
  group <- factor(group, levels = c("low", "high"))
  data.survdiff <- survdiff(my.surv ~ group)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  CancerType="PAAD"
  n = n + 1
  ptable$ID[n] <- i
  ptable$pvalue[n] <- p.val
  ptable$HR[n] <- HR
  ptable$CIlow[n] <- low95
  ptable$CIup[n] <- up95
  ptable$CancerType[n] <- CancerType
  HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
  CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
  svsort <- svdata[order(svdata[,i]),]
  pl[[i]]<-ggsurvplot(fit, data = survival_dat ,
                      #ggtheme = theme_bw(), #??Ҫ??????????????
                      conf.int = T, #???????????????????????????Ͱ?F?ĳ?T
                      #conf.int.style = "step",#?????????????ͣ????ɸ?Ϊribbon
                      censor = F, #????ʾ?۲?ֵ???ڵ?λ??
                      palette = c("#e64b35","#1b90be"), #?ߵ???ɫ??Ӧ?ߡ???
                      xlab="Time(Months)",
                      legend = c(0.8, 0.9),
                      legend.title = i,#??????д??ͼ????Ŀ??λ??
                      font.legend = 11,#ͼ??????????С
                      #font.title = 12,font.x = 10,font.y = 10,#??????????????С
                      #??ͼ???ϱ????ߵͷֽ????ı???��????????sample??��
                      legend.labs=c(paste0(">",round(svsort[fit$n[2],i],2),"(",fit$n[1],")"),
                                    paste0("<",round(svsort[fit$n[2],i],2),"(",fit$n[2],")")),
                      #????pvalue??HR??95% CI
                      #̫С??p value??Ϊp < 0.001
                      pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                                 paste("p = ",round(p.val,3), sep = "")),
                                   HR, CI, sep = "\n"))
  ggsave(paste0(i,".PAAD.pdf"),width = 4,height = 4)}
length(pl)
F22=t(ptable)#write.csv(ptable, "pvalue_output_BLCA.csv", quote = F, row.names = F)


#---------------------------------------------PCPG------------------------------------
svdata <- read.table("3score.time.PCPG.txt",sep="\t",header = T,row.names = 1)
dim(svdata)
pct = 0.1
expr <- svdata[,3:ncol(svdata)]
#expr <- log2(expr + 1)
colnames(expr) <- make.names(colnames(expr))
svdata <- cbind.data.frame(svdata[,1:2],expr)
res.cut <- surv_cutpoint(svdata, time = "futime", 
                         event = "fustat", 
                         variables = names(svdata)[3:ncol(svdata)], 
                         minprop = 0.3) #Ĭ??????sample???ܵ???30%
res.cat <- surv_categorize(res.cut)
res.cat$futime=res.cat$futime/30   #####ʱ?�?????????
my.surv <- Surv(res.cat$futime, res.cat$fustat)
pl <- pl2 <- list() # ??ʼ?????????ߺ??¼??????????ߵ??б?
ptable <- data.frame(matrix(NA,ncol(svdata)-2,5))
colnames(ptable) <- c("ID", "pvalue", "HR", "CIlow", "CIup")
n = 0
for (i in colnames(res.cat)[3:ncol(svdata)]) {
  group <- res.cat[,i] 
  survival_dat <- data.frame(group = group)
  fit <- survfit(my.surv ~ group)
  group <- factor(group, levels = c("low", "high"))
  data.survdiff <- survdiff(my.surv ~ group)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  CancerType="PCPG"
  n = n + 1
  ptable$ID[n] <- i
  ptable$pvalue[n] <- p.val
  ptable$HR[n] <- HR
  ptable$CIlow[n] <- low95
  ptable$CIup[n] <- up95
  ptable$CancerType[n] <- CancerType
  HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
  CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
  svsort <- svdata[order(svdata[,i]),]
  pl[[i]]<-ggsurvplot(fit, data = survival_dat ,
                      #ggtheme = theme_bw(), #??Ҫ??????????????
                      conf.int = T, #???????????????????????????Ͱ?F?ĳ?T
                      #conf.int.style = "step",#?????????????ͣ????ɸ?Ϊribbon
                      censor = F, #????ʾ?۲?ֵ???ڵ?λ??
                      palette = c("#e64b35","#1b90be"), #?ߵ???ɫ??Ӧ?ߡ???
                      xlab="Time(Months)",
                      legend = c(0.8, 0.9),
                      legend.title = i,#??????д??ͼ????Ŀ??λ??
                      font.legend = 11,#ͼ??????????С
                      #font.title = 12,font.x = 10,font.y = 10,#??????????????С
                      #??ͼ???ϱ????ߵͷֽ????ı???��????????sample??��
                      legend.labs=c(paste0(">",round(svsort[fit$n[2],i],2),"(",fit$n[1],")"),
                                    paste0("<",round(svsort[fit$n[2],i],2),"(",fit$n[2],")")),
                      #????pvalue??HR??95% CI
                      #̫С??p value??Ϊp < 0.001
                      pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                                 paste("p = ",round(p.val,3), sep = "")),
                                   HR, CI, sep = "\n"))
  ggsave(paste0(i,".PCPG.pdf"),width = 4,height = 4)}
length(pl)
F23=t(ptable)#write.csv(ptable, "pvalue_output_BLCA.csv", quote = F, row.names = F)


#---------------------------------------------PRAD------------------------------------
svdata <- read.table("3score.time.PRAD.txt",sep="\t",header = T,row.names = 1)
dim(svdata)
pct = 0.1
expr <- svdata[,3:ncol(svdata)]
#expr <- log2(expr + 1)
colnames(expr) <- make.names(colnames(expr))
svdata <- cbind.data.frame(svdata[,1:2],expr)
res.cut <- surv_cutpoint(svdata, time = "futime", 
                         event = "fustat", 
                         variables = names(svdata)[3:ncol(svdata)], 
                         minprop = 0.3) #Ĭ??????sample???ܵ???30%
res.cat <- surv_categorize(res.cut)
res.cat$futime=res.cat$futime/30   #####ʱ?�?????????
my.surv <- Surv(res.cat$futime, res.cat$fustat)
pl <- pl2 <- list() # ??ʼ?????????ߺ??¼??????????ߵ??б?
ptable <- data.frame(matrix(NA,ncol(svdata)-2,5))
colnames(ptable) <- c("ID", "pvalue", "HR", "CIlow", "CIup")
n = 0
for (i in colnames(res.cat)[3:ncol(svdata)]) {
  group <- res.cat[,i] 
  survival_dat <- data.frame(group = group)
  fit <- survfit(my.surv ~ group)
  group <- factor(group, levels = c("low", "high"))
  data.survdiff <- survdiff(my.surv ~ group)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  CancerType="PRAD"
  n = n + 1
  ptable$ID[n] <- i
  ptable$pvalue[n] <- p.val
  ptable$HR[n] <- HR
  ptable$CIlow[n] <- low95
  ptable$CIup[n] <- up95
  ptable$CancerType[n] <- CancerType
  HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
  CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
  svsort <- svdata[order(svdata[,i]),]
  pl[[i]]<-ggsurvplot(fit, data = survival_dat ,
                      #ggtheme = theme_bw(), #??Ҫ??????????????
                      conf.int = T, #???????????????????????????Ͱ?F?ĳ?T
                      #conf.int.style = "step",#?????????????ͣ????ɸ?Ϊribbon
                      censor = F, #????ʾ?۲?ֵ???ڵ?λ??
                      palette = c("#e64b35","#1b90be"), #?ߵ???ɫ??Ӧ?ߡ???
                      xlab="Time(Months)",
                      legend = c(0.8, 0.9),
                      legend.title = i,#??????д??ͼ????Ŀ??λ??
                      font.legend = 11,#ͼ??????????С
                      #font.title = 12,font.x = 10,font.y = 10,#??????????????С
                      #??ͼ???ϱ????ߵͷֽ????ı???��????????sample??��
                      legend.labs=c(paste0(">",round(svsort[fit$n[2],i],2),"(",fit$n[1],")"),
                                    paste0("<",round(svsort[fit$n[2],i],2),"(",fit$n[2],")")),
                      #????pvalue??HR??95% CI
                      #̫С??p value??Ϊp < 0.001
                      pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                                 paste("p = ",round(p.val,3), sep = "")),
                                   HR, CI, sep = "\n"))
  ggsave(paste0(i,".PRAD.pdf"),width = 4,height = 4)}
length(pl)
F24=t(ptable)#write.csv(ptable, "pvalue_output_BLCA.csv", quote = F, row.names = F)


#---------------------------------------------READ------------------------------------
svdata <- read.table("3score.time.READ.txt",sep="\t",header = T,row.names = 1)
dim(svdata)
pct = 0.1
expr <- svdata[,3:ncol(svdata)]
#expr <- log2(expr + 1)
colnames(expr) <- make.names(colnames(expr))
svdata <- cbind.data.frame(svdata[,1:2],expr)
res.cut <- surv_cutpoint(svdata, time = "futime", 
                         event = "fustat", 
                         variables = names(svdata)[3:ncol(svdata)], 
                         minprop = 0.3) #Ĭ??????sample???ܵ???30%
res.cat <- surv_categorize(res.cut)
res.cat$futime=res.cat$futime/30   #####ʱ?�?????????
my.surv <- Surv(res.cat$futime, res.cat$fustat)
pl <- pl2 <- list() # ??ʼ?????????ߺ??¼??????????ߵ??б?
ptable <- data.frame(matrix(NA,ncol(svdata)-2,5))
colnames(ptable) <- c("ID", "pvalue", "HR", "CIlow", "CIup")
n = 0
for (i in colnames(res.cat)[3:ncol(svdata)]) {
  group <- res.cat[,i] 
  survival_dat <- data.frame(group = group)
  fit <- survfit(my.surv ~ group)
  group <- factor(group, levels = c("low", "high"))
  data.survdiff <- survdiff(my.surv ~ group)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  CancerType="READ"
  n = n + 1
  ptable$ID[n] <- i
  ptable$pvalue[n] <- p.val
  ptable$HR[n] <- HR
  ptable$CIlow[n] <- low95
  ptable$CIup[n] <- up95
  ptable$CancerType[n] <- CancerType
  HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
  CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
  svsort <- svdata[order(svdata[,i]),]
  pl[[i]]<-ggsurvplot(fit, data = survival_dat ,
                      #ggtheme = theme_bw(), #??Ҫ??????????????
                      conf.int = T, #???????????????????????????Ͱ?F?ĳ?T
                      #conf.int.style = "step",#?????????????ͣ????ɸ?Ϊribbon
                      censor = F, #????ʾ?۲?ֵ???ڵ?λ??
                      palette = c("#e64b35","#1b90be"), #?ߵ???ɫ??Ӧ?ߡ???
                      xlab="Time(Months)",
                      legend = c(0.8, 0.9),
                      legend.title = i,#??????д??ͼ????Ŀ??λ??
                      font.legend = 11,#ͼ??????????С
                      #font.title = 12,font.x = 10,font.y = 10,#??????????????С
                      #??ͼ???ϱ????ߵͷֽ????ı???��????????sample??��
                      legend.labs=c(paste0(">",round(svsort[fit$n[2],i],2),"(",fit$n[1],")"),
                                    paste0("<",round(svsort[fit$n[2],i],2),"(",fit$n[2],")")),
                      #????pvalue??HR??95% CI
                      #̫С??p value??Ϊp < 0.001
                      pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                                 paste("p = ",round(p.val,3), sep = "")),
                                   HR, CI, sep = "\n"))
  ggsave(paste0(i,".READ.pdf"),width = 4,height = 4)}
length(pl)
F25=t(ptable)#write.csv(ptable, "pvalue_output_BLCA.csv", quote = F, row.names = F)


#---------------------------------------------SARC------------------------------------
svdata <- read.table("3score.time.SARC.txt",sep="\t",header = T,row.names = 1)
dim(svdata)
pct = 0.1
expr <- svdata[,3:ncol(svdata)]
#expr <- log2(expr + 1)
colnames(expr) <- make.names(colnames(expr))
svdata <- cbind.data.frame(svdata[,1:2],expr)
res.cut <- surv_cutpoint(svdata, time = "futime", 
                         event = "fustat", 
                         variables = names(svdata)[3:ncol(svdata)], 
                         minprop = 0.3) #Ĭ??????sample???ܵ???30%
res.cat <- surv_categorize(res.cut)
res.cat$futime=res.cat$futime/30   #####ʱ?�?????????
my.surv <- Surv(res.cat$futime, res.cat$fustat)
pl <- pl2 <- list() # ??ʼ?????????ߺ??¼??????????ߵ??б?
ptable <- data.frame(matrix(NA,ncol(svdata)-2,5))
colnames(ptable) <- c("ID", "pvalue", "HR", "CIlow", "CIup")
n = 0
for (i in colnames(res.cat)[3:ncol(svdata)]) {
  group <- res.cat[,i] 
  survival_dat <- data.frame(group = group)
  fit <- survfit(my.surv ~ group)
  group <- factor(group, levels = c("low", "high"))
  data.survdiff <- survdiff(my.surv ~ group)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  CancerType="SARC"
  n = n + 1
  ptable$ID[n] <- i
  ptable$pvalue[n] <- p.val
  ptable$HR[n] <- HR
  ptable$CIlow[n] <- low95
  ptable$CIup[n] <- up95
  ptable$CancerType[n] <- CancerType
  HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
  CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
  svsort <- svdata[order(svdata[,i]),]
  pl[[i]]<-ggsurvplot(fit, data = survival_dat ,
                      #ggtheme = theme_bw(), #??Ҫ??????????????
                      conf.int = T, #???????????????????????????Ͱ?F?ĳ?T
                      #conf.int.style = "step",#?????????????ͣ????ɸ?Ϊribbon
                      censor = F, #????ʾ?۲?ֵ???ڵ?λ??
                      palette = c("#e64b35","#1b90be"), #?ߵ???ɫ??Ӧ?ߡ???
                      xlab="Time(Months)",
                      legend = c(0.8, 0.9),
                      legend.title = i,#??????д??ͼ????Ŀ??λ??
                      font.legend = 11,#ͼ??????????С
                      #font.title = 12,font.x = 10,font.y = 10,#??????????????С
                      #??ͼ???ϱ????ߵͷֽ????ı???��????????sample??��
                      legend.labs=c(paste0(">",round(svsort[fit$n[2],i],2),"(",fit$n[1],")"),
                                    paste0("<",round(svsort[fit$n[2],i],2),"(",fit$n[2],")")),
                      #????pvalue??HR??95% CI
                      #̫С??p value??Ϊp < 0.001
                      pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                                 paste("p = ",round(p.val,3), sep = "")),
                                   HR, CI, sep = "\n"))
  ggsave(paste0(i,".SARC.pdf"),width = 4,height = 4)}
length(pl)
F26=t(ptable)#write.csv(ptable, "pvalue_output_BLCA.csv", quote = F, row.names = F)


#---------------------------------------------SKCM------------------------------------
svdata <- read.table("3score.time.SKCM.txt",sep="\t",header = T,row.names = 1)
dim(svdata)
pct = 0.1
expr <- svdata[,3:ncol(svdata)]
#expr <- log2(expr + 1)
colnames(expr) <- make.names(colnames(expr))
svdata <- cbind.data.frame(svdata[,1:2],expr)
res.cut <- surv_cutpoint(svdata, time = "futime", 
                         event = "fustat", 
                         variables = names(svdata)[3:ncol(svdata)], 
                         minprop = 0.3) #Ĭ??????sample???ܵ???30%
res.cat <- surv_categorize(res.cut)
res.cat$futime=res.cat$futime/30   #####ʱ?�?????????
my.surv <- Surv(res.cat$futime, res.cat$fustat)
pl <- pl2 <- list() # ??ʼ?????????ߺ??¼??????????ߵ??б?
ptable <- data.frame(matrix(NA,ncol(svdata)-2,5))
colnames(ptable) <- c("ID", "pvalue", "HR", "CIlow", "CIup")
n = 0
for (i in colnames(res.cat)[3:ncol(svdata)]) {
  group <- res.cat[,i] 
  survival_dat <- data.frame(group = group)
  fit <- survfit(my.surv ~ group)
  group <- factor(group, levels = c("low", "high"))
  data.survdiff <- survdiff(my.surv ~ group)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  CancerType="SKCM"
  n = n + 1
  ptable$ID[n] <- i
  ptable$pvalue[n] <- p.val
  ptable$HR[n] <- HR
  ptable$CIlow[n] <- low95
  ptable$CIup[n] <- up95
  ptable$CancerType[n] <- CancerType
  HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
  CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
  svsort <- svdata[order(svdata[,i]),]
  pl[[i]]<-ggsurvplot(fit, data = survival_dat ,
                      #ggtheme = theme_bw(), #??Ҫ??????????????
                      conf.int = T, #???????????????????????????Ͱ?F?ĳ?T
                      #conf.int.style = "step",#?????????????ͣ????ɸ?Ϊribbon
                      censor = F, #????ʾ?۲?ֵ???ڵ?λ??
                      palette = c("#e64b35","#1b90be"), #?ߵ???ɫ??Ӧ?ߡ???
                      xlab="Time(Months)",
                      legend = c(0.8, 0.9),
                      legend.title = i,#??????д??ͼ????Ŀ??λ??
                      font.legend = 11,#ͼ??????????С
                      #font.title = 12,font.x = 10,font.y = 10,#??????????????С
                      #??ͼ???ϱ????ߵͷֽ????ı???��????????sample??��
                      legend.labs=c(paste0(">",round(svsort[fit$n[2],i],2),"(",fit$n[1],")"),
                                    paste0("<",round(svsort[fit$n[2],i],2),"(",fit$n[2],")")),
                      #????pvalue??HR??95% CI
                      #̫С??p value??Ϊp < 0.001
                      pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                                 paste("p = ",round(p.val,3), sep = "")),
                                   HR, CI, sep = "\n"))
  ggsave(paste0(i,".SKCM.pdf"),width = 4,height = 4)}
length(pl)
F27=t(ptable)#write.csv(ptable, "pvalue_output_BLCA.csv", quote = F, row.names = F)


#---------------------------------------------STAD------------------------------------
svdata <- read.table("3score.time.STAD.txt",sep="\t",header = T,row.names = 1)
dim(svdata)
pct = 0.1
expr <- svdata[,3:ncol(svdata)]
#expr <- log2(expr + 1)
colnames(expr) <- make.names(colnames(expr))
svdata <- cbind.data.frame(svdata[,1:2],expr)
res.cut <- surv_cutpoint(svdata, time = "futime", 
                         event = "fustat", 
                         variables = names(svdata)[3:ncol(svdata)], 
                         minprop = 0.3) #Ĭ??????sample???ܵ???30%
res.cat <- surv_categorize(res.cut)
res.cat$futime=res.cat$futime/30   #####ʱ?�?????????
my.surv <- Surv(res.cat$futime, res.cat$fustat)
pl <- pl2 <- list() # ??ʼ?????????ߺ??¼??????????ߵ??б?
ptable <- data.frame(matrix(NA,ncol(svdata)-2,5))
colnames(ptable) <- c("ID", "pvalue", "HR", "CIlow", "CIup")
n = 0
for (i in colnames(res.cat)[3:ncol(svdata)]) {
  group <- res.cat[,i] 
  survival_dat <- data.frame(group = group)
  fit <- survfit(my.surv ~ group)
  group <- factor(group, levels = c("low", "high"))
  data.survdiff <- survdiff(my.surv ~ group)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  CancerType="STAD"
  n = n + 1
  ptable$ID[n] <- i
  ptable$pvalue[n] <- p.val
  ptable$HR[n] <- HR
  ptable$CIlow[n] <- low95
  ptable$CIup[n] <- up95
  ptable$CancerType[n] <- CancerType
  HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
  CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
  svsort <- svdata[order(svdata[,i]),]
  pl[[i]]<-ggsurvplot(fit, data = survival_dat ,
                      #ggtheme = theme_bw(), #??Ҫ??????????????
                      conf.int = T, #???????????????????????????Ͱ?F?ĳ?T
                      #conf.int.style = "step",#?????????????ͣ????ɸ?Ϊribbon
                      censor = F, #????ʾ?۲?ֵ???ڵ?λ??
                      palette = c("#e64b35","#1b90be"), #?ߵ???ɫ??Ӧ?ߡ???
                      xlab="Time(Months)",
                      legend = c(0.8, 0.9),
                      legend.title = i,#??????д??ͼ????Ŀ??λ??
                      font.legend = 11,#ͼ??????????С
                      #font.title = 12,font.x = 10,font.y = 10,#??????????????С
                      #??ͼ???ϱ????ߵͷֽ????ı???��????????sample??��
                      legend.labs=c(paste0(">",round(svsort[fit$n[2],i],2),"(",fit$n[1],")"),
                                    paste0("<",round(svsort[fit$n[2],i],2),"(",fit$n[2],")")),
                      #????pvalue??HR??95% CI
                      #̫С??p value??Ϊp < 0.001
                      pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                                 paste("p = ",round(p.val,3), sep = "")),
                                   HR, CI, sep = "\n"))
  ggsave(paste0(i,".STAD.pdf"),width = 4,height = 4)}
length(pl)
F28=t(ptable)#write.csv(ptable, "pvalue_output_BLCA.csv", quote = F, row.names = F)


#---------------------------------------------TGCT------------------------------------
svdata <- read.table("3score.time.TGCT.txt",sep="\t",header = T,row.names = 1)
dim(svdata)
pct = 0.1
expr <- svdata[,3:ncol(svdata)]
#expr <- log2(expr + 1)
colnames(expr) <- make.names(colnames(expr))
svdata <- cbind.data.frame(svdata[,1:2],expr)
res.cut <- surv_cutpoint(svdata, time = "futime", 
                         event = "fustat", 
                         variables = names(svdata)[3:ncol(svdata)], 
                         minprop = 0.3) #Ĭ??????sample???ܵ???30%
res.cat <- surv_categorize(res.cut)
res.cat$futime=res.cat$futime/30   #####ʱ?�?????????
my.surv <- Surv(res.cat$futime, res.cat$fustat)
pl <- pl2 <- list() # ??ʼ?????????ߺ??¼??????????ߵ??б?
ptable <- data.frame(matrix(NA,ncol(svdata)-2,5))
colnames(ptable) <- c("ID", "pvalue", "HR", "CIlow", "CIup")
n = 0
for (i in colnames(res.cat)[3:ncol(svdata)]) {
  group <- res.cat[,i] 
  survival_dat <- data.frame(group = group)
  fit <- survfit(my.surv ~ group)
  group <- factor(group, levels = c("low", "high"))
  data.survdiff <- survdiff(my.surv ~ group)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  CancerType="TGCT"
  n = n + 1
  ptable$ID[n] <- i
  ptable$pvalue[n] <- p.val
  ptable$HR[n] <- HR
  ptable$CIlow[n] <- low95
  ptable$CIup[n] <- up95
  ptable$CancerType[n] <- CancerType
  HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
  CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
  svsort <- svdata[order(svdata[,i]),]
  pl[[i]]<-ggsurvplot(fit, data = survival_dat ,
                      #ggtheme = theme_bw(), #??Ҫ??????????????
                      conf.int = T, #???????????????????????????Ͱ?F?ĳ?T
                      #conf.int.style = "step",#?????????????ͣ????ɸ?Ϊribbon
                      censor = F, #????ʾ?۲?ֵ???ڵ?λ??
                      palette = c("#e64b35","#1b90be"), #?ߵ???ɫ??Ӧ?ߡ???
                      xlab="Time(Months)",
                      legend = c(0.8, 0.9),
                      legend.title = i,#??????д??ͼ????Ŀ??λ??
                      font.legend = 11,#ͼ??????????С
                      #font.title = 12,font.x = 10,font.y = 10,#??????????????С
                      #??ͼ???ϱ????ߵͷֽ????ı???��????????sample??��
                      legend.labs=c(paste0(">",round(svsort[fit$n[2],i],2),"(",fit$n[1],")"),
                                    paste0("<",round(svsort[fit$n[2],i],2),"(",fit$n[2],")")),
                      #????pvalue??HR??95% CI
                      #̫С??p value??Ϊp < 0.001
                      pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                                 paste("p = ",round(p.val,3), sep = "")),
                                   HR, CI, sep = "\n"))
  ggsave(paste0(i,".TGCT.pdf"),width = 4,height = 4)}
length(pl)
F29=t(ptable)#write.csv(ptable, "pvalue_output_BLCA.csv", quote = F, row.names = F)


#---------------------------------------------THCA------------------------------------
svdata <- read.table("3score.time.THCA.txt",sep="\t",header = T,row.names = 1)
dim(svdata)
pct = 0.1
expr <- svdata[,3:ncol(svdata)]
#expr <- log2(expr + 1)
colnames(expr) <- make.names(colnames(expr))
svdata <- cbind.data.frame(svdata[,1:2],expr)
res.cut <- surv_cutpoint(svdata, time = "futime", 
                         event = "fustat", 
                         variables = names(svdata)[3:ncol(svdata)], 
                         minprop = 0.3) #Ĭ??????sample???ܵ???30%
res.cat <- surv_categorize(res.cut)
res.cat$futime=res.cat$futime/30   #####ʱ?�?????????
my.surv <- Surv(res.cat$futime, res.cat$fustat)
pl <- pl2 <- list() # ??ʼ?????????ߺ??¼??????????ߵ??б?
ptable <- data.frame(matrix(NA,ncol(svdata)-2,5))
colnames(ptable) <- c("ID", "pvalue", "HR", "CIlow", "CIup")
n = 0
for (i in colnames(res.cat)[3:ncol(svdata)]) {
  group <- res.cat[,i] 
  survival_dat <- data.frame(group = group)
  fit <- survfit(my.surv ~ group)
  group <- factor(group, levels = c("low", "high"))
  data.survdiff <- survdiff(my.surv ~ group)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  CancerType="THCA"
  n = n + 1
  ptable$ID[n] <- i
  ptable$pvalue[n] <- p.val
  ptable$HR[n] <- HR
  ptable$CIlow[n] <- low95
  ptable$CIup[n] <- up95
  ptable$CancerType[n] <- CancerType
  HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
  CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
  svsort <- svdata[order(svdata[,i]),]
  pl[[i]]<-ggsurvplot(fit, data = survival_dat ,
                      #ggtheme = theme_bw(), #??Ҫ??????????????
                      conf.int = T, #???????????????????????????Ͱ?F?ĳ?T
                      #conf.int.style = "step",#?????????????ͣ????ɸ?Ϊribbon
                      censor = F, #????ʾ?۲?ֵ???ڵ?λ??
                      palette = c("#e64b35","#1b90be"), #?ߵ???ɫ??Ӧ?ߡ???
                      xlab="Time(Months)",
                      legend = c(0.8, 0.9),
                      legend.title = i,#??????д??ͼ????Ŀ??λ??
                      font.legend = 11,#ͼ??????????С
                      #font.title = 12,font.x = 10,font.y = 10,#??????????????С
                      #??ͼ???ϱ????ߵͷֽ????ı???��????????sample??��
                      legend.labs=c(paste0(">",round(svsort[fit$n[2],i],2),"(",fit$n[1],")"),
                                    paste0("<",round(svsort[fit$n[2],i],2),"(",fit$n[2],")")),
                      #????pvalue??HR??95% CI
                      #̫С??p value??Ϊp < 0.001
                      pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                                 paste("p = ",round(p.val,3), sep = "")),
                                   HR, CI, sep = "\n"))
  ggsave(paste0(i,".THCA.pdf"),width = 4,height = 4)}
length(pl)
F30=t(ptable)#write.csv(ptable, "pvalue_output_BLCA.csv", quote = F, row.names = F)


#---------------------------------------------THYM------------------------------------
svdata <- read.table("3score.time.THYM.txt",sep="\t",header = T,row.names = 1)
dim(svdata)
pct = 0.1
expr <- svdata[,3:ncol(svdata)]
#expr <- log2(expr + 1)
colnames(expr) <- make.names(colnames(expr))
svdata <- cbind.data.frame(svdata[,1:2],expr)
res.cut <- surv_cutpoint(svdata, time = "futime", 
                         event = "fustat", 
                         variables = names(svdata)[3:ncol(svdata)], 
                         minprop = 0.3) #Ĭ??????sample???ܵ???30%
res.cat <- surv_categorize(res.cut)
res.cat$futime=res.cat$futime/30   #####ʱ?�?????????
my.surv <- Surv(res.cat$futime, res.cat$fustat)
pl <- pl2 <- list() # ??ʼ?????????ߺ??¼??????????ߵ??б?
ptable <- data.frame(matrix(NA,ncol(svdata)-2,5))
colnames(ptable) <- c("ID", "pvalue", "HR", "CIlow", "CIup")
n = 0
for (i in colnames(res.cat)[3:ncol(svdata)]) {
  group <- res.cat[,i] 
  survival_dat <- data.frame(group = group)
  fit <- survfit(my.surv ~ group)
  group <- factor(group, levels = c("low", "high"))
  data.survdiff <- survdiff(my.surv ~ group)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  CancerType="THYM"
  n = n + 1
  ptable$ID[n] <- i
  ptable$pvalue[n] <- p.val
  ptable$HR[n] <- HR
  ptable$CIlow[n] <- low95
  ptable$CIup[n] <- up95
  ptable$CancerType[n] <- CancerType
  HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
  CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
  svsort <- svdata[order(svdata[,i]),]
  pl[[i]]<-ggsurvplot(fit, data = survival_dat ,
                      #ggtheme = theme_bw(), #??Ҫ??????????????
                      conf.int = T, #???????????????????????????Ͱ?F?ĳ?T
                      #conf.int.style = "step",#?????????????ͣ????ɸ?Ϊribbon
                      censor = F, #????ʾ?۲?ֵ???ڵ?λ??
                      palette = c("#e64b35","#1b90be"), #?ߵ???ɫ??Ӧ?ߡ???
                      xlab="Time(Months)",
                      legend = c(0.8, 0.9),
                      legend.title = i,#??????д??ͼ????Ŀ??λ??
                      font.legend = 11,#ͼ??????????С
                      #font.title = 12,font.x = 10,font.y = 10,#??????????????С
                      #??ͼ???ϱ????ߵͷֽ????ı???��????????sample??��
                      legend.labs=c(paste0(">",round(svsort[fit$n[2],i],2),"(",fit$n[1],")"),
                                    paste0("<",round(svsort[fit$n[2],i],2),"(",fit$n[2],")")),
                      #????pvalue??HR??95% CI
                      #̫С??p value??Ϊp < 0.001
                      pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                                 paste("p = ",round(p.val,3), sep = "")),
                                   HR, CI, sep = "\n"))
  ggsave(paste0(i,".THYM.pdf"),width = 4,height = 4)}
length(pl)
F31=t(ptable)#write.csv(ptable, "pvalue_output_BLCA.csv", quote = F, row.names = F)


#---------------------------------------------UCEC------------------------------------
svdata <- read.table("3score.time.UCEC.txt",sep="\t",header = T,row.names = 1)
dim(svdata)
pct = 0.1
expr <- svdata[,3:ncol(svdata)]
#expr <- log2(expr + 1)
colnames(expr) <- make.names(colnames(expr))
svdata <- cbind.data.frame(svdata[,1:2],expr)
res.cut <- surv_cutpoint(svdata, time = "futime", 
                         event = "fustat", 
                         variables = names(svdata)[3:ncol(svdata)], 
                         minprop = 0.3) #Ĭ??????sample???ܵ???30%
res.cat <- surv_categorize(res.cut)
res.cat$futime=res.cat$futime/30   #####ʱ?�?????????
my.surv <- Surv(res.cat$futime, res.cat$fustat)
pl <- pl2 <- list() # ??ʼ?????????ߺ??¼??????????ߵ??б?
ptable <- data.frame(matrix(NA,ncol(svdata)-2,5))
colnames(ptable) <- c("ID", "pvalue", "HR", "CIlow", "CIup")
n = 0
for (i in colnames(res.cat)[3:ncol(svdata)]) {
  group <- res.cat[,i] 
  survival_dat <- data.frame(group = group)
  fit <- survfit(my.surv ~ group)
  group <- factor(group, levels = c("low", "high"))
  data.survdiff <- survdiff(my.surv ~ group)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  CancerType="UCEC"
  n = n + 1
  ptable$ID[n] <- i
  ptable$pvalue[n] <- p.val
  ptable$HR[n] <- HR
  ptable$CIlow[n] <- low95
  ptable$CIup[n] <- up95
  ptable$CancerType[n] <- CancerType
  HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
  CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
  svsort <- svdata[order(svdata[,i]),]
  pl[[i]]<-ggsurvplot(fit, data = survival_dat ,
                      #ggtheme = theme_bw(), #??Ҫ??????????????
                      conf.int = T, #???????????????????????????Ͱ?F?ĳ?T
                      #conf.int.style = "step",#?????????????ͣ????ɸ?Ϊribbon
                      censor = F, #????ʾ?۲?ֵ???ڵ?λ??
                      palette = c("#e64b35","#1b90be"), #?ߵ???ɫ??Ӧ?ߡ???
                      xlab="Time(Months)",
                      legend = c(0.8, 0.9),
                      legend.title = i,#??????д??ͼ????Ŀ??λ??
                      font.legend = 11,#ͼ??????????С
                      #font.title = 12,font.x = 10,font.y = 10,#??????????????С
                      #??ͼ???ϱ????ߵͷֽ????ı???��????????sample??��
                      legend.labs=c(paste0(">",round(svsort[fit$n[2],i],2),"(",fit$n[1],")"),
                                    paste0("<",round(svsort[fit$n[2],i],2),"(",fit$n[2],")")),
                      #????pvalue??HR??95% CI
                      #̫С??p value??Ϊp < 0.001
                      pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                                 paste("p = ",round(p.val,3), sep = "")),
                                   HR, CI, sep = "\n"))
  ggsave(paste0(i,".UCEC.pdf"),width = 4,height = 4)}
length(pl)
F32=t(ptable)#write.csv(ptable, "pvalue_output_BLCA.csv", quote = F, row.names = F)


#---------------------------------------------UCS------------------------------------
svdata <- read.table("3score.time.UCS.txt",sep="\t",header = T,row.names = 1)
dim(svdata)
pct = 0.1
expr <- svdata[,3:ncol(svdata)]
#expr <- log2(expr + 1)
colnames(expr) <- make.names(colnames(expr))
svdata <- cbind.data.frame(svdata[,1:2],expr)
res.cut <- surv_cutpoint(svdata, time = "futime", 
                         event = "fustat", 
                         variables = names(svdata)[3:ncol(svdata)], 
                         minprop = 0.3) #Ĭ??????sample???ܵ???30%
res.cat <- surv_categorize(res.cut)
res.cat$futime=res.cat$futime/30   #####ʱ?�?????????
my.surv <- Surv(res.cat$futime, res.cat$fustat)
pl <- pl2 <- list() # ??ʼ?????????ߺ??¼??????????ߵ??б?
ptable <- data.frame(matrix(NA,ncol(svdata)-2,5))
colnames(ptable) <- c("ID", "pvalue", "HR", "CIlow", "CIup")
n = 0
for (i in colnames(res.cat)[3:ncol(svdata)]) {
  group <- res.cat[,i] 
  survival_dat <- data.frame(group = group)
  fit <- survfit(my.surv ~ group)
  group <- factor(group, levels = c("low", "high"))
  data.survdiff <- survdiff(my.surv ~ group)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  CancerType="UCS"
  n = n + 1
  ptable$ID[n] <- i
  ptable$pvalue[n] <- p.val
  ptable$HR[n] <- HR
  ptable$CIlow[n] <- low95
  ptable$CIup[n] <- up95
  ptable$CancerType[n] <- CancerType
  HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
  CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
  svsort <- svdata[order(svdata[,i]),]
  pl[[i]]<-ggsurvplot(fit, data = survival_dat ,
                      #ggtheme = theme_bw(), #??Ҫ??????????????
                      conf.int = T, #???????????????????????????Ͱ?F?ĳ?T
                      #conf.int.style = "step",#?????????????ͣ????ɸ?Ϊribbon
                      censor = F, #????ʾ?۲?ֵ???ڵ?λ??
                      palette = c("#e64b35","#1b90be"), #?ߵ???ɫ??Ӧ?ߡ???
                      xlab="Time(Months)",
                      legend = c(0.8, 0.9),
                      legend.title = i,#??????д??ͼ????Ŀ??λ??
                      font.legend = 11,#ͼ??????????С
                      #font.title = 12,font.x = 10,font.y = 10,#??????????????С
                      #??ͼ???ϱ????ߵͷֽ????ı???��????????sample??��
                      legend.labs=c(paste0(">",round(svsort[fit$n[2],i],2),"(",fit$n[1],")"),
                                    paste0("<",round(svsort[fit$n[2],i],2),"(",fit$n[2],")")),
                      #????pvalue??HR??95% CI
                      #̫С??p value??Ϊp < 0.001
                      pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                                 paste("p = ",round(p.val,3), sep = "")),
                                   HR, CI, sep = "\n"))
  ggsave(paste0(i,".UCS.pdf"),width = 4,height = 4)}
length(pl)
F33=t(ptable)#write.csv(ptable, "pvalue_output_BLCA.csv", quote = F, row.names = F)


#---------------------------------------------UVM------------------------------------
svdata <- read.table("3score.time.UVM.txt",sep="\t",header = T,row.names = 1)
dim(svdata)
pct = 0.1
expr <- svdata[,3:ncol(svdata)]
#expr <- log2(expr + 1)
colnames(expr) <- make.names(colnames(expr))
svdata <- cbind.data.frame(svdata[,1:2],expr)
res.cut <- surv_cutpoint(svdata, time = "futime", 
                         event = "fustat", 
                         variables = names(svdata)[3:ncol(svdata)], 
                         minprop = 0.3) #Ĭ??????sample???ܵ???30%
res.cat <- surv_categorize(res.cut)
res.cat$futime=res.cat$futime/30   #####ʱ?�?????????
my.surv <- Surv(res.cat$futime, res.cat$fustat)
pl <- pl2 <- list() # ??ʼ?????????ߺ??¼??????????ߵ??б?
ptable <- data.frame(matrix(NA,ncol(svdata)-2,5))
colnames(ptable) <- c("ID", "pvalue", "HR", "CIlow", "CIup")
n = 0
for (i in colnames(res.cat)[3:ncol(svdata)]) {
  group <- res.cat[,i] 
  survival_dat <- data.frame(group = group)
  fit <- survfit(my.surv ~ group)
  group <- factor(group, levels = c("low", "high"))
  data.survdiff <- survdiff(my.surv ~ group)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  CancerType="UVM"
  n = n + 1
  ptable$ID[n] <- i
  ptable$pvalue[n] <- p.val
  ptable$HR[n] <- HR
  ptable$CIlow[n] <- low95
  ptable$CIup[n] <- up95
  ptable$CancerType[n] <- CancerType
  HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
  CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
  svsort <- svdata[order(svdata[,i]),]
  pl[[i]]<-ggsurvplot(fit, data = survival_dat ,
                      #ggtheme = theme_bw(), #??Ҫ??????????????
                      conf.int = T, #???????????????????????????Ͱ?F?ĳ?T
                      #conf.int.style = "step",#?????????????ͣ????ɸ?Ϊribbon
                      censor = F, #????ʾ?۲?ֵ???ڵ?λ??
                      palette = c("#e64b35","#1b90be"), #?ߵ???ɫ??Ӧ?ߡ???
                      xlab="Time(Months)",
                      legend = c(0.8, 0.9),
                      legend.title = i,#??????д??ͼ????Ŀ??λ??
                      font.legend = 11,#ͼ??????????С
                      #font.title = 12,font.x = 10,font.y = 10,#??????????????С
                      #??ͼ???ϱ????ߵͷֽ????ı???��????????sample??��
                      legend.labs=c(paste0(">",round(svsort[fit$n[2],i],2),"(",fit$n[1],")"),
                                    paste0("<",round(svsort[fit$n[2],i],2),"(",fit$n[2],")")),
                      #????pvalue??HR??95% CI
                      #̫С??p value??Ϊp < 0.001
                      pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                                 paste("p = ",round(p.val,3), sep = "")),
                                   HR, CI, sep = "\n"))
  ggsave(paste0(i,".UVM.pdf"),width = 4,height = 4)}
length(pl)
F01=t(ptable)#write.csv(ptable, "pvalue_output_BLCA.csv", quote = F, row.names = F)



ptable.all=t(cbind(F01,
                   F02,
                   F03,
                   F04,
                   F05,
                   F06,
                   F07,
                   F08,
                   F09,
                   F10,
                   F11,
                   F12,
                   F13,
                   F14,
                   F15,
                   F16,
                   F17,
                   F18,
                   F19,
                   F20,
                   F21,
                   F22,
                   F23,
                   F24,
                   F25,
                   F26,
                   F27,
                   F28,
                   F29,
                   F30,
                   F31,
                   F32,
                   F33
))
write.table(ptable.all,file="ptable.pan-cancer.3score.txt",sep="\t",row.names=F,quote=F)