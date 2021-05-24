#setwd("E:/MyProject_YanLab/Project4_CTCFbinding/logistics/Logistic_Methyl/windowsize1000.avg")
#setwd("E:/MyProject_YanLab/Project4_CTCFbinding/logistics/Logistic_Methyl/windowsize1000.max")
#setwd("E:/MyProject_YanLab/Project4_CTCFbinding/logistics/Logistic_Methyl/windowsize500.avg")
#setwd("E:/MyProject_YanLab/Project4_CTCFbinding/logistics/Logistic_Methyl/windowsize500.max")
#setwd("E:/MyProject_YanLab/Project4_CTCFbinding/logistics/Logistic_Methyl/windowsize300.avg")
#setwd("E:/MyProject_YanLab/Project4_CTCFbinding/logistics/Logistic_Methyl/windowsize300.max")

library(tidyverse)
library(tidyr)
library("ggplot2")
library("easyGgplot2")
library(pROC)

#pdf("windowsize1000.avg.pdf",width=10, height=10)
#pdf("windowsize1000.max.pdf",width=10, height=10)
#pdf("windowsize500.avg.pdf",width=10, height=10)
#pdf("windowsize500.max.pdf",width=10, height=10)
#pdf("windowsize300.avg.pdf",width=10, height=10)
#pdf("windowsize300.max.pdf",width=10, height=10)

# We treated CTCF binding as Y (binding->yes; not binding->no), and mGC% data as X.
# read data from file
# chr1
mydata.GCH<-read.table("chr1.GCH.bed",sep="\t")
# set colnames
colnames(mydata.GCH)<-c("chrom", "start", "end", "Methyl", "Repli", "ChIP", "DNase")
# what we care about in the model is the result (CTCT binding=yes), so we just sorted this factor.
mydata.GCH$Methyl<-factor(mydata.GCH$Methyl,levels = c("no","yes"),order=TRUE)
# construct Logistic model
model<-glm(Methyl~ChIP+Repli+DNase,data=mydata.GCH,family = "binomial")
summary(model)
# the function of step realizes stepwise regression method, trace=0 is for not print information during the running of step
model2<-step(object = model,trace = 0)
summary(model2)
# It was found from the results that the P values of all variables were less than 0.05, and the relatively important variables were retained through the significance test
anova(object = model2,test = "Chisq")
# Significance test of the model
# It can be seen from the above that as the variables are added to the model one by one from the first to the last, the model finally passes the significance test, indicating that the model composed of the above variables is meaningful and correct.
prob<-predict(object =model2,newdata=mydata.GCH,type = "response")
# Prediction accuracy of the model for test data (here, we select train data as test although there maybe is overfitting-related problems)
pred<-ifelse(prob>=0.5,"yes","no")
pred<-factor(pred,levels = c("no","yes"),order=TRUE)
mydata.GCH<-cbind(mydata.GCH,pred)
f<-table(mydata.GCH$Methyl,pred)
# draw ROC curve and calculate AUC
roc_curve <- roc(mydata.GCH$Methyl,prob)
x <- 1-roc_curve$specificities
y <- roc_curve$sensitivities
p <- ggplot(data = NULL, mapping = aes(x= x, y = y))
p + geom_line(colour = 'red') +geom_abline(intercept = 0, slope = 1)+ annotate('text', x = 0.4, y = 0.5, label =paste('AUC=',round(roc_curve$auc,2)))+ labs(x = '1-specificities',y = 'sensitivities', title = 'Chr1: ROC Curve (Methyl~ChIP+Repli+DNase)')
write.table(mydata.GCH,"chr1.result.txt",quote=F,sep="\t",row.names=F)
# chr2
mydata.GCH<-read.table("chr2.GCH.bed",sep="\t")
# set colnames
colnames(mydata.GCH)<-c("chrom", "start", "end", "Methyl", "Repli", "ChIP", "DNase")
# what we care about in the model is the result (CTCT binding=yes), so we just sorted this factor.
mydata.GCH$Methyl<-factor(mydata.GCH$Methyl,levels = c("no","yes"),order=TRUE)
# construct Logistic model
model<-glm(Methyl~ChIP+Repli+DNase,data=mydata.GCH,family = "binomial")
summary(model)
# the function of step realizes stepwise regression method, trace=0 is for not print information during the running of step
model2<-step(object = model,trace = 0)
summary(model2)
# It was found from the results that the P values of all variables were less than 0.05, and the relatively important variables were retained through the significance test
anova(object = model2,test = "Chisq")
# Significance test of the model
# It can be seen from the above that as the variables are added to the model one by one from the first to the last, the model finally passes the significance test, indicating that the model composed of the above variables is meaningful and correct.
prob<-predict(object =model2,newdata=mydata.GCH,type = "response")
# Prediction accuracy of the model for test data (here, we select train data as test although there maybe is overfitting-related problems)
pred<-ifelse(prob>=0.5,"yes","no")
pred<-factor(pred,levels = c("no","yes"),order=TRUE)
mydata.GCH<-cbind(mydata.GCH,pred)
f<-table(mydata.GCH$Methyl,pred)
# draw ROC curve and calculate AUC
roc_curve <- roc(mydata.GCH$Methyl,prob)
x <- 1-roc_curve$specificities
y <- roc_curve$sensitivities
p <- ggplot(data = NULL, mapping = aes(x= x, y = y))
p + geom_line(colour = 'red') +geom_abline(intercept = 0, slope = 1)+ annotate('text', x = 0.4, y = 0.5, label =paste('AUC=',round(roc_curve$auc,2)))+ labs(x = '1-specificities',y = 'sensitivities', title = 'chr2: ROC Curve (Methyl~ChIP+Repli+DNase)')
write.table(mydata.GCH,"chr2.result.txt",quote=F,sep="\t",row.names=F)
# chr3
mydata.GCH<-read.table("chr3.GCH.bed",sep="\t")
# set colnames
colnames(mydata.GCH)<-c("chrom", "start", "end", "Methyl", "Repli", "ChIP", "DNase")
# what we care about in the model is the result (CTCT binding=yes), so we just sorted this factor.
mydata.GCH$Methyl<-factor(mydata.GCH$Methyl,levels = c("no","yes"),order=TRUE)
# construct Logistic model
model<-glm(Methyl~ChIP+Repli+DNase,data=mydata.GCH,family = "binomial")
summary(model)
# the function of step realizes stepwise regression method, trace=0 is for not print information during the running of step
model2<-step(object = model,trace = 0)
summary(model2)
# It was found from the results that the P values of all variables were less than 0.05, and the relatively important variables were retained through the significance test
anova(object = model2,test = "Chisq")
# Significance test of the model
# It can be seen from the above that as the variables are added to the model one by one from the first to the last, the model finally passes the significance test, indicating that the model composed of the above variables is meaningful and correct.
prob<-predict(object =model2,newdata=mydata.GCH,type = "response")
# Prediction accuracy of the model for test data (here, we select train data as test although there maybe is overfitting-related problems)
pred<-ifelse(prob>=0.5,"yes","no")
pred<-factor(pred,levels = c("no","yes"),order=TRUE)
mydata.GCH<-cbind(mydata.GCH,pred)
f<-table(mydata.GCH$Methyl,pred)
# draw ROC curve and calculate AUC
roc_curve <- roc(mydata.GCH$Methyl,prob)
x <- 1-roc_curve$specificities
y <- roc_curve$sensitivities
p <- ggplot(data = NULL, mapping = aes(x= x, y = y))
p + geom_line(colour = 'red') +geom_abline(intercept = 0, slope = 1)+ annotate('text', x = 0.4, y = 0.5, label =paste('AUC=',round(roc_curve$auc,2)))+ labs(x = '1-specificities',y = 'sensitivities', title = 'chr3: ROC Curve (Methyl~ChIP+Repli+DNase)')
write.table(mydata.GCH,"chr3.result.txt",quote=F,sep="\t",row.names=F)
# chr4
mydata.GCH<-read.table("chr4.GCH.bed",sep="\t")
# set colnames
colnames(mydata.GCH)<-c("chrom", "start", "end", "Methyl", "Repli", "ChIP", "DNase")
# what we care about in the model is the result (CTCT binding=yes), so we just sorted this factor.
mydata.GCH$Methyl<-factor(mydata.GCH$Methyl,levels = c("no","yes"),order=TRUE)
# construct Logistic model
model<-glm(Methyl~ChIP+Repli+DNase,data=mydata.GCH,family = "binomial")
summary(model)
# the function of step realizes stepwise regression method, trace=0 is for not print information during the running of step
model2<-step(object = model,trace = 0)
summary(model2)
# It was found from the results that the P values of all variables were less than 0.05, and the relatively important variables were retained through the significance test
anova(object = model2,test = "Chisq")
# Significance test of the model
# It can be seen from the above that as the variables are added to the model one by one from the first to the last, the model finally passes the significance test, indicating that the model composed of the above variables is meaningful and correct.
prob<-predict(object =model2,newdata=mydata.GCH,type = "response")
# Prediction accuracy of the model for test data (here, we select train data as test although there maybe is overfitting-related problems)
pred<-ifelse(prob>=0.5,"yes","no")
pred<-factor(pred,levels = c("no","yes"),order=TRUE)
mydata.GCH<-cbind(mydata.GCH,pred)
f<-table(mydata.GCH$Methyl,pred)
# draw ROC curve and calculate AUC
roc_curve <- roc(mydata.GCH$Methyl,prob)
x <- 1-roc_curve$specificities
y <- roc_curve$sensitivities
p <- ggplot(data = NULL, mapping = aes(x= x, y = y))
p + geom_line(colour = 'red') +geom_abline(intercept = 0, slope = 1)+ annotate('text', x = 0.4, y = 0.5, label =paste('AUC=',round(roc_curve$auc,2)))+ labs(x = '1-specificities',y = 'sensitivities', title = 'chr4: ROC Curve (Methyl~ChIP+Repli+DNase)')
write.table(mydata.GCH,"chr4.result.txt",quote=F,sep="\t",row.names=F)
# chr5
mydata.GCH<-read.table("chr5.GCH.bed",sep="\t")
# set colnames
colnames(mydata.GCH)<-c("chrom", "start", "end", "Methyl", "Repli", "ChIP", "DNase")
# what we care about in the model is the result (CTCT binding=yes), so we just sorted this factor.
mydata.GCH$Methyl<-factor(mydata.GCH$Methyl,levels = c("no","yes"),order=TRUE)
# construct Logistic model
model<-glm(Methyl~ChIP+Repli+DNase,data=mydata.GCH,family = "binomial")
summary(model)
# the function of step realizes stepwise regression method, trace=0 is for not print information during the running of step
model2<-step(object = model,trace = 0)
summary(model2)
# It was found from the results that the P values of all variables were less than 0.05, and the relatively important variables were retained through the significance test
anova(object = model2,test = "Chisq")
# Significance test of the model
# It can be seen from the above that as the variables are added to the model one by one from the first to the last, the model finally passes the significance test, indicating that the model composed of the above variables is meaningful and correct.
prob<-predict(object =model2,newdata=mydata.GCH,type = "response")
# Prediction accuracy of the model for test data (here, we select train data as test although there maybe is overfitting-related problems)
pred<-ifelse(prob>=0.5,"yes","no")
pred<-factor(pred,levels = c("no","yes"),order=TRUE)
mydata.GCH<-cbind(mydata.GCH,pred)
f<-table(mydata.GCH$Methyl,pred)
# draw ROC curve and calculate AUC
roc_curve <- roc(mydata.GCH$Methyl,prob)
x <- 1-roc_curve$specificities
y <- roc_curve$sensitivities
p <- ggplot(data = NULL, mapping = aes(x= x, y = y))
p + geom_line(colour = 'red') +geom_abline(intercept = 0, slope = 1)+ annotate('text', x = 0.4, y = 0.5, label =paste('AUC=',round(roc_curve$auc,2)))+ labs(x = '1-specificities',y = 'sensitivities', title = 'chr5: ROC Curve (Methyl~ChIP+Repli+DNase)')
write.table(mydata.GCH,"chr5.result.txt",quote=F,sep="\t",row.names=F)
# chr6
mydata.GCH<-read.table("chr6.GCH.bed",sep="\t")
# set colnames
colnames(mydata.GCH)<-c("chrom", "start", "end", "Methyl", "Repli", "ChIP", "DNase")
# what we care about in the model is the result (CTCT binding=yes), so we just sorted this factor.
mydata.GCH$Methyl<-factor(mydata.GCH$Methyl,levels = c("no","yes"),order=TRUE)
# construct Logistic model
model<-glm(Methyl~ChIP+Repli+DNase,data=mydata.GCH,family = "binomial")
summary(model)
# the function of step realizes stepwise regression method, trace=0 is for not print information during the running of step
model2<-step(object = model,trace = 0)
summary(model2)
# It was found from the results that the P values of all variables were less than 0.05, and the relatively important variables were retained through the significance test
anova(object = model2,test = "Chisq")
# Significance test of the model
# It can be seen from the above that as the variables are added to the model one by one from the first to the last, the model finally passes the significance test, indicating that the model composed of the above variables is meaningful and correct.
prob<-predict(object =model2,newdata=mydata.GCH,type = "response")
# Prediction accuracy of the model for test data (here, we select train data as test although there maybe is overfitting-related problems)
pred<-ifelse(prob>=0.5,"yes","no")
pred<-factor(pred,levels = c("no","yes"),order=TRUE)
mydata.GCH<-cbind(mydata.GCH,pred)
f<-table(mydata.GCH$Methyl,pred)
# draw ROC curve and calculate AUC
roc_curve <- roc(mydata.GCH$Methyl,prob)
x <- 1-roc_curve$specificities
y <- roc_curve$sensitivities
p <- ggplot(data = NULL, mapping = aes(x= x, y = y))
p + geom_line(colour = 'red') +geom_abline(intercept = 0, slope = 1)+ annotate('text', x = 0.4, y = 0.5, label =paste('AUC=',round(roc_curve$auc,2)))+ labs(x = '1-specificities',y = 'sensitivities', title = 'chr6: ROC Curve (Methyl~ChIP+Repli+DNase)')
write.table(mydata.GCH,"chr6.result.txt",quote=F,sep="\t",row.names=F)
# chr7
mydata.GCH<-read.table("chr7.GCH.bed",sep="\t")
# set colnames
colnames(mydata.GCH)<-c("chrom", "start", "end", "Methyl", "Repli", "ChIP", "DNase")
# what we care about in the model is the result (CTCT binding=yes), so we just sorted this factor.
mydata.GCH$Methyl<-factor(mydata.GCH$Methyl,levels = c("no","yes"),order=TRUE)
# construct Logistic model
model<-glm(Methyl~ChIP+Repli+DNase,data=mydata.GCH,family = "binomial")
summary(model)
# the function of step realizes stepwise regression method, trace=0 is for not print information during the running of step
model2<-step(object = model,trace = 0)
summary(model2)
# It was found from the results that the P values of all variables were less than 0.05, and the relatively important variables were retained through the significance test
anova(object = model2,test = "Chisq")
# Significance test of the model
# It can be seen from the above that as the variables are added to the model one by one from the first to the last, the model finally passes the significance test, indicating that the model composed of the above variables is meaningful and correct.
prob<-predict(object =model2,newdata=mydata.GCH,type = "response")
# Prediction accuracy of the model for test data (here, we select train data as test although there maybe is overfitting-related problems)
pred<-ifelse(prob>=0.5,"yes","no")
pred<-factor(pred,levels = c("no","yes"),order=TRUE)
mydata.GCH<-cbind(mydata.GCH,pred)
f<-table(mydata.GCH$Methyl,pred)
# draw ROC curve and calculate AUC
roc_curve <- roc(mydata.GCH$Methyl,prob)
x <- 1-roc_curve$specificities
y <- roc_curve$sensitivities
p <- ggplot(data = NULL, mapping = aes(x= x, y = y))
p + geom_line(colour = 'red') +geom_abline(intercept = 0, slope = 1)+ annotate('text', x = 0.4, y = 0.5, label =paste('AUC=',round(roc_curve$auc,2)))+ labs(x = '1-specificities',y = 'sensitivities', title = 'chr7: ROC Curve (Methyl~ChIP+Repli+DNase)')
write.table(mydata.GCH,"chr7.result.txt",quote=F,sep="\t",row.names=F)
# chr8
mydata.GCH<-read.table("chr8.GCH.bed",sep="\t")
# set colnames
colnames(mydata.GCH)<-c("chrom", "start", "end", "Methyl", "Repli", "ChIP", "DNase")
# what we care about in the model is the result (CTCT binding=yes), so we just sorted this factor.
mydata.GCH$Methyl<-factor(mydata.GCH$Methyl,levels = c("no","yes"),order=TRUE)
# construct Logistic model
model<-glm(Methyl~ChIP+Repli+DNase,data=mydata.GCH,family = "binomial")
summary(model)
# the function of step realizes stepwise regression method, trace=0 is for not print information during the running of step
model2<-step(object = model,trace = 0)
summary(model2)
# It was found from the results that the P values of all variables were less than 0.05, and the relatively important variables were retained through the significance test
anova(object = model2,test = "Chisq")
# Significance test of the model
# It can be seen from the above that as the variables are added to the model one by one from the first to the last, the model finally passes the significance test, indicating that the model composed of the above variables is meaningful and correct.
prob<-predict(object =model2,newdata=mydata.GCH,type = "response")
# Prediction accuracy of the model for test data (here, we select train data as test although there maybe is overfitting-related problems)
pred<-ifelse(prob>=0.5,"yes","no")
pred<-factor(pred,levels = c("no","yes"),order=TRUE)
mydata.GCH<-cbind(mydata.GCH,pred)
f<-table(mydata.GCH$Methyl,pred)
# draw ROC curve and calculate AUC
roc_curve <- roc(mydata.GCH$Methyl,prob)
x <- 1-roc_curve$specificities
y <- roc_curve$sensitivities
p <- ggplot(data = NULL, mapping = aes(x= x, y = y))
p + geom_line(colour = 'red') +geom_abline(intercept = 0, slope = 1)+ annotate('text', x = 0.4, y = 0.5, label =paste('AUC=',round(roc_curve$auc,2)))+ labs(x = '1-specificities',y = 'sensitivities', title = 'chr8: ROC Curve (Methyl~ChIP+Repli+DNase)')
write.table(mydata.GCH,"chr8.result.txt",quote=F,sep="\t",row.names=F)
# chr9
mydata.GCH<-read.table("chr9.GCH.bed",sep="\t")
# set colnames
colnames(mydata.GCH)<-c("chrom", "start", "end", "Methyl", "Repli", "ChIP", "DNase")
# what we care about in the model is the result (CTCT binding=yes), so we just sorted this factor.
mydata.GCH$Methyl<-factor(mydata.GCH$Methyl,levels = c("no","yes"),order=TRUE)
# construct Logistic model
model<-glm(Methyl~ChIP+Repli+DNase,data=mydata.GCH,family = "binomial")
summary(model)
# the function of step realizes stepwise regression method, trace=0 is for not print information during the running of step
model2<-step(object = model,trace = 0)
summary(model2)
# It was found from the results that the P values of all variables were less than 0.05, and the relatively important variables were retained through the significance test
anova(object = model2,test = "Chisq")
# Significance test of the model
# It can be seen from the above that as the variables are added to the model one by one from the first to the last, the model finally passes the significance test, indicating that the model composed of the above variables is meaningful and correct.
prob<-predict(object =model2,newdata=mydata.GCH,type = "response")
# Prediction accuracy of the model for test data (here, we select train data as test although there maybe is overfitting-related problems)
pred<-ifelse(prob>=0.5,"yes","no")
pred<-factor(pred,levels = c("no","yes"),order=TRUE)
mydata.GCH<-cbind(mydata.GCH,pred)
f<-table(mydata.GCH$Methyl,pred)
# draw ROC curve and calculate AUC
roc_curve <- roc(mydata.GCH$Methyl,prob)
x <- 1-roc_curve$specificities
y <- roc_curve$sensitivities
p <- ggplot(data = NULL, mapping = aes(x= x, y = y))
p + geom_line(colour = 'red') +geom_abline(intercept = 0, slope = 1)+ annotate('text', x = 0.4, y = 0.5, label =paste('AUC=',round(roc_curve$auc,2)))+ labs(x = '1-specificities',y = 'sensitivities', title = 'chr9: ROC Curve (Methyl~ChIP+Repli+DNase)')
write.table(mydata.GCH,"chr9.result.txt",quote=F,sep="\t",row.names=F)
# chr10
mydata.GCH<-read.table("chr10.GCH.bed",sep="\t")
# set colnames
colnames(mydata.GCH)<-c("chrom", "start", "end", "Methyl", "Repli", "ChIP", "DNase")
# what we care about in the model is the result (CTCT binding=yes), so we just sorted this factor.
mydata.GCH$Methyl<-factor(mydata.GCH$Methyl,levels = c("no","yes"),order=TRUE)
# construct Logistic model
model<-glm(Methyl~ChIP+Repli+DNase,data=mydata.GCH,family = "binomial")
summary(model)
# the function of step realizes stepwise regression method, trace=0 is for not print information during the running of step
model2<-step(object = model,trace = 0)
summary(model2)
# It was found from the results that the P values of all variables were less than 0.05, and the relatively important variables were retained through the significance test
anova(object = model2,test = "Chisq")
# Significance test of the model
# It can be seen from the above that as the variables are added to the model one by one from the first to the last, the model finally passes the significance test, indicating that the model composed of the above variables is meaningful and correct.
prob<-predict(object =model2,newdata=mydata.GCH,type = "response")
# Prediction accuracy of the model for test data (here, we select train data as test although there maybe is overfitting-related problems)
pred<-ifelse(prob>=0.5,"yes","no")
pred<-factor(pred,levels = c("no","yes"),order=TRUE)
mydata.GCH<-cbind(mydata.GCH,pred)
f<-table(mydata.GCH$Methyl,pred)
# draw ROC curve and calculate AUC
roc_curve <- roc(mydata.GCH$Methyl,prob)
x <- 1-roc_curve$specificities
y <- roc_curve$sensitivities
p <- ggplot(data = NULL, mapping = aes(x= x, y = y))
p + geom_line(colour = 'red') +geom_abline(intercept = 0, slope = 1)+ annotate('text', x = 0.4, y = 0.5, label =paste('AUC=',round(roc_curve$auc,2)))+ labs(x = '1-specificities',y = 'sensitivities', title = 'chr10: ROC Curve (Methyl~ChIP+Repli+DNase)')
write.table(mydata.GCH,"chr10.result.txt",quote=F,sep="\t",row.names=F)
# chr11
mydata.GCH<-read.table("chr11.GCH.bed",sep="\t")
# set colnames
colnames(mydata.GCH)<-c("chrom", "start", "end", "Methyl", "Repli", "ChIP", "DNase")
# what we care about in the model is the result (CTCT binding=yes), so we just sorted this factor.
mydata.GCH$Methyl<-factor(mydata.GCH$Methyl,levels = c("no","yes"),order=TRUE)
# construct Logistic model
model<-glm(Methyl~ChIP+Repli+DNase,data=mydata.GCH,family = "binomial")
summary(model)
# the function of step realizes stepwise regression method, trace=0 is for not print information during the running of step
model2<-step(object = model,trace = 0)
summary(model2)
# It was found from the results that the P values of all variables were less than 0.05, and the relatively important variables were retained through the significance test
anova(object = model2,test = "Chisq")
# Significance test of the model
# It can be seen from the above that as the variables are added to the model one by one from the first to the last, the model finally passes the significance test, indicating that the model composed of the above variables is meaningful and correct.
prob<-predict(object =model2,newdata=mydata.GCH,type = "response")
# Prediction accuracy of the model for test data (here, we select train data as test although there maybe is overfitting-related problems)
pred<-ifelse(prob>=0.5,"yes","no")
pred<-factor(pred,levels = c("no","yes"),order=TRUE)
mydata.GCH<-cbind(mydata.GCH,pred)
f<-table(mydata.GCH$Methyl,pred)
# draw ROC curve and calculate AUC
roc_curve <- roc(mydata.GCH$Methyl,prob)
x <- 1-roc_curve$specificities
y <- roc_curve$sensitivities
p <- ggplot(data = NULL, mapping = aes(x= x, y = y))
p + geom_line(colour = 'red') +geom_abline(intercept = 0, slope = 1)+ annotate('text', x = 0.4, y = 0.5, label =paste('AUC=',round(roc_curve$auc,2)))+ labs(x = '1-specificities',y = 'sensitivities', title = 'chr11: ROC Curve (Methyl~ChIP+Repli+DNase)')
write.table(mydata.GCH,"chr11.result.txt",quote=F,sep="\t",row.names=F)
# chr12
mydata.GCH<-read.table("chr12.GCH.bed",sep="\t")
# set colnames
colnames(mydata.GCH)<-c("chrom", "start", "end", "Methyl", "Repli", "ChIP", "DNase")
# what we care about in the model is the result (CTCT binding=yes), so we just sorted this factor.
mydata.GCH$Methyl<-factor(mydata.GCH$Methyl,levels = c("no","yes"),order=TRUE)
# construct Logistic model
model<-glm(Methyl~ChIP+Repli+DNase,data=mydata.GCH,family = "binomial")
summary(model)
# the function of step realizes stepwise regression method, trace=0 is for not print information during the running of step
model2<-step(object = model,trace = 0)
summary(model2)
# It was found from the results that the P values of all variables were less than 0.05, and the relatively important variables were retained through the significance test
anova(object = model2,test = "Chisq")
# Significance test of the model
# It can be seen from the above that as the variables are added to the model one by one from the first to the last, the model finally passes the significance test, indicating that the model composed of the above variables is meaningful and correct.
prob<-predict(object =model2,newdata=mydata.GCH,type = "response")
# Prediction accuracy of the model for test data (here, we select train data as test although there maybe is overfitting-related problems)
pred<-ifelse(prob>=0.5,"yes","no")
pred<-factor(pred,levels = c("no","yes"),order=TRUE)
mydata.GCH<-cbind(mydata.GCH,pred)
f<-table(mydata.GCH$Methyl,pred)
# draw ROC curve and calculate AUC
roc_curve <- roc(mydata.GCH$Methyl,prob)
x <- 1-roc_curve$specificities
y <- roc_curve$sensitivities
p <- ggplot(data = NULL, mapping = aes(x= x, y = y))
p + geom_line(colour = 'red') +geom_abline(intercept = 0, slope = 1)+ annotate('text', x = 0.4, y = 0.5, label =paste('AUC=',round(roc_curve$auc,2)))+ labs(x = '1-specificities',y = 'sensitivities', title = 'chr12: ROC Curve (Methyl~ChIP+Repli+DNase)')
write.table(mydata.GCH,"chr12.result.txt",quote=F,sep="\t",row.names=F)
# chr13
mydata.GCH<-read.table("chr13.GCH.bed",sep="\t")
# set colnames
colnames(mydata.GCH)<-c("chrom", "start", "end", "Methyl", "Repli", "ChIP", "DNase")
# what we care about in the model is the result (CTCT binding=yes), so we just sorted this factor.
mydata.GCH$Methyl<-factor(mydata.GCH$Methyl,levels = c("no","yes"),order=TRUE)
# construct Logistic model
model<-glm(Methyl~ChIP+Repli+DNase,data=mydata.GCH,family = "binomial")
summary(model)
# the function of step realizes stepwise regression method, trace=0 is for not print information during the running of step
model2<-step(object = model,trace = 0)
summary(model2)
# It was found from the results that the P values of all variables were less than 0.05, and the relatively important variables were retained through the significance test
anova(object = model2,test = "Chisq")
# Significance test of the model
# It can be seen from the above that as the variables are added to the model one by one from the first to the last, the model finally passes the significance test, indicating that the model composed of the above variables is meaningful and correct.
prob<-predict(object =model2,newdata=mydata.GCH,type = "response")
# Prediction accuracy of the model for test data (here, we select train data as test although there maybe is overfitting-related problems)
pred<-ifelse(prob>=0.5,"yes","no")
pred<-factor(pred,levels = c("no","yes"),order=TRUE)
mydata.GCH<-cbind(mydata.GCH,pred)
f<-table(mydata.GCH$Methyl,pred)
# draw ROC curve and calculate AUC
roc_curve <- roc(mydata.GCH$Methyl,prob)
x <- 1-roc_curve$specificities
y <- roc_curve$sensitivities
p <- ggplot(data = NULL, mapping = aes(x= x, y = y))
p + geom_line(colour = 'red') +geom_abline(intercept = 0, slope = 1)+ annotate('text', x = 0.4, y = 0.5, label =paste('AUC=',round(roc_curve$auc,2)))+ labs(x = '1-specificities',y = 'sensitivities', title = 'chr13: ROC Curve (Methyl~ChIP+Repli+DNase)')
write.table(mydata.GCH,"chr13.result.txt",quote=F,sep="\t",row.names=F)
# chr14
mydata.GCH<-read.table("chr14.GCH.bed",sep="\t")
# set colnames
colnames(mydata.GCH)<-c("chrom", "start", "end", "Methyl", "Repli", "ChIP", "DNase")
# what we care about in the model is the result (CTCT binding=yes), so we just sorted this factor.
mydata.GCH$Methyl<-factor(mydata.GCH$Methyl,levels = c("no","yes"),order=TRUE)
# construct Logistic model
model<-glm(Methyl~ChIP+Repli+DNase,data=mydata.GCH,family = "binomial")
summary(model)
# the function of step realizes stepwise regression method, trace=0 is for not print information during the running of step
model2<-step(object = model,trace = 0)
summary(model2)
# It was found from the results that the P values of all variables were less than 0.05, and the relatively important variables were retained through the significance test
anova(object = model2,test = "Chisq")
# Significance test of the model
# It can be seen from the above that as the variables are added to the model one by one from the first to the last, the model finally passes the significance test, indicating that the model composed of the above variables is meaningful and correct.
prob<-predict(object =model2,newdata=mydata.GCH,type = "response")
# Prediction accuracy of the model for test data (here, we select train data as test although there maybe is overfitting-related problems)
pred<-ifelse(prob>=0.5,"yes","no")
pred<-factor(pred,levels = c("no","yes"),order=TRUE)
mydata.GCH<-cbind(mydata.GCH,pred)
f<-table(mydata.GCH$Methyl,pred)
# draw ROC curve and calculate AUC
roc_curve <- roc(mydata.GCH$Methyl,prob)
x <- 1-roc_curve$specificities
y <- roc_curve$sensitivities
p <- ggplot(data = NULL, mapping = aes(x= x, y = y))
p + geom_line(colour = 'red') +geom_abline(intercept = 0, slope = 1)+ annotate('text', x = 0.4, y = 0.5, label =paste('AUC=',round(roc_curve$auc,2)))+ labs(x = '1-specificities',y = 'sensitivities', title = 'chr14: ROC Curve (Methyl~ChIP+Repli+DNase)')
write.table(mydata.GCH,"chr14.result.txt",quote=F,sep="\t",row.names=F)
# chr15
mydata.GCH<-read.table("chr15.GCH.bed",sep="\t")
# set colnames
colnames(mydata.GCH)<-c("chrom", "start", "end", "Methyl", "Repli", "ChIP", "DNase")
# what we care about in the model is the result (CTCT binding=yes), so we just sorted this factor.
mydata.GCH$Methyl<-factor(mydata.GCH$Methyl,levels = c("no","yes"),order=TRUE)
# construct Logistic model
model<-glm(Methyl~ChIP+Repli+DNase,data=mydata.GCH,family = "binomial")
summary(model)
# the function of step realizes stepwise regression method, trace=0 is for not print information during the running of step
model2<-step(object = model,trace = 0)
summary(model2)
# It was found from the results that the P values of all variables were less than 0.05, and the relatively important variables were retained through the significance test
anova(object = model2,test = "Chisq")
# Significance test of the model
# It can be seen from the above that as the variables are added to the model one by one from the first to the last, the model finally passes the significance test, indicating that the model composed of the above variables is meaningful and correct.
prob<-predict(object =model2,newdata=mydata.GCH,type = "response")
# Prediction accuracy of the model for test data (here, we select train data as test although there maybe is overfitting-related problems)
pred<-ifelse(prob>=0.5,"yes","no")
pred<-factor(pred,levels = c("no","yes"),order=TRUE)
mydata.GCH<-cbind(mydata.GCH,pred)
f<-table(mydata.GCH$Methyl,pred)
# draw ROC curve and calculate AUC
roc_curve <- roc(mydata.GCH$Methyl,prob)
x <- 1-roc_curve$specificities
y <- roc_curve$sensitivities
p <- ggplot(data = NULL, mapping = aes(x= x, y = y))
p + geom_line(colour = 'red') +geom_abline(intercept = 0, slope = 1)+ annotate('text', x = 0.4, y = 0.5, label =paste('AUC=',round(roc_curve$auc,2)))+ labs(x = '1-specificities',y = 'sensitivities', title = 'chr15: ROC Curve (Methyl~ChIP+Repli+DNase)')
write.table(mydata.GCH,"chr15.result.txt",quote=F,sep="\t",row.names=F)
# chr16
mydata.GCH<-read.table("chr16.GCH.bed",sep="\t")
# set colnames
colnames(mydata.GCH)<-c("chrom", "start", "end", "Methyl", "Repli", "ChIP", "DNase")
# what we care about in the model is the result (CTCT binding=yes), so we just sorted this factor.
mydata.GCH$Methyl<-factor(mydata.GCH$Methyl,levels = c("no","yes"),order=TRUE)
# construct Logistic model
model<-glm(Methyl~ChIP+Repli+DNase,data=mydata.GCH,family = "binomial")
summary(model)
# the function of step realizes stepwise regression method, trace=0 is for not print information during the running of step
model2<-step(object = model,trace = 0)
summary(model2)
# It was found from the results that the P values of all variables were less than 0.05, and the relatively important variables were retained through the significance test
anova(object = model2,test = "Chisq")
# Significance test of the model
# It can be seen from the above that as the variables are added to the model one by one from the first to the last, the model finally passes the significance test, indicating that the model composed of the above variables is meaningful and correct.
prob<-predict(object =model2,newdata=mydata.GCH,type = "response")
# Prediction accuracy of the model for test data (here, we select train data as test although there maybe is overfitting-related problems)
pred<-ifelse(prob>=0.5,"yes","no")
pred<-factor(pred,levels = c("no","yes"),order=TRUE)
mydata.GCH<-cbind(mydata.GCH,pred)
f<-table(mydata.GCH$Methyl,pred)
# draw ROC curve and calculate AUC
roc_curve <- roc(mydata.GCH$Methyl,prob)
x <- 1-roc_curve$specificities
y <- roc_curve$sensitivities
p <- ggplot(data = NULL, mapping = aes(x= x, y = y))
p + geom_line(colour = 'red') +geom_abline(intercept = 0, slope = 1)+ annotate('text', x = 0.4, y = 0.5, label =paste('AUC=',round(roc_curve$auc,2)))+ labs(x = '1-specificities',y = 'sensitivities', title = 'chr16: ROC Curve (Methyl~ChIP+Repli+DNase)')
write.table(mydata.GCH,"chr16.result.txt",quote=F,sep="\t",row.names=F)
# chr17
mydata.GCH<-read.table("chr17.GCH.bed",sep="\t")
# set colnames
colnames(mydata.GCH)<-c("chrom", "start", "end", "Methyl", "Repli", "ChIP", "DNase")
# what we care about in the model is the result (CTCT binding=yes), so we just sorted this factor.
mydata.GCH$Methyl<-factor(mydata.GCH$Methyl,levels = c("no","yes"),order=TRUE)
# construct Logistic model
model<-glm(Methyl~ChIP+Repli+DNase,data=mydata.GCH,family = "binomial")
summary(model)
# the function of step realizes stepwise regression method, trace=0 is for not print information during the running of step
model2<-step(object = model,trace = 0)
summary(model2)
# It was found from the results that the P values of all variables were less than 0.05, and the relatively important variables were retained through the significance test
anova(object = model2,test = "Chisq")
# Significance test of the model
# It can be seen from the above that as the variables are added to the model one by one from the first to the last, the model finally passes the significance test, indicating that the model composed of the above variables is meaningful and correct.
prob<-predict(object =model2,newdata=mydata.GCH,type = "response")
# Prediction accuracy of the model for test data (here, we select train data as test although there maybe is overfitting-related problems)
pred<-ifelse(prob>=0.5,"yes","no")
pred<-factor(pred,levels = c("no","yes"),order=TRUE)
mydata.GCH<-cbind(mydata.GCH,pred)
f<-table(mydata.GCH$Methyl,pred)
# draw ROC curve and calculate AUC
roc_curve <- roc(mydata.GCH$Methyl,prob)
x <- 1-roc_curve$specificities
y <- roc_curve$sensitivities
p <- ggplot(data = NULL, mapping = aes(x= x, y = y))
p + geom_line(colour = 'red') +geom_abline(intercept = 0, slope = 1)+ annotate('text', x = 0.4, y = 0.5, label =paste('AUC=',round(roc_curve$auc,2)))+ labs(x = '1-specificities',y = 'sensitivities', title = 'chr17: ROC Curve (Methyl~ChIP+Repli+DNase)')
write.table(mydata.GCH,"chr17.result.txt",quote=F,sep="\t",row.names=F)
# chr18
mydata.GCH<-read.table("chr18.GCH.bed",sep="\t")
# set colnames
colnames(mydata.GCH)<-c("chrom", "start", "end", "Methyl", "Repli", "ChIP", "DNase")
# what we care about in the model is the result (CTCT binding=yes), so we just sorted this factor.
mydata.GCH$Methyl<-factor(mydata.GCH$Methyl,levels = c("no","yes"),order=TRUE)
# construct Logistic model
model<-glm(Methyl~ChIP+Repli+DNase,data=mydata.GCH,family = "binomial")
summary(model)
# the function of step realizes stepwise regression method, trace=0 is for not print information during the running of step
model2<-step(object = model,trace = 0)
summary(model2)
# It was found from the results that the P values of all variables were less than 0.05, and the relatively important variables were retained through the significance test
anova(object = model2,test = "Chisq")
# Significance test of the model
# It can be seen from the above that as the variables are added to the model one by one from the first to the last, the model finally passes the significance test, indicating that the model composed of the above variables is meaningful and correct.
prob<-predict(object =model2,newdata=mydata.GCH,type = "response")
# Prediction accuracy of the model for test data (here, we select train data as test although there maybe is overfitting-related problems)
pred<-ifelse(prob>=0.5,"yes","no")
pred<-factor(pred,levels = c("no","yes"),order=TRUE)
mydata.GCH<-cbind(mydata.GCH,pred)
f<-table(mydata.GCH$Methyl,pred)
# draw ROC curve and calculate AUC
roc_curve <- roc(mydata.GCH$Methyl,prob)
x <- 1-roc_curve$specificities
y <- roc_curve$sensitivities
p <- ggplot(data = NULL, mapping = aes(x= x, y = y))
p + geom_line(colour = 'red') +geom_abline(intercept = 0, slope = 1)+ annotate('text', x = 0.4, y = 0.5, label =paste('AUC=',round(roc_curve$auc,2)))+ labs(x = '1-specificities',y = 'sensitivities', title = 'chr18: ROC Curve (Methyl~ChIP+Repli+DNase)')
write.table(mydata.GCH,"chr18.result.txt",quote=F,sep="\t",row.names=F)
# chr19
mydata.GCH<-read.table("chr19.GCH.bed",sep="\t")
# set colnames
colnames(mydata.GCH)<-c("chrom", "start", "end", "Methyl", "Repli", "ChIP", "DNase")
# what we care about in the model is the result (CTCT binding=yes), so we just sorted this factor.
mydata.GCH$Methyl<-factor(mydata.GCH$Methyl,levels = c("no","yes"),order=TRUE)
# construct Logistic model
model<-glm(Methyl~ChIP+Repli+DNase,data=mydata.GCH,family = "binomial")
summary(model)
# the function of step realizes stepwise regression method, trace=0 is for not print information during the running of step
model2<-step(object = model,trace = 0)
summary(model2)
# It was found from the results that the P values of all variables were less than 0.05, and the relatively important variables were retained through the significance test
anova(object = model2,test = "Chisq")
# Significance test of the model
# It can be seen from the above that as the variables are added to the model one by one from the first to the last, the model finally passes the significance test, indicating that the model composed of the above variables is meaningful and correct.
prob<-predict(object =model2,newdata=mydata.GCH,type = "response")
# Prediction accuracy of the model for test data (here, we select train data as test although there maybe is overfitting-related problems)
pred<-ifelse(prob>=0.5,"yes","no")
pred<-factor(pred,levels = c("no","yes"),order=TRUE)
mydata.GCH<-cbind(mydata.GCH,pred)
f<-table(mydata.GCH$Methyl,pred)
# draw ROC curve and calculate AUC
roc_curve <- roc(mydata.GCH$Methyl,prob)
x <- 1-roc_curve$specificities
y <- roc_curve$sensitivities
p <- ggplot(data = NULL, mapping = aes(x= x, y = y))
p + geom_line(colour = 'red') +geom_abline(intercept = 0, slope = 1)+ annotate('text', x = 0.4, y = 0.5, label =paste('AUC=',round(roc_curve$auc,2)))+ labs(x = '1-specificities',y = 'sensitivities', title = 'chr19: ROC Curve (Methyl~ChIP+Repli+DNase)')
write.table(mydata.GCH,"chr19.result.txt",quote=F,sep="\t",row.names=F)
# chr20
mydata.GCH<-read.table("chr20.GCH.bed",sep="\t")
# set colnames
colnames(mydata.GCH)<-c("chrom", "start", "end", "Methyl", "Repli", "ChIP", "DNase")
# what we care about in the model is the result (CTCT binding=yes), so we just sorted this factor.
mydata.GCH$Methyl<-factor(mydata.GCH$Methyl,levels = c("no","yes"),order=TRUE)
# construct Logistic model
model<-glm(Methyl~ChIP+Repli+DNase,data=mydata.GCH,family = "binomial")
summary(model)
# the function of step realizes stepwise regression method, trace=0 is for not print information during the running of step
model2<-step(object = model,trace = 0)
summary(model2)
# It was found from the results that the P values of all variables were less than 0.05, and the relatively important variables were retained through the significance test
anova(object = model2,test = "Chisq")
# Significance test of the model
# It can be seen from the above that as the variables are added to the model one by one from the first to the last, the model finally passes the significance test, indicating that the model composed of the above variables is meaningful and correct.
prob<-predict(object =model2,newdata=mydata.GCH,type = "response")
# Prediction accuracy of the model for test data (here, we select train data as test although there maybe is overfitting-related problems)
pred<-ifelse(prob>=0.5,"yes","no")
pred<-factor(pred,levels = c("no","yes"),order=TRUE)
mydata.GCH<-cbind(mydata.GCH,pred)
f<-table(mydata.GCH$Methyl,pred)
# draw ROC curve and calculate AUC
roc_curve <- roc(mydata.GCH$Methyl,prob)
x <- 1-roc_curve$specificities
y <- roc_curve$sensitivities
p <- ggplot(data = NULL, mapping = aes(x= x, y = y))
p + geom_line(colour = 'red') +geom_abline(intercept = 0, slope = 1)+ annotate('text', x = 0.4, y = 0.5, label =paste('AUC=',round(roc_curve$auc,2)))+ labs(x = '1-specificities',y = 'sensitivities', title = 'chr20: ROC Curve (Methyl~ChIP+Repli+DNase)')
write.table(mydata.GCH,"chr20.result.txt",quote=F,sep="\t",row.names=F)
# chr21
mydata.GCH<-read.table("chr21.GCH.bed",sep="\t")
# set colnames
colnames(mydata.GCH)<-c("chrom", "start", "end", "Methyl", "Repli", "ChIP", "DNase")
# what we care about in the model is the result (CTCT binding=yes), so we just sorted this factor.
mydata.GCH$Methyl<-factor(mydata.GCH$Methyl,levels = c("no","yes"),order=TRUE)
# construct Logistic model
model<-glm(Methyl~ChIP+Repli+DNase,data=mydata.GCH,family = "binomial")
summary(model)
# the function of step realizes stepwise regression method, trace=0 is for not print information during the running of step
model2<-step(object = model,trace = 0)
summary(model2)
# It was found from the results that the P values of all variables were less than 0.05, and the relatively important variables were retained through the significance test
anova(object = model2,test = "Chisq")
# Significance test of the model
# It can be seen from the above that as the variables are added to the model one by one from the first to the last, the model finally passes the significance test, indicating that the model composed of the above variables is meaningful and correct.
prob<-predict(object =model2,newdata=mydata.GCH,type = "response")
# Prediction accuracy of the model for test data (here, we select train data as test although there maybe is overfitting-related problems)
pred<-ifelse(prob>=0.5,"yes","no")
pred<-factor(pred,levels = c("no","yes"),order=TRUE)
mydata.GCH<-cbind(mydata.GCH,pred)
f<-table(mydata.GCH$Methyl,pred)
# draw ROC curve and calculate AUC
roc_curve <- roc(mydata.GCH$Methyl,prob)
x <- 1-roc_curve$specificities
y <- roc_curve$sensitivities
p <- ggplot(data = NULL, mapping = aes(x= x, y = y))
p + geom_line(colour = 'red') +geom_abline(intercept = 0, slope = 1)+ annotate('text', x = 0.4, y = 0.5, label =paste('AUC=',round(roc_curve$auc,2)))+ labs(x = '1-specificities',y = 'sensitivities', title = 'chr21: ROC Curve (Methyl~ChIP+Repli+DNase)')
write.table(mydata.GCH,"chr21.result.txt",quote=F,sep="\t",row.names=F)
# chrX
mydata.GCH<-read.table("chrX.GCH.bed",sep="\t")
# set colnames
colnames(mydata.GCH)<-c("chrom", "start", "end", "Methyl", "Repli", "ChIP", "DNase")
# what we care about in the model is the result (CTCT binding=yes), so we just sorted this factor.
mydata.GCH$Methyl<-factor(mydata.GCH$Methyl,levels = c("no","yes"),order=TRUE)
# construct Logistic model
model<-glm(Methyl~ChIP+Repli+DNase,data=mydata.GCH,family = "binomial")
summary(model)
# the function of step realizes stepwise regression method, trace=0 is for not print information during the running of step
model2<-step(object = model,trace = 0)
summary(model2)
# It was found from the results that the P values of all variables were less than 0.05, and the relatively important variables were retained through the significance test
anova(object = model2,test = "Chisq")
# Significance test of the model
# It can be seen from the above that as the variables are added to the model one by one from the first to the last, the model finally passes the significance test, indicating that the model composed of the above variables is meaningful and correct.
prob<-predict(object =model2,newdata=mydata.GCH,type = "response")
# Prediction accuracy of the model for test data (here, we select train data as test although there maybe is overfitting-related problems)
pred<-ifelse(prob>=0.5,"yes","no")
pred<-factor(pred,levels = c("no","yes"),order=TRUE)
mydata.GCH<-cbind(mydata.GCH,pred)
f<-table(mydata.GCH$Methyl,pred)
# draw ROC curve and calculate AUC
roc_curve <- roc(mydata.GCH$Methyl,prob)
x <- 1-roc_curve$specificities
y <- roc_curve$sensitivities
p <- ggplot(data = NULL, mapping = aes(x= x, y = y))
p + geom_line(colour = 'red') +geom_abline(intercept = 0, slope = 1)+ annotate('text', x = 0.4, y = 0.5, label =paste('AUC=',round(roc_curve$auc,2)))+ labs(x = '1-specificities',y = 'sensitivities', title = 'chrX: ROC Curve (Methyl~ChIP+Repli+DNase)')
write.table(mydata.GCH,"chrX.result.txt",quote=F,sep="\t",row.names=F)
# chrY
mydata.GCH<-read.table("chrY.GCH.bed",sep="\t")
# set colnames
colnames(mydata.GCH)<-c("chrom", "start", "end", "Methyl", "Repli", "ChIP", "DNase")
# what we care about in the model is the result (CTCT binding=yes), so we just sorted this factor.
mydata.GCH$Methyl<-factor(mydata.GCH$Methyl,levels = c("no","yes"),order=TRUE)
# construct Logistic model
model<-glm(Methyl~ChIP+Repli+DNase,data=mydata.GCH,family = "binomial")
summary(model)
# the function of step realizes stepwise regression method, trace=0 is for not print information during the running of step
model2<-step(object = model,trace = 0)
summary(model2)
# It was found from the results that the P values of all variables were less than 0.05, and the relatively important variables were retained through the significance test
anova(object = model2,test = "Chisq")
# Significance test of the model
# It can be seen from the above that as the variables are added to the model one by one from the first to the last, the model finally passes the significance test, indicating that the model composed of the above variables is meaningful and correct.
prob<-predict(object =model2,newdata=mydata.GCH,type = "response")
# Prediction accuracy of the model for test data (here, we select train data as test although there maybe is overfitting-related problems)
pred<-ifelse(prob>=0.5,"yes","no")
pred<-factor(pred,levels = c("no","yes"),order=TRUE)
mydata.GCH<-cbind(mydata.GCH,pred)
f<-table(mydata.GCH$Methyl,pred)
# draw ROC curve and calculate AUC
roc_curve <- roc(mydata.GCH$Methyl,prob)
x <- 1-roc_curve$specificities
y <- roc_curve$sensitivities
p <- ggplot(data = NULL, mapping = aes(x= x, y = y))
p + geom_line(colour = 'red') +geom_abline(intercept = 0, slope = 1)+ annotate('text', x = 0.4, y = 0.5, label =paste('AUC=',round(roc_curve$auc,2)))+ labs(x = '1-specificities',y = 'sensitivities', title = 'chrY: ROC Curve (Methyl~ChIP+Repli+DNase)')
write.table(mydata.GCH,"chrY.result.txt",quote=F,sep="\t",row.names=F)
dev.off()
