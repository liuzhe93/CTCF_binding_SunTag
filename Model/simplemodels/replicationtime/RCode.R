setwd("E:/MyProject_YanLab/Project4_CTCFbinding/Model/simplemodels/replicationtime/windowsize1000.max")
#setwd("E:/MyProject_YanLab/Project4_CTCFbinding/Model/simplemodels/replicationtime/windowsize1000.avg")
#setwd("E:/MyProject_YanLab/Project4_CTCFbinding/Model/simplemodels/replicationtime/windowsize500.max")
#setwd("E:/MyProject_YanLab/Project4_CTCFbinding/Model/simplemodels/replicationtime/windowsize500.avg")
#setwd("E:/MyProject_YanLab/Project4_CTCFbinding/Model/simplemodels/replicationtime/windowsize300.max")
#setwd("E:/MyProject_YanLab/Project4_CTCFbinding/Model/simplemodels/replicationtime/windowsize300.avg")

library(tidyverse)
library(tidyr)
library("ggplot2")
library("easyGgplot2")

pdf("merge1-x.pdf",width=9.14, height=9.11)
# chr1+2
# Section A: For GCH
# read data from file
mydata.GCH_1<-read.table("chr1.GCH.bed",sep="\t")
# set column names
mydata.GCH_1<-cbind(mydata.GCH_1[,1:3],log(mydata.GCH_1[,4]),mydata.GCH_1[,5])
colnames(mydata.GCH_1)<-c("chrom", "start", "end", "Methyl_transform", "Repli")
# remove x-axis cutoff < -10
mydata.GCH_1<-mydata.GCH_1[with(mydata.GCH_1,mydata.GCH_1$Methyl_transform>-10),]
# Estimation of regression parameters
# fit linear models and carry out regression
model_1<-lm(Methyl_transform~Repli,data=mydata.GCH_1)
# Obtain predicted and residual values and draw figures
mydata.GCH_1$predicted <- predict(model_1)
mydata.GCH_1$residuals <- residuals(model_1)
suma<-summary(model_1)
# R2 represents the proportion of variance, in the outcome variable y, that may be predicted by knowing the value of the x variables. 
# An R2 value close to 1 indicates that the model explains a large portion of the variance in the outcome variable.
result_r.squared<-round(suma$r.squared,3)
# the coefficient of RT
a_2<-round(model_1$coefficients[2],3)
a_3<-round(model_1$coefficients[1],3)
n=dim(mydata.GCH_1)[1]
# combine
if(a_3<0){result1<-paste0("Methyl_transform=",a_2,"XRepli",a_3,"; R^2=",result_r.squared,"; n=",n)}else{result1<-paste0("Methyl_transform=",a_2,"XRepli+",a_3,"; R^2=",result_r.squared,"; n=",n)}

p1=ggplot()+
	geom_hex(mapping=aes(x=mydata.GCH_1$Methyl,y=mydata.GCH_1$predicted),color="black")+
	labs(x = "Observed", y = "Model-predicted", title = "chr1: The Scatter Plot of mGpC(JY608-bg)",subtitle=result1)+
	theme(plot.title = element_text(hjust = 0.5,vjust = 1),plot.subtitle = element_text(hjust = 0.5,size=6))+
	expand_limits(x = 0)
	
# Section B: For WCG
# read data from file
mydata.WCG_2<-read.table("chr1.WCG.bed",sep="\t")
# set column names
mydata.WCG_2<-cbind(mydata.WCG_2[,1:3],log(mydata.WCG_2[,4]),mydata.WCG_2[,5])
colnames(mydata.WCG_2)<-c("chrom", "start", "end", "Methyl_transform", "Repli")
# remove x-axis cutoff < -10
mydata.WCG_2<-mydata.WCG_2[with(mydata.WCG_2,mydata.WCG_2$Methyl_transform>-10),]
# Estimation of regression parameters
# fit linear models and carry out regression
model_2<-lm(Methyl_transform~Repli,data=mydata.WCG_2)
# Obtain predicted and residual values and draw figures
mydata.WCG_2$predicted <- predict(model_2)
mydata.WCG_2$residuals <- residuals(model_2)
suma<-summary(model_2)
# R2 represents the proportion of variance, in the outcome variable y, that may be predicted by knowing the value of the x variables. 
# An R2 value close to 1 indicates that the model explains a large portion of the variance in the outcome variable.
result_r.squared<-round(suma$r.squared,3)
# the coefficient of RT
a_2<-round(model_2$coefficients[2],3)
a_3<-round(model_2$coefficients[1],3)
n=dim(mydata.WCG_2)[1]
# combine
if(a_3<0){result2<-paste0("Methyl_transform=",a_2,"XRepli",a_3,"; R^2=",result_r.squared,"; n=",n)}else{result2<-paste0("Methyl_transform=",a_2,"XRepli+",a_3,"; R^2=",result_r.squared,"; n=",n)}

p2=ggplot()+
	geom_hex(mapping=aes(x=mydata.WCG_2$Methyl,y=mydata.WCG_2$predicted),color="black")+
	labs(x = "Observed", y = "Model-predicted", title = "chr1: The Scatter Plot of mCpG(JY608-bg)",subtitle=result2)+
	theme(plot.title = element_text(hjust = 0.5,vjust = 1),plot.subtitle = element_text(hjust = 0.5,size=6))+
	expand_limits(x = 0)
	
# Section A: For GCH
# read data from file
mydata.GCH_3<-read.table("chr2.GCH.bed",sep="\t")
# set column names
mydata.GCH_3<-cbind(mydata.GCH_3[,1:3],log(mydata.GCH_3[,4]),mydata.GCH_3[,5])
colnames(mydata.GCH_3)<-c("chrom", "start", "end", "Methyl_transform", "Repli")
# remove x-axis cutoff < -10
mydata.GCH_3<-mydata.GCH_3[with(mydata.GCH_3,mydata.GCH_3$Methyl_transform>-10),]
# Estimation of regression parameters
# fit linear models and carry out regression
model_3<-lm(Methyl_transform~Repli,data=mydata.GCH_3)
# Obtain predicted and residual values and draw figures
mydata.GCH_3$predicted <- predict(model_3)
mydata.GCH_3$residuals <- residuals(model_3)
suma<-summary(model_3)
# R2 represents the proportion of variance, in the outcome variable y, that may be predicted by knowing the value of the x variables. 
# An R2 value close to 1 indicates that the model explains a large portion of the variance in the outcome variable.
result_r.squared<-round(suma$r.squared,3)
# the coefficient of RT
a_2<-round(model_3$coefficients[2],3)
a_3<-round(model_3$coefficients[1],3)
n=dim(mydata.GCH_3)[1]
# combine
if(a_3<0){result3<-paste0("Methyl_transform=",a_2,"XRepli",a_3,"; R^2=",result_r.squared,"; n=",n)}else{result3<-paste0("Methyl_transform=",a_2,"XRepli+",a_3,"; R^2=",result_r.squared,"; n=",n)}

p3=ggplot()+
	geom_hex(mapping=aes(x=mydata.GCH_3$Methyl,y=mydata.GCH_3$predicted),color="black")+
	labs(x = "Observed", y = "Model-predicted", title = "chr2: The Scatter Plot of mGpC(JY608-bg)",subtitle=result3)+
	theme(plot.title = element_text(hjust = 0.5,vjust = 1),plot.subtitle = element_text(hjust = 0.5,size=6))+
	expand_limits(x = 0)
	
# Section B: For WCG
# read data from file
mydata.WCG_4<-read.table("chr2.WCG.bed",sep="\t")
# set column names
mydata.WCG_4<-cbind(mydata.WCG_4[,1:3],log(mydata.WCG_4[,4]),mydata.WCG_4[,5])
colnames(mydata.WCG_4)<-c("chrom", "start", "end", "Methyl_transform", "Repli")
# remove x-axis cutoff < -10
mydata.WCG_4<-mydata.WCG_4[with(mydata.WCG_4,mydata.WCG_4$Methyl_transform>-10),]
# Estimation of regression parameters
# fit linear models and carry out regression
model_4<-lm(Methyl_transform~Repli,data=mydata.WCG_4)
# Obtain predicted and residual values and draw figures
mydata.WCG_4$predicted <- predict(model_4)
mydata.WCG_4$residuals <- residuals(model_4)
suma<-summary(model_4)
# R2 represents the proportion of variance, in the outcome variable y, that may be predicted by knowing the value of the x variables. 
# An R2 value close to 1 indicates that the model explains a large portion of the variance in the outcome variable.
result_r.squared<-round(suma$r.squared,3)
# the coefficient of RT
a_2<-round(model_4$coefficients[2],3)
a_3<-round(model_4$coefficients[1],3)
n=dim(mydata.WCG_4)[1]
# combine
if(a_3<0){result4<-paste0("Methyl_transform=",a_2,"XRepli",a_3,"; R^2=",result_r.squared,"; n=",n)}else{result4<-paste0("Methyl_transform=",a_2,"XRepli+",a_3,"; R^2=",result_r.squared,"; n=",n)}

p4=ggplot()+
	geom_hex(mapping=aes(x=mydata.WCG_4$Methyl,y=mydata.WCG_4$predicted),color="black")+
	labs(x = "Observed", y = "Model-predicted", title = "chr2: The Scatter Plot of mCpG(JY608-bg)",subtitle=result4)+
	theme(plot.title = element_text(hjust = 0.5,vjust = 1),plot.subtitle = element_text(hjust = 0.5,size=6))+
	expand_limits(x = 0)
	
ggplot2.multiplot(p1,p2,p3,p4,cols=2)		



# chr3+4
# Section A: For GCH
# read data from file
mydata.GCH_1<-read.table("chr3.GCH.bed",sep="\t")
# set column names
mydata.GCH_1<-cbind(mydata.GCH_1[,1:3],log(mydata.GCH_1[,4]),mydata.GCH_1[,5])
colnames(mydata.GCH_1)<-c("chrom", "start", "end", "Methyl_transform", "Repli")
# remove x-axis cutoff < -10
mydata.GCH_1<-mydata.GCH_1[with(mydata.GCH_1,mydata.GCH_1$Methyl_transform>-10),]
# Estimation of regression parameters
# fit linear models and carry out regression
model_1<-lm(Methyl_transform~Repli,data=mydata.GCH_1)
# Obtain predicted and residual values and draw figures
mydata.GCH_1$predicted <- predict(model_1)
mydata.GCH_1$residuals <- residuals(model_1)
suma<-summary(model_1)
# R2 represents the proportion of variance, in the outcome variable y, that may be predicted by knowing the value of the x variables. 
# An R2 value close to 1 indicates that the model explains a large portion of the variance in the outcome variable.
result_r.squared<-round(suma$r.squared,3)
# the coefficient of RT
a_2<-round(model_1$coefficients[2],3)
a_3<-round(model_1$coefficients[1],3)
n=dim(mydata.GCH_1)[1]
# combine
if(a_3<0){result1<-paste0("Methyl_transform=",a_2,"XRepli",a_3,"; R^2=",result_r.squared,"; n=",n)}else{result1<-paste0("Methyl_transform=",a_2,"XRepli+",a_3,"; R^2=",result_r.squared,"; n=",n)}

p1=ggplot()+
	geom_hex(mapping=aes(x=mydata.GCH_1$Methyl,y=mydata.GCH_1$predicted),color="black")+
	labs(x = "Observed", y = "Model-predicted", title = "chr3: The Scatter Plot of mGpC(JY608-bg)",subtitle=result1)+
	theme(plot.title = element_text(hjust = 0.5,vjust = 1),plot.subtitle = element_text(hjust = 0.5,size=6))+
	expand_limits(x = 0)
	
# Section B: For WCG
# read data from file
mydata.WCG_2<-read.table("chr3.WCG.bed",sep="\t")
# set column names
mydata.WCG_2<-cbind(mydata.WCG_2[,1:3],log(mydata.WCG_2[,4]),mydata.WCG_2[,5])
colnames(mydata.WCG_2)<-c("chrom", "start", "end", "Methyl_transform", "Repli")
# remove x-axis cutoff < -10
mydata.WCG_2<-mydata.WCG_2[with(mydata.WCG_2,mydata.WCG_2$Methyl_transform>-10),]
# Estimation of regression parameters
# fit linear models and carry out regression
model_2<-lm(Methyl_transform~Repli,data=mydata.WCG_2)
# Obtain predicted and residual values and draw figures
mydata.WCG_2$predicted <- predict(model_2)
mydata.WCG_2$residuals <- residuals(model_2)
suma<-summary(model_2)
# R2 represents the proportion of variance, in the outcome variable y, that may be predicted by knowing the value of the x variables. 
# An R2 value close to 1 indicates that the model explains a large portion of the variance in the outcome variable.
result_r.squared<-round(suma$r.squared,3)
# the coefficient of RT
a_2<-round(model_2$coefficients[2],3)
a_3<-round(model_2$coefficients[1],3)
n=dim(mydata.WCG_2)[1]
# combine
if(a_3<0){result2<-paste0("Methyl_transform=",a_2,"XRepli",a_3,"; R^2=",result_r.squared,"; n=",n)}else{result2<-paste0("Methyl_transform=",a_2,"XRepli+",a_3,"; R^2=",result_r.squared,"; n=",n)}

p2=ggplot()+
	geom_hex(mapping=aes(x=mydata.WCG_2$Methyl,y=mydata.WCG_2$predicted),color="black")+
	labs(x = "Observed", y = "Model-predicted", title = "chr3: The Scatter Plot of mCpG(JY608-bg)",subtitle=result2)+
	theme(plot.title = element_text(hjust = 0.5,vjust = 1),plot.subtitle = element_text(hjust = 0.5,size=6))+
	expand_limits(x = 0)
	
# Section A: For GCH
# read data from file
mydata.GCH_3<-read.table("chr4.GCH.bed",sep="\t")
# set column names
mydata.GCH_3<-cbind(mydata.GCH_3[,1:3],log(mydata.GCH_3[,4]),mydata.GCH_3[,5])
colnames(mydata.GCH_3)<-c("chrom", "start", "end", "Methyl_transform", "Repli")
# remove x-axis cutoff < -10
mydata.GCH_3<-mydata.GCH_3[with(mydata.GCH_3,mydata.GCH_3$Methyl_transform>-10),]
# Estimation of regression parameters
# fit linear models and carry out regression
model_3<-lm(Methyl_transform~Repli,data=mydata.GCH_3)
# Obtain predicted and residual values and draw figures
mydata.GCH_3$predicted <- predict(model_3)
mydata.GCH_3$residuals <- residuals(model_3)
suma<-summary(model_3)
# R2 represents the proportion of variance, in the outcome variable y, that may be predicted by knowing the value of the x variables. 
# An R2 value close to 1 indicates that the model explains a large portion of the variance in the outcome variable.
result_r.squared<-round(suma$r.squared,3)
# the coefficient of RT
a_2<-round(model_3$coefficients[2],3)
a_3<-round(model_3$coefficients[1],3)
n=dim(mydata.GCH_3)[1]
# combine
if(a_3<0){result3<-paste0("Methyl_transform=",a_2,"XRepli",a_3,"; R^2=",result_r.squared,"; n=",n)}else{result3<-paste0("Methyl_transform=",a_2,"XRepli+",a_3,"; R^2=",result_r.squared,"; n=",n)}

p3=ggplot()+
	geom_hex(mapping=aes(x=mydata.GCH_3$Methyl,y=mydata.GCH_3$predicted),color="black")+
	labs(x = "Observed", y = "Model-predicted", title = "chr4: The Scatter Plot of mGpC(JY608-bg)",subtitle=result3)+
	theme(plot.title = element_text(hjust = 0.5,vjust = 1),plot.subtitle = element_text(hjust = 0.5,size=6))+
	expand_limits(x = 0)
	
# Section B: For WCG
# read data from file
mydata.WCG_4<-read.table("chr4.WCG.bed",sep="\t")
# set column names
mydata.WCG_4<-cbind(mydata.WCG_4[,1:3],log(mydata.WCG_4[,4]),mydata.WCG_4[,5])
colnames(mydata.WCG_4)<-c("chrom", "start", "end", "Methyl_transform", "Repli")
# remove x-axis cutoff < -10
mydata.WCG_4<-mydata.WCG_4[with(mydata.WCG_4,mydata.WCG_4$Methyl_transform>-10),]
# Estimation of regression parameters
# fit linear models and carry out regression
model_4<-lm(Methyl_transform~Repli,data=mydata.WCG_4)
# Obtain predicted and residual values and draw figures
mydata.WCG_4$predicted <- predict(model_4)
mydata.WCG_4$residuals <- residuals(model_4)
suma<-summary(model_4)
# R2 represents the proportion of variance, in the outcome variable y, that may be predicted by knowing the value of the x variables. 
# An R2 value close to 1 indicates that the model explains a large portion of the variance in the outcome variable.
result_r.squared<-round(suma$r.squared,3)
# the coefficient of RT
a_2<-round(model_4$coefficients[2],3)
a_3<-round(model_4$coefficients[1],3)
n=dim(mydata.WCG_4)[1]
# combine
if(a_3<0){result4<-paste0("Methyl_transform=",a_2,"XRepli",a_3,"; R^2=",result_r.squared,"; n=",n)}else{result4<-paste0("Methyl_transform=",a_2,"XRepli+",a_3,"; R^2=",result_r.squared,"; n=",n)}

p4=ggplot()+
	geom_hex(mapping=aes(x=mydata.WCG_4$Methyl,y=mydata.WCG_4$predicted),color="black")+
	labs(x = "Observed", y = "Model-predicted", title = "chr4: The Scatter Plot of mCpG(JY608-bg)",subtitle=result4)+
	theme(plot.title = element_text(hjust = 0.5,vjust = 1),plot.subtitle = element_text(hjust = 0.5,size=6))+
	expand_limits(x = 0)
	
ggplot2.multiplot(p1,p2,p3,p4,cols=2)		



# chr5+6
# Section A: For GCH
# read data from file
mydata.GCH_1<-read.table("chr5.GCH.bed",sep="\t")
# set column names
mydata.GCH_1<-cbind(mydata.GCH_1[,1:3],log(mydata.GCH_1[,4]),mydata.GCH_1[,5])
colnames(mydata.GCH_1)<-c("chrom", "start", "end", "Methyl_transform", "Repli")
# remove x-axis cutoff < -10
mydata.GCH_1<-mydata.GCH_1[with(mydata.GCH_1,mydata.GCH_1$Methyl_transform>-10),]
# Estimation of regression parameters
# fit linear models and carry out regression
model_1<-lm(Methyl_transform~Repli,data=mydata.GCH_1)
# Obtain predicted and residual values and draw figures
mydata.GCH_1$predicted <- predict(model_1)
mydata.GCH_1$residuals <- residuals(model_1)
suma<-summary(model_1)
# R2 represents the proportion of variance, in the outcome variable y, that may be predicted by knowing the value of the x variables. 
# An R2 value close to 1 indicates that the model explains a large portion of the variance in the outcome variable.
result_r.squared<-round(suma$r.squared,3)
# the coefficient of RT
a_2<-round(model_1$coefficients[2],3)
a_3<-round(model_1$coefficients[1],3)
n=dim(mydata.GCH_1)[1]
# combine
if(a_3<0){result1<-paste0("Methyl_transform=",a_2,"XRepli",a_3,"; R^2=",result_r.squared,"; n=",n)}else{result1<-paste0("Methyl_transform=",a_2,"XRepli+",a_3,"; R^2=",result_r.squared,"; n=",n)}

p1=ggplot()+
	geom_hex(mapping=aes(x=mydata.GCH_1$Methyl,y=mydata.GCH_1$predicted),color="black")+
	labs(x = "Observed", y = "Model-predicted", title = "chr5: The Scatter Plot of mGpC(JY608-bg)",subtitle=result1)+
	theme(plot.title = element_text(hjust = 0.5,vjust = 1),plot.subtitle = element_text(hjust = 0.5,size=6))+
	expand_limits(x = 0)
	
# Section B: For WCG
# read data from file
mydata.WCG_2<-read.table("chr5.WCG.bed",sep="\t")
# set column names
mydata.WCG_2<-cbind(mydata.WCG_2[,1:3],log(mydata.WCG_2[,4]),mydata.WCG_2[,5])
colnames(mydata.WCG_2)<-c("chrom", "start", "end", "Methyl_transform", "Repli")
# remove x-axis cutoff < -10
mydata.WCG_2<-mydata.WCG_2[with(mydata.WCG_2,mydata.WCG_2$Methyl_transform>-10),]
# Estimation of regression parameters
# fit linear models and carry out regression
model_2<-lm(Methyl_transform~Repli,data=mydata.WCG_2)
# Obtain predicted and residual values and draw figures
mydata.WCG_2$predicted <- predict(model_2)
mydata.WCG_2$residuals <- residuals(model_2)
suma<-summary(model_2)
# R2 represents the proportion of variance, in the outcome variable y, that may be predicted by knowing the value of the x variables. 
# An R2 value close to 1 indicates that the model explains a large portion of the variance in the outcome variable.
result_r.squared<-round(suma$r.squared,3)
# the coefficient of RT
a_2<-round(model_2$coefficients[2],3)
a_3<-round(model_2$coefficients[1],3)
n=dim(mydata.WCG_2)[1]
# combine
if(a_3<0){result2<-paste0("Methyl_transform=",a_2,"XRepli",a_3,"; R^2=",result_r.squared,"; n=",n)}else{result2<-paste0("Methyl_transform=",a_2,"XRepli+",a_3,"; R^2=",result_r.squared,"; n=",n)}

p2=ggplot()+
	geom_hex(mapping=aes(x=mydata.WCG_2$Methyl,y=mydata.WCG_2$predicted),color="black")+
	labs(x = "Observed", y = "Model-predicted", title = "chr5: The Scatter Plot of mCpG(JY608-bg)",subtitle=result2)+
	theme(plot.title = element_text(hjust = 0.5,vjust = 1),plot.subtitle = element_text(hjust = 0.5,size=6))+
	expand_limits(x = 0)
	
# Section A: For GCH
# read data from file
mydata.GCH_3<-read.table("chr6.GCH.bed",sep="\t")
# set column names
mydata.GCH_3<-cbind(mydata.GCH_3[,1:3],log(mydata.GCH_3[,4]),mydata.GCH_3[,5])
colnames(mydata.GCH_3)<-c("chrom", "start", "end", "Methyl_transform", "Repli")
# remove x-axis cutoff < -10
mydata.GCH_3<-mydata.GCH_3[with(mydata.GCH_3,mydata.GCH_3$Methyl_transform>-10),]
# Estimation of regression parameters
# fit linear models and carry out regression
model_3<-lm(Methyl_transform~Repli,data=mydata.GCH_3)
# Obtain predicted and residual values and draw figures
mydata.GCH_3$predicted <- predict(model_3)
mydata.GCH_3$residuals <- residuals(model_3)
suma<-summary(model_3)
# R2 represents the proportion of variance, in the outcome variable y, that may be predicted by knowing the value of the x variables. 
# An R2 value close to 1 indicates that the model explains a large portion of the variance in the outcome variable.
result_r.squared<-round(suma$r.squared,3)
# the coefficient of RT
a_2<-round(model_3$coefficients[2],3)
a_3<-round(model_3$coefficients[1],3)
n=dim(mydata.GCH_3)[1]
# combine
if(a_3<0){result3<-paste0("Methyl_transform=",a_2,"XRepli",a_3,"; R^2=",result_r.squared,"; n=",n)}else{result3<-paste0("Methyl_transform=",a_2,"XRepli+",a_3,"; R^2=",result_r.squared,"; n=",n)}

p3=ggplot()+
	geom_hex(mapping=aes(x=mydata.GCH_3$Methyl,y=mydata.GCH_3$predicted),color="black")+
	labs(x = "Observed", y = "Model-predicted", title = "chr6: The Scatter Plot of mGpC(JY608-bg)",subtitle=result3)+
	theme(plot.title = element_text(hjust = 0.5,vjust = 1),plot.subtitle = element_text(hjust = 0.5,size=6))+
	expand_limits(x = 0)
	
# Section B: For WCG
# read data from file
mydata.WCG_4<-read.table("chr6.WCG.bed",sep="\t")
# set column names
mydata.WCG_4<-cbind(mydata.WCG_4[,1:3],log(mydata.WCG_4[,4]),mydata.WCG_4[,5])
colnames(mydata.WCG_4)<-c("chrom", "start", "end", "Methyl_transform", "Repli")
# remove x-axis cutoff < -10
mydata.WCG_4<-mydata.WCG_4[with(mydata.WCG_4,mydata.WCG_4$Methyl_transform>-10),]
# Estimation of regression parameters
# fit linear models and carry out regression
model_4<-lm(Methyl_transform~Repli,data=mydata.WCG_4)
# Obtain predicted and residual values and draw figures
mydata.WCG_4$predicted <- predict(model_4)
mydata.WCG_4$residuals <- residuals(model_4)
suma<-summary(model_4)
# R2 represents the proportion of variance, in the outcome variable y, that may be predicted by knowing the value of the x variables. 
# An R2 value close to 1 indicates that the model explains a large portion of the variance in the outcome variable.
result_r.squared<-round(suma$r.squared,3)
# the coefficient of RT
a_2<-round(model_4$coefficients[2],3)
a_3<-round(model_4$coefficients[1],3)
n=dim(mydata.WCG_4)[1]
# combine
if(a_3<0){result4<-paste0("Methyl_transform=",a_2,"XRepli",a_3,"; R^2=",result_r.squared,"; n=",n)}else{result4<-paste0("Methyl_transform=",a_2,"XRepli+",a_3,"; R^2=",result_r.squared,"; n=",n)}

p4=ggplot()+
	geom_hex(mapping=aes(x=mydata.WCG_4$Methyl,y=mydata.WCG_4$predicted),color="black")+
	labs(x = "Observed", y = "Model-predicted", title = "chr6: The Scatter Plot of mCpG(JY608-bg)",subtitle=result4)+
	theme(plot.title = element_text(hjust = 0.5,vjust = 1),plot.subtitle = element_text(hjust = 0.5,size=6))+
	expand_limits(x = 0)
	
ggplot2.multiplot(p1,p2,p3,p4,cols=2)		



# chr7+8
# Section A: For GCH
# read data from file
mydata.GCH_1<-read.table("chr7.GCH.bed",sep="\t")
# set column names
mydata.GCH_1<-cbind(mydata.GCH_1[,1:3],log(mydata.GCH_1[,4]),mydata.GCH_1[,5])
colnames(mydata.GCH_1)<-c("chrom", "start", "end", "Methyl_transform", "Repli")
# remove x-axis cutoff < -10
mydata.GCH_1<-mydata.GCH_1[with(mydata.GCH_1,mydata.GCH_1$Methyl_transform>-10),]
# Estimation of regression parameters
# fit linear models and carry out regression
model_1<-lm(Methyl_transform~Repli,data=mydata.GCH_1)
# Obtain predicted and residual values and draw figures
mydata.GCH_1$predicted <- predict(model_1)
mydata.GCH_1$residuals <- residuals(model_1)
suma<-summary(model_1)
# R2 represents the proportion of variance, in the outcome variable y, that may be predicted by knowing the value of the x variables. 
# An R2 value close to 1 indicates that the model explains a large portion of the variance in the outcome variable.
result_r.squared<-round(suma$r.squared,3)
# the coefficient of RT
a_2<-round(model_1$coefficients[2],3)
a_3<-round(model_1$coefficients[1],3)
n=dim(mydata.GCH_1)[1]
# combine
if(a_3<0){result1<-paste0("Methyl_transform=",a_2,"XRepli",a_3,"; R^2=",result_r.squared,"; n=",n)}else{result1<-paste0("Methyl_transform=",a_2,"XRepli+",a_3,"; R^2=",result_r.squared,"; n=",n)}

p1=ggplot()+
	geom_hex(mapping=aes(x=mydata.GCH_1$Methyl,y=mydata.GCH_1$predicted),color="black")+
	labs(x = "Observed", y = "Model-predicted", title = "chr7: The Scatter Plot of mGpC(JY608-bg)",subtitle=result1)+
	theme(plot.title = element_text(hjust = 0.5,vjust = 1),plot.subtitle = element_text(hjust = 0.5,size=6))+
	expand_limits(x = 0)
	
# Section B: For WCG
# read data from file
mydata.WCG_2<-read.table("chr7.WCG.bed",sep="\t")
# set column names
mydata.WCG_2<-cbind(mydata.WCG_2[,1:3],log(mydata.WCG_2[,4]),mydata.WCG_2[,5])
colnames(mydata.WCG_2)<-c("chrom", "start", "end", "Methyl_transform", "Repli")
# remove x-axis cutoff < -10
mydata.WCG_2<-mydata.WCG_2[with(mydata.WCG_2,mydata.WCG_2$Methyl_transform>-10),]
# Estimation of regression parameters
# fit linear models and carry out regression
model_2<-lm(Methyl_transform~Repli,data=mydata.WCG_2)
# Obtain predicted and residual values and draw figures
mydata.WCG_2$predicted <- predict(model_2)
mydata.WCG_2$residuals <- residuals(model_2)
suma<-summary(model_2)
# R2 represents the proportion of variance, in the outcome variable y, that may be predicted by knowing the value of the x variables. 
# An R2 value close to 1 indicates that the model explains a large portion of the variance in the outcome variable.
result_r.squared<-round(suma$r.squared,3)
# the coefficient of RT
a_2<-round(model_2$coefficients[2],3)
a_3<-round(model_2$coefficients[1],3)
n=dim(mydata.WCG_2)[1]
# combine
if(a_3<0){result2<-paste0("Methyl_transform=",a_2,"XRepli",a_3,"; R^2=",result_r.squared,"; n=",n)}else{result2<-paste0("Methyl_transform=",a_2,"XRepli+",a_3,"; R^2=",result_r.squared,"; n=",n)}

p2=ggplot()+
	geom_hex(mapping=aes(x=mydata.WCG_2$Methyl,y=mydata.WCG_2$predicted),color="black")+
	labs(x = "Observed", y = "Model-predicted", title = "chr7: The Scatter Plot of mCpG(JY608-bg)",subtitle=result2)+
	theme(plot.title = element_text(hjust = 0.5,vjust = 1),plot.subtitle = element_text(hjust = 0.5,size=6))+
	expand_limits(x = 0)
	
# Section A: For GCH
# read data from file
mydata.GCH_3<-read.table("chr8.GCH.bed",sep="\t")
# set column names
mydata.GCH_3<-cbind(mydata.GCH_3[,1:3],log(mydata.GCH_3[,4]),mydata.GCH_3[,5])
colnames(mydata.GCH_3)<-c("chrom", "start", "end", "Methyl_transform", "Repli")
# remove x-axis cutoff < -10
mydata.GCH_3<-mydata.GCH_3[with(mydata.GCH_3,mydata.GCH_3$Methyl_transform>-10),]
# Estimation of regression parameters
# fit linear models and carry out regression
model_3<-lm(Methyl_transform~Repli,data=mydata.GCH_3)
# Obtain predicted and residual values and draw figures
mydata.GCH_3$predicted <- predict(model_3)
mydata.GCH_3$residuals <- residuals(model_3)
suma<-summary(model_3)
# R2 represents the proportion of variance, in the outcome variable y, that may be predicted by knowing the value of the x variables. 
# An R2 value close to 1 indicates that the model explains a large portion of the variance in the outcome variable.
result_r.squared<-round(suma$r.squared,3)
# the coefficient of RT
a_2<-round(model_3$coefficients[2],3)
a_3<-round(model_3$coefficients[1],3)
n=dim(mydata.GCH_3)[1]
# combine
if(a_3<0){result3<-paste0("Methyl_transform=",a_2,"XRepli",a_3,"; R^2=",result_r.squared,"; n=",n)}else{result3<-paste0("Methyl_transform=",a_2,"XRepli+",a_3,"; R^2=",result_r.squared,"; n=",n)}

p3=ggplot()+
	geom_hex(mapping=aes(x=mydata.GCH_3$Methyl,y=mydata.GCH_3$predicted),color="black")+
	labs(x = "Observed", y = "Model-predicted", title = "chr8: The Scatter Plot of mGpC(JY608-bg)",subtitle=result3)+
	theme(plot.title = element_text(hjust = 0.5,vjust = 1),plot.subtitle = element_text(hjust = 0.5,size=6))+
	expand_limits(x = 0)
	
# Section B: For WCG
# read data from file
mydata.WCG_4<-read.table("chr8.WCG.bed",sep="\t")
# set column names
mydata.WCG_4<-cbind(mydata.WCG_4[,1:3],log(mydata.WCG_4[,4]),mydata.WCG_4[,5])
colnames(mydata.WCG_4)<-c("chrom", "start", "end", "Methyl_transform", "Repli")
# remove x-axis cutoff < -10
mydata.WCG_4<-mydata.WCG_4[with(mydata.WCG_4,mydata.WCG_4$Methyl_transform>-10),]
# Estimation of regression parameters
# fit linear models and carry out regression
model_4<-lm(Methyl_transform~Repli,data=mydata.WCG_4)
# Obtain predicted and residual values and draw figures
mydata.WCG_4$predicted <- predict(model_4)
mydata.WCG_4$residuals <- residuals(model_4)
suma<-summary(model_4)
# R2 represents the proportion of variance, in the outcome variable y, that may be predicted by knowing the value of the x variables. 
# An R2 value close to 1 indicates that the model explains a large portion of the variance in the outcome variable.
result_r.squared<-round(suma$r.squared,3)
# the coefficient of RT
a_2<-round(model_4$coefficients[2],3)
a_3<-round(model_4$coefficients[1],3)
n=dim(mydata.WCG_4)[1]
# combine
if(a_3<0){result4<-paste0("Methyl_transform=",a_2,"XRepli",a_3,"; R^2=",result_r.squared,"; n=",n)}else{result4<-paste0("Methyl_transform=",a_2,"XRepli+",a_3,"; R^2=",result_r.squared,"; n=",n)}

p4=ggplot()+
	geom_hex(mapping=aes(x=mydata.WCG_4$Methyl,y=mydata.WCG_4$predicted),color="black")+
	labs(x = "Observed", y = "Model-predicted", title = "chr8: The Scatter Plot of mCpG(JY608-bg)",subtitle=result4)+
	theme(plot.title = element_text(hjust = 0.5,vjust = 1),plot.subtitle = element_text(hjust = 0.5,size=6))+
	expand_limits(x = 0)
	
ggplot2.multiplot(p1,p2,p3,p4,cols=2)		



# chr9+10
# Section A: For GCH
# read data from file
mydata.GCH_1<-read.table("chr9.GCH.bed",sep="\t")
# set column names
mydata.GCH_1<-cbind(mydata.GCH_1[,1:3],log(mydata.GCH_1[,4]),mydata.GCH_1[,5])
colnames(mydata.GCH_1)<-c("chrom", "start", "end", "Methyl_transform", "Repli")
# remove x-axis cutoff < -10
mydata.GCH_1<-mydata.GCH_1[with(mydata.GCH_1,mydata.GCH_1$Methyl_transform>-10),]
# Estimation of regression parameters
# fit linear models and carry out regression
model_1<-lm(Methyl_transform~Repli,data=mydata.GCH_1)
# Obtain predicted and residual values and draw figures
mydata.GCH_1$predicted <- predict(model_1)
mydata.GCH_1$residuals <- residuals(model_1)
suma<-summary(model_1)
# R2 represents the proportion of variance, in the outcome variable y, that may be predicted by knowing the value of the x variables. 
# An R2 value close to 1 indicates that the model explains a large portion of the variance in the outcome variable.
result_r.squared<-round(suma$r.squared,3)
# the coefficient of RT
a_2<-round(model_1$coefficients[2],3)
a_3<-round(model_1$coefficients[1],3)
n=dim(mydata.GCH_1)[1]
# combine
if(a_3<0){result1<-paste0("Methyl_transform=",a_2,"XRepli",a_3,"; R^2=",result_r.squared,"; n=",n)}else{result1<-paste0("Methyl_transform=",a_2,"XRepli+",a_3,"; R^2=",result_r.squared,"; n=",n)}

p1=ggplot()+
	geom_hex(mapping=aes(x=mydata.GCH_1$Methyl,y=mydata.GCH_1$predicted),color="black")+
	labs(x = "Observed", y = "Model-predicted", title = "chr9: The Scatter Plot of mGpC(JY608-bg)",subtitle=result1)+
	theme(plot.title = element_text(hjust = 0.5,vjust = 1),plot.subtitle = element_text(hjust = 0.5,size=6))+
	expand_limits(x = 0)
	
# Section B: For WCG
# read data from file
mydata.WCG_2<-read.table("chr9.WCG.bed",sep="\t")
# set column names
mydata.WCG_2<-cbind(mydata.WCG_2[,1:3],log(mydata.WCG_2[,4]),mydata.WCG_2[,5])
colnames(mydata.WCG_2)<-c("chrom", "start", "end", "Methyl_transform", "Repli")
# remove x-axis cutoff < -10
mydata.WCG_2<-mydata.WCG_2[with(mydata.WCG_2,mydata.WCG_2$Methyl_transform>-10),]
# Estimation of regression parameters
# fit linear models and carry out regression
model_2<-lm(Methyl_transform~Repli,data=mydata.WCG_2)
# Obtain predicted and residual values and draw figures
mydata.WCG_2$predicted <- predict(model_2)
mydata.WCG_2$residuals <- residuals(model_2)
suma<-summary(model_2)
# R2 represents the proportion of variance, in the outcome variable y, that may be predicted by knowing the value of the x variables. 
# An R2 value close to 1 indicates that the model explains a large portion of the variance in the outcome variable.
result_r.squared<-round(suma$r.squared,3)
# the coefficient of RT
a_2<-round(model_2$coefficients[2],3)
a_3<-round(model_2$coefficients[1],3)
n=dim(mydata.WCG_2)[1]
# combine
if(a_3<0){result2<-paste0("Methyl_transform=",a_2,"XRepli",a_3,"; R^2=",result_r.squared,"; n=",n)}else{result2<-paste0("Methyl_transform=",a_2,"XRepli+",a_3,"; R^2=",result_r.squared,"; n=",n)}

p2=ggplot()+
	geom_hex(mapping=aes(x=mydata.WCG_2$Methyl,y=mydata.WCG_2$predicted),color="black")+
	labs(x = "Observed", y = "Model-predicted", title = "chr9: The Scatter Plot of mCpG(JY608-bg)",subtitle=result2)+
	theme(plot.title = element_text(hjust = 0.5,vjust = 1),plot.subtitle = element_text(hjust = 0.5,size=6))+
	expand_limits(x = 0)
	
# Section A: For GCH
# read data from file
mydata.GCH_3<-read.table("chr10.GCH.bed",sep="\t")
# set column names
mydata.GCH_3<-cbind(mydata.GCH_3[,1:3],log(mydata.GCH_3[,4]),mydata.GCH_3[,5])
colnames(mydata.GCH_3)<-c("chrom", "start", "end", "Methyl_transform", "Repli")
# remove x-axis cutoff < -10
mydata.GCH_3<-mydata.GCH_3[with(mydata.GCH_3,mydata.GCH_3$Methyl_transform>-10),]
# Estimation of regression parameters
# fit linear models and carry out regression
model_3<-lm(Methyl_transform~Repli,data=mydata.GCH_3)
# Obtain predicted and residual values and draw figures
mydata.GCH_3$predicted <- predict(model_3)
mydata.GCH_3$residuals <- residuals(model_3)
suma<-summary(model_3)
# R2 represents the proportion of variance, in the outcome variable y, that may be predicted by knowing the value of the x variables. 
# An R2 value close to 1 indicates that the model explains a large portion of the variance in the outcome variable.
result_r.squared<-round(suma$r.squared,3)
# the coefficient of RT
a_2<-round(model_3$coefficients[2],3)
a_3<-round(model_3$coefficients[1],3)
n=dim(mydata.GCH_3)[1]
# combine
if(a_3<0){result3<-paste0("Methyl_transform=",a_2,"XRepli",a_3,"; R^2=",result_r.squared,"; n=",n)}else{result3<-paste0("Methyl_transform=",a_2,"XRepli+",a_3,"; R^2=",result_r.squared,"; n=",n)}

p3=ggplot()+
	geom_hex(mapping=aes(x=mydata.GCH_3$Methyl,y=mydata.GCH_3$predicted),color="black")+
	labs(x = "Observed", y = "Model-predicted", title = "chr10: The Scatter Plot of mGpC(JY608-bg)",subtitle=result3)+
	theme(plot.title = element_text(hjust = 0.5,vjust = 1),plot.subtitle = element_text(hjust = 0.5,size=6))+
	expand_limits(x = 0)
	
# Section B: For WCG
# read data from file
mydata.WCG_4<-read.table("chr10.WCG.bed",sep="\t")
# set column names
mydata.WCG_4<-cbind(mydata.WCG_4[,1:3],log(mydata.WCG_4[,4]),mydata.WCG_4[,5])
colnames(mydata.WCG_4)<-c("chrom", "start", "end", "Methyl_transform", "Repli")
# remove x-axis cutoff < -10
mydata.WCG_4<-mydata.WCG_4[with(mydata.WCG_4,mydata.WCG_4$Methyl_transform>-10),]
# Estimation of regression parameters
# fit linear models and carry out regression
model_4<-lm(Methyl_transform~Repli,data=mydata.WCG_4)
# Obtain predicted and residual values and draw figures
mydata.WCG_4$predicted <- predict(model_4)
mydata.WCG_4$residuals <- residuals(model_4)
suma<-summary(model_4)
# R2 represents the proportion of variance, in the outcome variable y, that may be predicted by knowing the value of the x variables. 
# An R2 value close to 1 indicates that the model explains a large portion of the variance in the outcome variable.
result_r.squared<-round(suma$r.squared,3)
# the coefficient of RT
a_2<-round(model_4$coefficients[2],3)
a_3<-round(model_4$coefficients[1],3)
n=dim(mydata.WCG_4)[1]
# combine
if(a_3<0){result4<-paste0("Methyl_transform=",a_2,"XRepli",a_3,"; R^2=",result_r.squared,"; n=",n)}else{result4<-paste0("Methyl_transform=",a_2,"XRepli+",a_3,"; R^2=",result_r.squared,"; n=",n)}

p4=ggplot()+
	geom_hex(mapping=aes(x=mydata.WCG_4$Methyl,y=mydata.WCG_4$predicted),color="black")+
	labs(x = "Observed", y = "Model-predicted", title = "chr10: The Scatter Plot of mCpG(JY608-bg)",subtitle=result4)+
	theme(plot.title = element_text(hjust = 0.5,vjust = 1),plot.subtitle = element_text(hjust = 0.5,size=6))+
	expand_limits(x = 0)
	
ggplot2.multiplot(p1,p2,p3,p4,cols=2)		




# chr11+12
# Section A: For GCH
# read data from file
mydata.GCH_1<-read.table("chr11.GCH.bed",sep="\t")
# set column names
mydata.GCH_1<-cbind(mydata.GCH_1[,1:3],log(mydata.GCH_1[,4]),mydata.GCH_1[,5])
colnames(mydata.GCH_1)<-c("chrom", "start", "end", "Methyl_transform", "Repli")
# remove x-axis cutoff < -10
mydata.GCH_1<-mydata.GCH_1[with(mydata.GCH_1,mydata.GCH_1$Methyl_transform>-10),]
# Estimation of regression parameters
# fit linear models and carry out regression
model_1<-lm(Methyl_transform~Repli,data=mydata.GCH_1)
# Obtain predicted and residual values and draw figures
mydata.GCH_1$predicted <- predict(model_1)
mydata.GCH_1$residuals <- residuals(model_1)
suma<-summary(model_1)
# R2 represents the proportion of variance, in the outcome variable y, that may be predicted by knowing the value of the x variables. 
# An R2 value close to 1 indicates that the model explains a large portion of the variance in the outcome variable.
result_r.squared<-round(suma$r.squared,3)
# the coefficient of RT
a_2<-round(model_1$coefficients[2],3)
a_3<-round(model_1$coefficients[1],3)
n=dim(mydata.GCH_1)[1]
# combine
if(a_3<0){result1<-paste0("Methyl_transform=",a_2,"XRepli",a_3,"; R^2=",result_r.squared,"; n=",n)}else{result1<-paste0("Methyl_transform=",a_2,"XRepli+",a_3,"; R^2=",result_r.squared,"; n=",n)}

p1=ggplot()+
	geom_hex(mapping=aes(x=mydata.GCH_1$Methyl,y=mydata.GCH_1$predicted),color="black")+
	labs(x = "Observed", y = "Model-predicted", title = "chr11: The Scatter Plot of mGpC(JY608-bg)",subtitle=result1)+
	theme(plot.title = element_text(hjust = 0.5,vjust = 1),plot.subtitle = element_text(hjust = 0.5,size=6))+
	expand_limits(x = 0)
	
# Section B: For WCG
# read data from file
mydata.WCG_2<-read.table("chr11.WCG.bed",sep="\t")
# set column names
mydata.WCG_2<-cbind(mydata.WCG_2[,1:3],log(mydata.WCG_2[,4]),mydata.WCG_2[,5])
colnames(mydata.WCG_2)<-c("chrom", "start", "end", "Methyl_transform", "Repli")
# remove x-axis cutoff < -10
mydata.WCG_2<-mydata.WCG_2[with(mydata.WCG_2,mydata.WCG_2$Methyl_transform>-10),]
# Estimation of regression parameters
# fit linear models and carry out regression
model_2<-lm(Methyl_transform~Repli,data=mydata.WCG_2)
# Obtain predicted and residual values and draw figures
mydata.WCG_2$predicted <- predict(model_2)
mydata.WCG_2$residuals <- residuals(model_2)
suma<-summary(model_2)
# R2 represents the proportion of variance, in the outcome variable y, that may be predicted by knowing the value of the x variables. 
# An R2 value close to 1 indicates that the model explains a large portion of the variance in the outcome variable.
result_r.squared<-round(suma$r.squared,3)
# the coefficient of RT
a_2<-round(model_2$coefficients[2],3)
a_3<-round(model_2$coefficients[1],3)
n=dim(mydata.WCG_2)[1]
# combine
if(a_3<0){result2<-paste0("Methyl_transform=",a_2,"XRepli",a_3,"; R^2=",result_r.squared,"; n=",n)}else{result2<-paste0("Methyl_transform=",a_2,"XRepli+",a_3,"; R^2=",result_r.squared,"; n=",n)}

p2=ggplot()+
	geom_hex(mapping=aes(x=mydata.WCG_2$Methyl,y=mydata.WCG_2$predicted),color="black")+
	labs(x = "Observed", y = "Model-predicted", title = "chr11: The Scatter Plot of mCpG(JY608-bg)",subtitle=result2)+
	theme(plot.title = element_text(hjust = 0.5,vjust = 1),plot.subtitle = element_text(hjust = 0.5,size=6))+
	expand_limits(x = 0)
	
# Section A: For GCH
# read data from file
mydata.GCH_3<-read.table("chr12.GCH.bed",sep="\t")
# set column names
mydata.GCH_3<-cbind(mydata.GCH_3[,1:3],log(mydata.GCH_3[,4]),mydata.GCH_3[,5])
colnames(mydata.GCH_3)<-c("chrom", "start", "end", "Methyl_transform", "Repli")
# remove x-axis cutoff < -10
mydata.GCH_3<-mydata.GCH_3[with(mydata.GCH_3,mydata.GCH_3$Methyl_transform>-10),]
# Estimation of regression parameters
# fit linear models and carry out regression
model_3<-lm(Methyl_transform~Repli,data=mydata.GCH_3)
# Obtain predicted and residual values and draw figures
mydata.GCH_3$predicted <- predict(model_3)
mydata.GCH_3$residuals <- residuals(model_3)
suma<-summary(model_3)
# R2 represents the proportion of variance, in the outcome variable y, that may be predicted by knowing the value of the x variables. 
# An R2 value close to 1 indicates that the model explains a large portion of the variance in the outcome variable.
result_r.squared<-round(suma$r.squared,3)
# the coefficient of RT
a_2<-round(model_3$coefficients[2],3)
a_3<-round(model_3$coefficients[1],3)
n=dim(mydata.GCH_3)[1]
# combine
if(a_3<0){result3<-paste0("Methyl_transform=",a_2,"XRepli",a_3,"; R^2=",result_r.squared,"; n=",n)}else{result3<-paste0("Methyl_transform=",a_2,"XRepli+",a_3,"; R^2=",result_r.squared,"; n=",n)}

p3=ggplot()+
	geom_hex(mapping=aes(x=mydata.GCH_3$Methyl,y=mydata.GCH_3$predicted),color="black")+
	labs(x = "Observed", y = "Model-predicted", title = "chr12: The Scatter Plot of mGpC(JY608-bg)",subtitle=result3)+
	theme(plot.title = element_text(hjust = 0.5,vjust = 1),plot.subtitle = element_text(hjust = 0.5,size=6))+
	expand_limits(x = 0)
	
# Section B: For WCG
# read data from file
mydata.WCG_4<-read.table("chr12.WCG.bed",sep="\t")
# set column names
mydata.WCG_4<-cbind(mydata.WCG_4[,1:3],log(mydata.WCG_4[,4]),mydata.WCG_4[,5])
colnames(mydata.WCG_4)<-c("chrom", "start", "end", "Methyl_transform", "Repli")
# remove x-axis cutoff < -10
mydata.WCG_4<-mydata.WCG_4[with(mydata.WCG_4,mydata.WCG_4$Methyl_transform>-10),]
# Estimation of regression parameters
# fit linear models and carry out regression
model_4<-lm(Methyl_transform~Repli,data=mydata.WCG_4)
# Obtain predicted and residual values and draw figures
mydata.WCG_4$predicted <- predict(model_4)
mydata.WCG_4$residuals <- residuals(model_4)
suma<-summary(model_4)
# R2 represents the proportion of variance, in the outcome variable y, that may be predicted by knowing the value of the x variables. 
# An R2 value close to 1 indicates that the model explains a large portion of the variance in the outcome variable.
result_r.squared<-round(suma$r.squared,3)
# the coefficient of RT
a_2<-round(model_4$coefficients[2],3)
a_3<-round(model_4$coefficients[1],3)
n=dim(mydata.WCG_4)[1]
# combine
if(a_3<0){result4<-paste0("Methyl_transform=",a_2,"XRepli",a_3,"; R^2=",result_r.squared,"; n=",n)}else{result4<-paste0("Methyl_transform=",a_2,"XRepli+",a_3,"; R^2=",result_r.squared,"; n=",n)}

p4=ggplot()+
	geom_hex(mapping=aes(x=mydata.WCG_4$Methyl,y=mydata.WCG_4$predicted),color="black")+
	labs(x = "Observed", y = "Model-predicted", title = "chr12: The Scatter Plot of mCpG(JY608-bg)",subtitle=result4)+
	theme(plot.title = element_text(hjust = 0.5,vjust = 1),plot.subtitle = element_text(hjust = 0.5,size=6))+
	expand_limits(x = 0)
	
ggplot2.multiplot(p1,p2,p3,p4,cols=2)		



# chr13+14
# Section A: For GCH
# read data from file
mydata.GCH_1<-read.table("chr13.GCH.bed",sep="\t")
# set column names
mydata.GCH_1<-cbind(mydata.GCH_1[,1:3],log(mydata.GCH_1[,4]),mydata.GCH_1[,5])
colnames(mydata.GCH_1)<-c("chrom", "start", "end", "Methyl_transform", "Repli")
# remove x-axis cutoff < -10
mydata.GCH_1<-mydata.GCH_1[with(mydata.GCH_1,mydata.GCH_1$Methyl_transform>-10),]
# Estimation of regression parameters
# fit linear models and carry out regression
model_1<-lm(Methyl_transform~Repli,data=mydata.GCH_1)
# Obtain predicted and residual values and draw figures
mydata.GCH_1$predicted <- predict(model_1)
mydata.GCH_1$residuals <- residuals(model_1)
suma<-summary(model_1)
# R2 represents the proportion of variance, in the outcome variable y, that may be predicted by knowing the value of the x variables. 
# An R2 value close to 1 indicates that the model explains a large portion of the variance in the outcome variable.
result_r.squared<-round(suma$r.squared,3)
# the coefficient of RT
a_2<-round(model_1$coefficients[2],3)
a_3<-round(model_1$coefficients[1],3)
n=dim(mydata.GCH_1)[1]
# combine
if(a_3<0){result1<-paste0("Methyl_transform=",a_2,"XRepli",a_3,"; R^2=",result_r.squared,"; n=",n)}else{result1<-paste0("Methyl_transform=",a_2,"XRepli+",a_3,"; R^2=",result_r.squared,"; n=",n)}

p1=ggplot()+
	geom_hex(mapping=aes(x=mydata.GCH_1$Methyl,y=mydata.GCH_1$predicted),color="black")+
	labs(x = "Observed", y = "Model-predicted", title = "chr13: The Scatter Plot of mGpC(JY608-bg)",subtitle=result1)+
	theme(plot.title = element_text(hjust = 0.5,vjust = 1),plot.subtitle = element_text(hjust = 0.5,size=6))+
	expand_limits(x = 0)
	
# Section B: For WCG
# read data from file
mydata.WCG_2<-read.table("chr13.WCG.bed",sep="\t")
# set column names
mydata.WCG_2<-cbind(mydata.WCG_2[,1:3],log(mydata.WCG_2[,4]),mydata.WCG_2[,5])
colnames(mydata.WCG_2)<-c("chrom", "start", "end", "Methyl_transform", "Repli")
# remove x-axis cutoff < -10
mydata.WCG_2<-mydata.WCG_2[with(mydata.WCG_2,mydata.WCG_2$Methyl_transform>-10),]
# Estimation of regression parameters
# fit linear models and carry out regression
model_2<-lm(Methyl_transform~Repli,data=mydata.WCG_2)
# Obtain predicted and residual values and draw figures
mydata.WCG_2$predicted <- predict(model_2)
mydata.WCG_2$residuals <- residuals(model_2)
suma<-summary(model_2)
# R2 represents the proportion of variance, in the outcome variable y, that may be predicted by knowing the value of the x variables. 
# An R2 value close to 1 indicates that the model explains a large portion of the variance in the outcome variable.
result_r.squared<-round(suma$r.squared,3)
# the coefficient of RT
a_2<-round(model_2$coefficients[2],3)
a_3<-round(model_2$coefficients[1],3)
n=dim(mydata.WCG_2)[1]
# combine
if(a_3<0){result2<-paste0("Methyl_transform=",a_2,"XRepli",a_3,"; R^2=",result_r.squared,"; n=",n)}else{result2<-paste0("Methyl_transform=",a_2,"XRepli+",a_3,"; R^2=",result_r.squared,"; n=",n)}

p2=ggplot()+
	geom_hex(mapping=aes(x=mydata.WCG_2$Methyl,y=mydata.WCG_2$predicted),color="black")+
	labs(x = "Observed", y = "Model-predicted", title = "chr13: The Scatter Plot of mCpG(JY608-bg)",subtitle=result2)+
	theme(plot.title = element_text(hjust = 0.5,vjust = 1),plot.subtitle = element_text(hjust = 0.5,size=6))+
	expand_limits(x = 0)
	
# Section A: For GCH
# read data from file
mydata.GCH_3<-read.table("chr14.GCH.bed",sep="\t")
# set column names
mydata.GCH_3<-cbind(mydata.GCH_3[,1:3],log(mydata.GCH_3[,4]),mydata.GCH_3[,5])
colnames(mydata.GCH_3)<-c("chrom", "start", "end", "Methyl_transform", "Repli")
# remove x-axis cutoff < -10
mydata.GCH_3<-mydata.GCH_3[with(mydata.GCH_3,mydata.GCH_3$Methyl_transform>-10),]
# Estimation of regression parameters
# fit linear models and carry out regression
model_3<-lm(Methyl_transform~Repli,data=mydata.GCH_3)
# Obtain predicted and residual values and draw figures
mydata.GCH_3$predicted <- predict(model_3)
mydata.GCH_3$residuals <- residuals(model_3)
suma<-summary(model_3)
# R2 represents the proportion of variance, in the outcome variable y, that may be predicted by knowing the value of the x variables. 
# An R2 value close to 1 indicates that the model explains a large portion of the variance in the outcome variable.
result_r.squared<-round(suma$r.squared,3)
# the coefficient of RT
a_2<-round(model_3$coefficients[2],3)
a_3<-round(model_3$coefficients[1],3)
n=dim(mydata.GCH_3)[1]
# combine
if(a_3<0){result3<-paste0("Methyl_transform=",a_2,"XRepli",a_3,"; R^2=",result_r.squared,"; n=",n)}else{result3<-paste0("Methyl_transform=",a_2,"XRepli+",a_3,"; R^2=",result_r.squared,"; n=",n)}

p3=ggplot()+
	geom_hex(mapping=aes(x=mydata.GCH_3$Methyl,y=mydata.GCH_3$predicted),color="black")+
	labs(x = "Observed", y = "Model-predicted", title = "chr14: The Scatter Plot of mGpC(JY608-bg)",subtitle=result3)+
	theme(plot.title = element_text(hjust = 0.5,vjust = 1),plot.subtitle = element_text(hjust = 0.5,size=6))+
	expand_limits(x = 0)
	
# Section B: For WCG
# read data from file
mydata.WCG_4<-read.table("chr14.WCG.bed",sep="\t")
# set column names
mydata.WCG_4<-cbind(mydata.WCG_4[,1:3],log(mydata.WCG_4[,4]),mydata.WCG_4[,5])
colnames(mydata.WCG_4)<-c("chrom", "start", "end", "Methyl_transform", "Repli")
# remove x-axis cutoff < -10
mydata.WCG_4<-mydata.WCG_4[with(mydata.WCG_4,mydata.WCG_4$Methyl_transform>-10),]
# Estimation of regression parameters
# fit linear models and carry out regression
model_4<-lm(Methyl_transform~Repli,data=mydata.WCG_4)
# Obtain predicted and residual values and draw figures
mydata.WCG_4$predicted <- predict(model_4)
mydata.WCG_4$residuals <- residuals(model_4)
suma<-summary(model_4)
# R2 represents the proportion of variance, in the outcome variable y, that may be predicted by knowing the value of the x variables. 
# An R2 value close to 1 indicates that the model explains a large portion of the variance in the outcome variable.
result_r.squared<-round(suma$r.squared,3)
# the coefficient of RT
a_2<-round(model_4$coefficients[2],3)
a_3<-round(model_4$coefficients[1],3)
n=dim(mydata.WCG_4)[1]
# combine
if(a_3<0){result4<-paste0("Methyl_transform=",a_2,"XRepli",a_3,"; R^2=",result_r.squared,"; n=",n)}else{result4<-paste0("Methyl_transform=",a_2,"XRepli+",a_3,"; R^2=",result_r.squared,"; n=",n)}

p4=ggplot()+
	geom_hex(mapping=aes(x=mydata.WCG_4$Methyl,y=mydata.WCG_4$predicted),color="black")+
	labs(x = "Observed", y = "Model-predicted", title = "chr14: The Scatter Plot of mCpG(JY608-bg)",subtitle=result4)+
	theme(plot.title = element_text(hjust = 0.5,vjust = 1),plot.subtitle = element_text(hjust = 0.5,size=6))+
	expand_limits(x = 0)
	
ggplot2.multiplot(p1,p2,p3,p4,cols=2)		



# chr15+16
# Section A: For GCH
# read data from file
mydata.GCH_1<-read.table("chr15.GCH.bed",sep="\t")
# set column names
mydata.GCH_1<-cbind(mydata.GCH_1[,1:3],log(mydata.GCH_1[,4]),mydata.GCH_1[,5])
colnames(mydata.GCH_1)<-c("chrom", "start", "end", "Methyl_transform", "Repli")
# remove x-axis cutoff < -10
mydata.GCH_1<-mydata.GCH_1[with(mydata.GCH_1,mydata.GCH_1$Methyl_transform>-10),]
# Estimation of regression parameters
# fit linear models and carry out regression
model_1<-lm(Methyl_transform~Repli,data=mydata.GCH_1)
# Obtain predicted and residual values and draw figures
mydata.GCH_1$predicted <- predict(model_1)
mydata.GCH_1$residuals <- residuals(model_1)
suma<-summary(model_1)
# R2 represents the proportion of variance, in the outcome variable y, that may be predicted by knowing the value of the x variables. 
# An R2 value close to 1 indicates that the model explains a large portion of the variance in the outcome variable.
result_r.squared<-round(suma$r.squared,3)
# the coefficient of RT
a_2<-round(model_1$coefficients[2],3)
a_3<-round(model_1$coefficients[1],3)
n=dim(mydata.GCH_1)[1]
# combine
if(a_3<0){result1<-paste0("Methyl_transform=",a_2,"XRepli",a_3,"; R^2=",result_r.squared,"; n=",n)}else{result1<-paste0("Methyl_transform=",a_2,"XRepli+",a_3,"; R^2=",result_r.squared,"; n=",n)}

p1=ggplot()+
	geom_hex(mapping=aes(x=mydata.GCH_1$Methyl,y=mydata.GCH_1$predicted),color="black")+
	labs(x = "Observed", y = "Model-predicted", title = "chr15: The Scatter Plot of mGpC(JY608-bg)",subtitle=result1)+
	theme(plot.title = element_text(hjust = 0.5,vjust = 1),plot.subtitle = element_text(hjust = 0.5,size=6))+
	expand_limits(x = 0)
	
# Section B: For WCG
# read data from file
mydata.WCG_2<-read.table("chr15.WCG.bed",sep="\t")
# set column names
mydata.WCG_2<-cbind(mydata.WCG_2[,1:3],log(mydata.WCG_2[,4]),mydata.WCG_2[,5])
colnames(mydata.WCG_2)<-c("chrom", "start", "end", "Methyl_transform", "Repli")
# remove x-axis cutoff < -10
mydata.WCG_2<-mydata.WCG_2[with(mydata.WCG_2,mydata.WCG_2$Methyl_transform>-10),]
# Estimation of regression parameters
# fit linear models and carry out regression
model_2<-lm(Methyl_transform~Repli,data=mydata.WCG_2)
# Obtain predicted and residual values and draw figures
mydata.WCG_2$predicted <- predict(model_2)
mydata.WCG_2$residuals <- residuals(model_2)
suma<-summary(model_2)
# R2 represents the proportion of variance, in the outcome variable y, that may be predicted by knowing the value of the x variables. 
# An R2 value close to 1 indicates that the model explains a large portion of the variance in the outcome variable.
result_r.squared<-round(suma$r.squared,3)
# the coefficient of RT
a_2<-round(model_2$coefficients[2],3)
a_3<-round(model_2$coefficients[1],3)
n=dim(mydata.WCG_2)[1]
# combine
if(a_3<0){result2<-paste0("Methyl_transform=",a_2,"XRepli",a_3,"; R^2=",result_r.squared,"; n=",n)}else{result2<-paste0("Methyl_transform=",a_2,"XRepli+",a_3,"; R^2=",result_r.squared,"; n=",n)}

p2=ggplot()+
	geom_hex(mapping=aes(x=mydata.WCG_2$Methyl,y=mydata.WCG_2$predicted),color="black")+
	labs(x = "Observed", y = "Model-predicted", title = "chr15: The Scatter Plot of mCpG(JY608-bg)",subtitle=result2)+
	theme(plot.title = element_text(hjust = 0.5,vjust = 1),plot.subtitle = element_text(hjust = 0.5,size=6))+
	expand_limits(x = 0)
	
# Section A: For GCH
# read data from file
mydata.GCH_3<-read.table("chr16.GCH.bed",sep="\t")
# set column names
mydata.GCH_3<-cbind(mydata.GCH_3[,1:3],log(mydata.GCH_3[,4]),mydata.GCH_3[,5])
colnames(mydata.GCH_3)<-c("chrom", "start", "end", "Methyl_transform", "Repli")
# remove x-axis cutoff < -10
mydata.GCH_3<-mydata.GCH_3[with(mydata.GCH_3,mydata.GCH_3$Methyl_transform>-10),]
# Estimation of regression parameters
# fit linear models and carry out regression
model_3<-lm(Methyl_transform~Repli,data=mydata.GCH_3)
# Obtain predicted and residual values and draw figures
mydata.GCH_3$predicted <- predict(model_3)
mydata.GCH_3$residuals <- residuals(model_3)
suma<-summary(model_3)
# R2 represents the proportion of variance, in the outcome variable y, that may be predicted by knowing the value of the x variables. 
# An R2 value close to 1 indicates that the model explains a large portion of the variance in the outcome variable.
result_r.squared<-round(suma$r.squared,3)
# the coefficient of RT
a_2<-round(model_3$coefficients[2],3)
a_3<-round(model_3$coefficients[1],3)
n=dim(mydata.GCH_3)[1]
# combine
if(a_3<0){result3<-paste0("Methyl_transform=",a_2,"XRepli",a_3,"; R^2=",result_r.squared,"; n=",n)}else{result3<-paste0("Methyl_transform=",a_2,"XRepli+",a_3,"; R^2=",result_r.squared,"; n=",n)}

p3=ggplot()+
	geom_hex(mapping=aes(x=mydata.GCH_3$Methyl,y=mydata.GCH_3$predicted),color="black")+
	labs(x = "Observed", y = "Model-predicted", title = "chr16: The Scatter Plot of mGpC(JY608-bg)",subtitle=result3)+
	theme(plot.title = element_text(hjust = 0.5,vjust = 1),plot.subtitle = element_text(hjust = 0.5,size=6))+
	expand_limits(x = 0)
	
# Section B: For WCG
# read data from file
mydata.WCG_4<-read.table("chr16.WCG.bed",sep="\t")
# set column names
mydata.WCG_4<-cbind(mydata.WCG_4[,1:3],log(mydata.WCG_4[,4]),mydata.WCG_4[,5])
colnames(mydata.WCG_4)<-c("chrom", "start", "end", "Methyl_transform", "Repli")
# remove x-axis cutoff < -10
mydata.WCG_4<-mydata.WCG_4[with(mydata.WCG_4,mydata.WCG_4$Methyl_transform>-10),]
# Estimation of regression parameters
# fit linear models and carry out regression
model_4<-lm(Methyl_transform~Repli,data=mydata.WCG_4)
# Obtain predicted and residual values and draw figures
mydata.WCG_4$predicted <- predict(model_4)
mydata.WCG_4$residuals <- residuals(model_4)
suma<-summary(model_4)
# R2 represents the proportion of variance, in the outcome variable y, that may be predicted by knowing the value of the x variables. 
# An R2 value close to 1 indicates that the model explains a large portion of the variance in the outcome variable.
result_r.squared<-round(suma$r.squared,3)
# the coefficient of RT
a_2<-round(model_4$coefficients[2],3)
a_3<-round(model_4$coefficients[1],3)
n=dim(mydata.WCG_4)[1]
# combine
if(a_3<0){result4<-paste0("Methyl_transform=",a_2,"XRepli",a_3,"; R^2=",result_r.squared,"; n=",n)}else{result4<-paste0("Methyl_transform=",a_2,"XRepli+",a_3,"; R^2=",result_r.squared,"; n=",n)}

p4=ggplot()+
	geom_hex(mapping=aes(x=mydata.WCG_4$Methyl,y=mydata.WCG_4$predicted),color="black")+
	labs(x = "Observed", y = "Model-predicted", title = "chr16: The Scatter Plot of mCpG(JY608-bg)",subtitle=result4)+
	theme(plot.title = element_text(hjust = 0.5,vjust = 1),plot.subtitle = element_text(hjust = 0.5,size=6))+
	expand_limits(x = 0)
	
ggplot2.multiplot(p1,p2,p3,p4,cols=2)		



# chr17+18
# Section A: For GCH
# read data from file
mydata.GCH_1<-read.table("chr17.GCH.bed",sep="\t")
# set column names
mydata.GCH_1<-cbind(mydata.GCH_1[,1:3],log(mydata.GCH_1[,4]),mydata.GCH_1[,5])
colnames(mydata.GCH_1)<-c("chrom", "start", "end", "Methyl_transform", "Repli")
# remove x-axis cutoff < -10
mydata.GCH_1<-mydata.GCH_1[with(mydata.GCH_1,mydata.GCH_1$Methyl_transform>-10),]
# Estimation of regression parameters
# fit linear models and carry out regression
model_1<-lm(Methyl_transform~Repli,data=mydata.GCH_1)
# Obtain predicted and residual values and draw figures
mydata.GCH_1$predicted <- predict(model_1)
mydata.GCH_1$residuals <- residuals(model_1)
suma<-summary(model_1)
# R2 represents the proportion of variance, in the outcome variable y, that may be predicted by knowing the value of the x variables. 
# An R2 value close to 1 indicates that the model explains a large portion of the variance in the outcome variable.
result_r.squared<-round(suma$r.squared,3)
# the coefficient of RT
a_2<-round(model_1$coefficients[2],3)
a_3<-round(model_1$coefficients[1],3)
n=dim(mydata.GCH_1)[1]
# combine
if(a_3<0){result1<-paste0("Methyl_transform=",a_2,"XRepli",a_3,"; R^2=",result_r.squared,"; n=",n)}else{result1<-paste0("Methyl_transform=",a_2,"XRepli+",a_3,"; R^2=",result_r.squared,"; n=",n)}

p1=ggplot()+
	geom_hex(mapping=aes(x=mydata.GCH_1$Methyl,y=mydata.GCH_1$predicted),color="black")+
	labs(x = "Observed", y = "Model-predicted", title = "chr17: The Scatter Plot of mGpC(JY608-bg)",subtitle=result1)+
	theme(plot.title = element_text(hjust = 0.5,vjust = 1),plot.subtitle = element_text(hjust = 0.5,size=6))+
	expand_limits(x = 0)
	
# Section B: For WCG
# read data from file
mydata.WCG_2<-read.table("chr17.WCG.bed",sep="\t")
# set column names
mydata.WCG_2<-cbind(mydata.WCG_2[,1:3],log(mydata.WCG_2[,4]),mydata.WCG_2[,5])
colnames(mydata.WCG_2)<-c("chrom", "start", "end", "Methyl_transform", "Repli")
# remove x-axis cutoff < -10
mydata.WCG_2<-mydata.WCG_2[with(mydata.WCG_2,mydata.WCG_2$Methyl_transform>-10),]
# Estimation of regression parameters
# fit linear models and carry out regression
model_2<-lm(Methyl_transform~Repli,data=mydata.WCG_2)
# Obtain predicted and residual values and draw figures
mydata.WCG_2$predicted <- predict(model_2)
mydata.WCG_2$residuals <- residuals(model_2)
suma<-summary(model_2)
# R2 represents the proportion of variance, in the outcome variable y, that may be predicted by knowing the value of the x variables. 
# An R2 value close to 1 indicates that the model explains a large portion of the variance in the outcome variable.
result_r.squared<-round(suma$r.squared,3)
# the coefficient of RT
a_2<-round(model_2$coefficients[2],3)
a_3<-round(model_2$coefficients[1],3)
n=dim(mydata.WCG_2)[1]
# combine
if(a_3<0){result2<-paste0("Methyl_transform=",a_2,"XRepli",a_3,"; R^2=",result_r.squared,"; n=",n)}else{result2<-paste0("Methyl_transform=",a_2,"XRepli+",a_3,"; R^2=",result_r.squared,"; n=",n)}

p2=ggplot()+
	geom_hex(mapping=aes(x=mydata.WCG_2$Methyl,y=mydata.WCG_2$predicted),color="black")+
	labs(x = "Observed", y = "Model-predicted", title = "chr17: The Scatter Plot of mCpG(JY608-bg)",subtitle=result2)+
	theme(plot.title = element_text(hjust = 0.5,vjust = 1),plot.subtitle = element_text(hjust = 0.5,size=6))+
	expand_limits(x = 0)
	
# Section A: For GCH
# read data from file
mydata.GCH_3<-read.table("chr18.GCH.bed",sep="\t")
# set column names
mydata.GCH_3<-cbind(mydata.GCH_3[,1:3],log(mydata.GCH_3[,4]),mydata.GCH_3[,5])
colnames(mydata.GCH_3)<-c("chrom", "start", "end", "Methyl_transform", "Repli")
# remove x-axis cutoff < -10
mydata.GCH_3<-mydata.GCH_3[with(mydata.GCH_3,mydata.GCH_3$Methyl_transform>-10),]
# Estimation of regression parameters
# fit linear models and carry out regression
model_3<-lm(Methyl_transform~Repli,data=mydata.GCH_3)
# Obtain predicted and residual values and draw figures
mydata.GCH_3$predicted <- predict(model_3)
mydata.GCH_3$residuals <- residuals(model_3)
suma<-summary(model_3)
# R2 represents the proportion of variance, in the outcome variable y, that may be predicted by knowing the value of the x variables. 
# An R2 value close to 1 indicates that the model explains a large portion of the variance in the outcome variable.
result_r.squared<-round(suma$r.squared,3)
# the coefficient of RT
a_2<-round(model_3$coefficients[2],3)
a_3<-round(model_3$coefficients[1],3)
n=dim(mydata.GCH_3)[1]
# combine
if(a_3<0){result3<-paste0("Methyl_transform=",a_2,"XRepli",a_3,"; R^2=",result_r.squared,"; n=",n)}else{result3<-paste0("Methyl_transform=",a_2,"XRepli+",a_3,"; R^2=",result_r.squared,"; n=",n)}

p3=ggplot()+
	geom_hex(mapping=aes(x=mydata.GCH_3$Methyl,y=mydata.GCH_3$predicted),color="black")+
	labs(x = "Observed", y = "Model-predicted", title = "chr18: The Scatter Plot of mGpC(JY608-bg)",subtitle=result3)+
	theme(plot.title = element_text(hjust = 0.5,vjust = 1),plot.subtitle = element_text(hjust = 0.5,size=6))+
	expand_limits(x = 0)
	
# Section B: For WCG
# read data from file
mydata.WCG_4<-read.table("chr18.WCG.bed",sep="\t")
# set column names
mydata.WCG_4<-cbind(mydata.WCG_4[,1:3],log(mydata.WCG_4[,4]),mydata.WCG_4[,5])
colnames(mydata.WCG_4)<-c("chrom", "start", "end", "Methyl_transform", "Repli")
# remove x-axis cutoff < -10
mydata.WCG_4<-mydata.WCG_4[with(mydata.WCG_4,mydata.WCG_4$Methyl_transform>-10),]
# Estimation of regression parameters
# fit linear models and carry out regression
model_4<-lm(Methyl_transform~Repli,data=mydata.WCG_4)
# Obtain predicted and residual values and draw figures
mydata.WCG_4$predicted <- predict(model_4)
mydata.WCG_4$residuals <- residuals(model_4)
suma<-summary(model_4)
# R2 represents the proportion of variance, in the outcome variable y, that may be predicted by knowing the value of the x variables. 
# An R2 value close to 1 indicates that the model explains a large portion of the variance in the outcome variable.
result_r.squared<-round(suma$r.squared,3)
# the coefficient of RT
a_2<-round(model_4$coefficients[2],3)
a_3<-round(model_4$coefficients[1],3)
n=dim(mydata.WCG_4)[1]
# combine
if(a_3<0){result4<-paste0("Methyl_transform=",a_2,"XRepli",a_3,"; R^2=",result_r.squared,"; n=",n)}else{result4<-paste0("Methyl_transform=",a_2,"XRepli+",a_3,"; R^2=",result_r.squared,"; n=",n)}

p4=ggplot()+
	geom_hex(mapping=aes(x=mydata.WCG_4$Methyl,y=mydata.WCG_4$predicted),color="black")+
	labs(x = "Observed", y = "Model-predicted", title = "chr18: The Scatter Plot of mCpG(JY608-bg)",subtitle=result4)+
	theme(plot.title = element_text(hjust = 0.5,vjust = 1),plot.subtitle = element_text(hjust = 0.5,size=6))+
	expand_limits(x = 0)
	
ggplot2.multiplot(p1,p2,p3,p4,cols=2)		



# chr19+20
# Section A: For GCH
# read data from file
mydata.GCH_1<-read.table("chr19.GCH.bed",sep="\t")
# set column names
mydata.GCH_1<-cbind(mydata.GCH_1[,1:3],log(mydata.GCH_1[,4]),mydata.GCH_1[,5])
colnames(mydata.GCH_1)<-c("chrom", "start", "end", "Methyl_transform", "Repli")
# remove x-axis cutoff < -10
mydata.GCH_1<-mydata.GCH_1[with(mydata.GCH_1,mydata.GCH_1$Methyl_transform>-10),]
# Estimation of regression parameters
# fit linear models and carry out regression
model_1<-lm(Methyl_transform~Repli,data=mydata.GCH_1)
# Obtain predicted and residual values and draw figures
mydata.GCH_1$predicted <- predict(model_1)
mydata.GCH_1$residuals <- residuals(model_1)
suma<-summary(model_1)
# R2 represents the proportion of variance, in the outcome variable y, that may be predicted by knowing the value of the x variables. 
# An R2 value close to 1 indicates that the model explains a large portion of the variance in the outcome variable.
result_r.squared<-round(suma$r.squared,3)
# the coefficient of RT
a_2<-round(model_1$coefficients[2],3)
a_3<-round(model_1$coefficients[1],3)
n=dim(mydata.GCH_1)[1]
# combine
if(a_3<0){result1<-paste0("Methyl_transform=",a_2,"XRepli",a_3,"; R^2=",result_r.squared,"; n=",n)}else{result1<-paste0("Methyl_transform=",a_2,"XRepli+",a_3,"; R^2=",result_r.squared,"; n=",n)}

p1=ggplot()+
	geom_hex(mapping=aes(x=mydata.GCH_1$Methyl,y=mydata.GCH_1$predicted),color="black")+
	labs(x = "Observed", y = "Model-predicted", title = "chr19: The Scatter Plot of mGpC(JY608-bg)",subtitle=result1)+
	theme(plot.title = element_text(hjust = 0.5,vjust = 1),plot.subtitle = element_text(hjust = 0.5,size=6))+
	expand_limits(x = 0)
	
# Section B: For WCG
# read data from file
mydata.WCG_2<-read.table("chr19.WCG.bed",sep="\t")
# set column names
mydata.WCG_2<-cbind(mydata.WCG_2[,1:3],log(mydata.WCG_2[,4]),mydata.WCG_2[,5])
colnames(mydata.WCG_2)<-c("chrom", "start", "end", "Methyl_transform", "Repli")
# remove x-axis cutoff < -10
mydata.WCG_2<-mydata.WCG_2[with(mydata.WCG_2,mydata.WCG_2$Methyl_transform>-10),]
# Estimation of regression parameters
# fit linear models and carry out regression
model_2<-lm(Methyl_transform~Repli,data=mydata.WCG_2)
# Obtain predicted and residual values and draw figures
mydata.WCG_2$predicted <- predict(model_2)
mydata.WCG_2$residuals <- residuals(model_2)
suma<-summary(model_2)
# R2 represents the proportion of variance, in the outcome variable y, that may be predicted by knowing the value of the x variables. 
# An R2 value close to 1 indicates that the model explains a large portion of the variance in the outcome variable.
result_r.squared<-round(suma$r.squared,3)
# the coefficient of RT
a_2<-round(model_2$coefficients[2],3)
a_3<-round(model_2$coefficients[1],3)
n=dim(mydata.WCG_2)[1]
# combine
if(a_3<0){result2<-paste0("Methyl_transform=",a_2,"XRepli",a_3,"; R^2=",result_r.squared,"; n=",n)}else{result2<-paste0("Methyl_transform=",a_2,"XRepli+",a_3,"; R^2=",result_r.squared,"; n=",n)}

p2=ggplot()+
	geom_hex(mapping=aes(x=mydata.WCG_2$Methyl,y=mydata.WCG_2$predicted),color="black")+
	labs(x = "Observed", y = "Model-predicted", title = "chr19: The Scatter Plot of mCpG(JY608-bg)",subtitle=result2)+
	theme(plot.title = element_text(hjust = 0.5,vjust = 1),plot.subtitle = element_text(hjust = 0.5,size=6))+
	expand_limits(x = 0)
	
# Section A: For GCH
# read data from file
mydata.GCH_3<-read.table("chr20.GCH.bed",sep="\t")
# set column names
mydata.GCH_3<-cbind(mydata.GCH_3[,1:3],log(mydata.GCH_3[,4]),mydata.GCH_3[,5])
colnames(mydata.GCH_3)<-c("chrom", "start", "end", "Methyl_transform", "Repli")
# remove x-axis cutoff < -10
mydata.GCH_3<-mydata.GCH_3[with(mydata.GCH_3,mydata.GCH_3$Methyl_transform>-10),]
# Estimation of regression parameters
# fit linear models and carry out regression
model_3<-lm(Methyl_transform~Repli,data=mydata.GCH_3)
# Obtain predicted and residual values and draw figures
mydata.GCH_3$predicted <- predict(model_3)
mydata.GCH_3$residuals <- residuals(model_3)
suma<-summary(model_3)
# R2 represents the proportion of variance, in the outcome variable y, that may be predicted by knowing the value of the x variables. 
# An R2 value close to 1 indicates that the model explains a large portion of the variance in the outcome variable.
result_r.squared<-round(suma$r.squared,3)
# the coefficient of RT
a_2<-round(model_3$coefficients[2],3)
a_3<-round(model_3$coefficients[1],3)
n=dim(mydata.GCH_3)[1]
# combine
if(a_3<0){result3<-paste0("Methyl_transform=",a_2,"XRepli",a_3,"; R^2=",result_r.squared,"; n=",n)}else{result3<-paste0("Methyl_transform=",a_2,"XRepli+",a_3,"; R^2=",result_r.squared,"; n=",n)}

p3=ggplot()+
	geom_hex(mapping=aes(x=mydata.GCH_3$Methyl,y=mydata.GCH_3$predicted),color="black")+
	labs(x = "Observed", y = "Model-predicted", title = "chr20: The Scatter Plot of mGpC(JY608-bg)",subtitle=result3)+
	theme(plot.title = element_text(hjust = 0.5,vjust = 1),plot.subtitle = element_text(hjust = 0.5,size=6))+
	expand_limits(x = 0)
	
# Section B: For WCG
# read data from file
mydata.WCG_4<-read.table("chr20.WCG.bed",sep="\t")
# set column names
mydata.WCG_4<-cbind(mydata.WCG_4[,1:3],log(mydata.WCG_4[,4]),mydata.WCG_4[,5])
colnames(mydata.WCG_4)<-c("chrom", "start", "end", "Methyl_transform", "Repli")
# remove x-axis cutoff < -10
mydata.WCG_4<-mydata.WCG_4[with(mydata.WCG_4,mydata.WCG_4$Methyl_transform>-10),]
# Estimation of regression parameters
# fit linear models and carry out regression
model_4<-lm(Methyl_transform~Repli,data=mydata.WCG_4)
# Obtain predicted and residual values and draw figures
mydata.WCG_4$predicted <- predict(model_4)
mydata.WCG_4$residuals <- residuals(model_4)
suma<-summary(model_4)
# R2 represents the proportion of variance, in the outcome variable y, that may be predicted by knowing the value of the x variables. 
# An R2 value close to 1 indicates that the model explains a large portion of the variance in the outcome variable.
result_r.squared<-round(suma$r.squared,3)
# the coefficient of RT
a_2<-round(model_4$coefficients[2],3)
a_3<-round(model_4$coefficients[1],3)
n=dim(mydata.WCG_4)[1]
# combine
if(a_3<0){result4<-paste0("Methyl_transform=",a_2,"XRepli",a_3,"; R^2=",result_r.squared,"; n=",n)}else{result4<-paste0("Methyl_transform=",a_2,"XRepli+",a_3,"; R^2=",result_r.squared,"; n=",n)}

p4=ggplot()+
	geom_hex(mapping=aes(x=mydata.WCG_4$Methyl,y=mydata.WCG_4$predicted),color="black")+
	labs(x = "Observed", y = "Model-predicted", title = "chr20: The Scatter Plot of mCpG(JY608-bg)",subtitle=result4)+
	theme(plot.title = element_text(hjust = 0.5,vjust = 1),plot.subtitle = element_text(hjust = 0.5,size=6))+
	expand_limits(x = 0)
	
ggplot2.multiplot(p1,p2,p3,p4,cols=2)		



# chr21+x
# Section A: For GCH
# read data from file
mydata.GCH_1<-read.table("chr21.GCH.bed",sep="\t")
# set column names
mydata.GCH_1<-cbind(mydata.GCH_1[,1:3],log(mydata.GCH_1[,4]),mydata.GCH_1[,5])
colnames(mydata.GCH_1)<-c("chrom", "start", "end", "Methyl_transform", "Repli")
# remove x-axis cutoff < -10
mydata.GCH_1<-mydata.GCH_1[with(mydata.GCH_1,mydata.GCH_1$Methyl_transform>-10),]
# Estimation of regression parameters
# fit linear models and carry out regression
model_1<-lm(Methyl_transform~Repli,data=mydata.GCH_1)
# Obtain predicted and residual values and draw figures
mydata.GCH_1$predicted <- predict(model_1)
mydata.GCH_1$residuals <- residuals(model_1)
suma<-summary(model_1)
# R2 represents the proportion of variance, in the outcome variable y, that may be predicted by knowing the value of the x variables. 
# An R2 value close to 1 indicates that the model explains a large portion of the variance in the outcome variable.
result_r.squared<-round(suma$r.squared,3)
# the coefficient of RT
a_2<-round(model_1$coefficients[2],3)
a_3<-round(model_1$coefficients[1],3)
n=dim(mydata.GCH_1)[1]
# combine
if(a_3<0){result1<-paste0("Methyl_transform=",a_2,"XRepli",a_3,"; R^2=",result_r.squared,"; n=",n)}else{result1<-paste0("Methyl_transform=",a_2,"XRepli+",a_3,"; R^2=",result_r.squared,"; n=",n)}

p1=ggplot()+
	geom_hex(mapping=aes(x=mydata.GCH_1$Methyl,y=mydata.GCH_1$predicted),color="black")+
	labs(x = "Observed", y = "Model-predicted", title = "chr21: The Scatter Plot of mGpC(JY608-bg)",subtitle=result1)+
	theme(plot.title = element_text(hjust = 0.5,vjust = 1),plot.subtitle = element_text(hjust = 0.5,size=6))+
	expand_limits(x = 0)
	
# Section B: For WCG
# read data from file
mydata.WCG_2<-read.table("chr21.WCG.bed",sep="\t")
# set column names
mydata.WCG_2<-cbind(mydata.WCG_2[,1:3],log(mydata.WCG_2[,4]),mydata.WCG_2[,5])
colnames(mydata.WCG_2)<-c("chrom", "start", "end", "Methyl_transform", "Repli")
# remove x-axis cutoff < -10
mydata.WCG_2<-mydata.WCG_2[with(mydata.WCG_2,mydata.WCG_2$Methyl_transform>-10),]
# Estimation of regression parameters
# fit linear models and carry out regression
model_2<-lm(Methyl_transform~Repli,data=mydata.WCG_2)
# Obtain predicted and residual values and draw figures
mydata.WCG_2$predicted <- predict(model_2)
mydata.WCG_2$residuals <- residuals(model_2)
suma<-summary(model_2)
# R2 represents the proportion of variance, in the outcome variable y, that may be predicted by knowing the value of the x variables. 
# An R2 value close to 1 indicates that the model explains a large portion of the variance in the outcome variable.
result_r.squared<-round(suma$r.squared,3)
# the coefficient of RT
a_2<-round(model_2$coefficients[2],3)
a_3<-round(model_2$coefficients[1],3)
n=dim(mydata.WCG_2)[1]
# combine
if(a_3<0){result2<-paste0("Methyl_transform=",a_2,"XRepli",a_3,"; R^2=",result_r.squared,"; n=",n)}else{result2<-paste0("Methyl_transform=",a_2,"XRepli+",a_3,"; R^2=",result_r.squared,"; n=",n)}

p2=ggplot()+
	geom_hex(mapping=aes(x=mydata.WCG_2$Methyl,y=mydata.WCG_2$predicted),color="black")+
	labs(x = "Observed", y = "Model-predicted", title = "chr21: The Scatter Plot of mCpG(JY608-bg)",subtitle=result2)+
	theme(plot.title = element_text(hjust = 0.5,vjust = 1),plot.subtitle = element_text(hjust = 0.5,size=6))+
	expand_limits(x = 0)
	
# Section A: For GCH
# read data from file
mydata.GCH_3<-read.table("chrX.GCH.bed",sep="\t")
# set column names
mydata.GCH_3<-cbind(mydata.GCH_3[,1:3],log(mydata.GCH_3[,4]),mydata.GCH_3[,5])
colnames(mydata.GCH_3)<-c("chrom", "start", "end", "Methyl_transform", "Repli")
# remove x-axis cutoff < -10
mydata.GCH_3<-mydata.GCH_3[with(mydata.GCH_3,mydata.GCH_3$Methyl_transform>-10),]
# Estimation of regression parameters
# fit linear models and carry out regression
model_3<-lm(Methyl_transform~Repli,data=mydata.GCH_3)
# Obtain predicted and residual values and draw figures
mydata.GCH_3$predicted <- predict(model_3)
mydata.GCH_3$residuals <- residuals(model_3)
suma<-summary(model_3)
# R2 represents the proportion of variance, in the outcome variable y, that may be predicted by knowing the value of the x variables. 
# An R2 value close to 1 indicates that the model explains a large portion of the variance in the outcome variable.
result_r.squared<-round(suma$r.squared,3)
# the coefficient of RT
a_2<-round(model_3$coefficients[2],3)
a_3<-round(model_3$coefficients[1],3)
n=dim(mydata.GCH_3)[1]
# combine
if(a_3<0){result3<-paste0("Methyl_transform=",a_2,"XRepli",a_3,"; R^2=",result_r.squared,"; n=",n)}else{result3<-paste0("Methyl_transform=",a_2,"XRepli+",a_3,"; R^2=",result_r.squared,"; n=",n)}

p3=ggplot()+
	geom_hex(mapping=aes(x=mydata.GCH_3$Methyl,y=mydata.GCH_3$predicted),color="black")+
	labs(x = "Observed", y = "Model-predicted", title = "chrX: The Scatter Plot of mGpC(JY608-bg)",subtitle=result3)+
	theme(plot.title = element_text(hjust = 0.5,vjust = 1),plot.subtitle = element_text(hjust = 0.5,size=6))+
	expand_limits(x = 0)
	
# Section B: For WCG
# read data from file
mydata.WCG_4<-read.table("chrX.WCG.bed",sep="\t")
# set column names
mydata.WCG_4<-cbind(mydata.WCG_4[,1:3],log(mydata.WCG_4[,4]),mydata.WCG_4[,5])
colnames(mydata.WCG_4)<-c("chrom", "start", "end", "Methyl_transform", "Repli")
# remove x-axis cutoff < -10
mydata.WCG_4<-mydata.WCG_4[with(mydata.WCG_4,mydata.WCG_4$Methyl_transform>-10),]
# Estimation of regression parameters
# fit linear models and carry out regression
model_4<-lm(Methyl_transform~Repli,data=mydata.WCG_4)
# Obtain predicted and residual values and draw figures
mydata.WCG_4$predicted <- predict(model_4)
mydata.WCG_4$residuals <- residuals(model_4)
suma<-summary(model_4)
# R2 represents the proportion of variance, in the outcome variable y, that may be predicted by knowing the value of the x variables. 
# An R2 value close to 1 indicates that the model explains a large portion of the variance in the outcome variable.
result_r.squared<-round(suma$r.squared,3)
# the coefficient of RT
a_2<-round(model_4$coefficients[2],3)
a_3<-round(model_4$coefficients[1],3)
n=dim(mydata.WCG_4)[1]
# combine
if(a_3<0){result4<-paste0("Methyl_transform=",a_2,"XRepli",a_3,"; R^2=",result_r.squared,"; n=",n)}else{result4<-paste0("Methyl_transform=",a_2,"XRepli+",a_3,"; R^2=",result_r.squared,"; n=",n)}

p4=ggplot()+
	geom_hex(mapping=aes(x=mydata.WCG_4$Methyl,y=mydata.WCG_4$predicted),color="black")+
	labs(x = "Observed", y = "Model-predicted", title = "chrX: The Scatter Plot of mCpG(JY608-bg)",subtitle=result4)+
	theme(plot.title = element_text(hjust = 0.5,vjust = 1),plot.subtitle = element_text(hjust = 0.5,size=6))+
	expand_limits(x = 0)
	
ggplot2.multiplot(p1,p2,p3,p4,cols=2)		

# chry+x
# Section A: For GCH
# read data from file
mydata.GCH_1<-read.table("chrY.GCH.bed",sep="\t")
# set column names
mydata.GCH_1<-cbind(mydata.GCH_1[,1:3],log(mydata.GCH_1[,4]),mydata.GCH_1[,5])
colnames(mydata.GCH_1)<-c("chrom", "start", "end", "Methyl_transform", "Repli")
# remove x-axis cutoff < -10
mydata.GCH_1<-mydata.GCH_1[with(mydata.GCH_1,mydata.GCH_1$Methyl_transform>-10),]
# Estimation of regression parameters
# fit linear models and carry out regression
model_1<-lm(Methyl_transform~Repli,data=mydata.GCH_1)
# Obtain predicted and residual values and draw figures
mydata.GCH_1$predicted <- predict(model_1)
mydata.GCH_1$residuals <- residuals(model_1)
suma<-summary(model_1)
# R2 represents the proportion of variance, in the outcome variable y, that may be predicted by knowing the value of the x variables. 
# An R2 value close to 1 indicates that the model explains a large portion of the variance in the outcome variable.
result_r.squared<-round(suma$r.squared,3)
# the coefficient of RT
a_2<-round(model_1$coefficients[2],3)
a_3<-round(model_1$coefficients[1],3)
n=dim(mydata.GCH_1)[1]
# combine
if(a_3<0){result1<-paste0("Methyl_transform=",a_2,"XRepli",a_3,"; R^2=",result_r.squared,"; n=",n)}else{result1<-paste0("Methyl_transform=",a_2,"XRepli+",a_3,"; R^2=",result_r.squared,"; n=",n)}

p1=ggplot()+
	geom_hex(mapping=aes(x=mydata.GCH_1$Methyl,y=mydata.GCH_1$predicted),color="black")+
	labs(x = "Observed", y = "Model-predicted", title = "chrY: The Scatter Plot of mGpC(JY608-bg)",subtitle=result1)+
	theme(plot.title = element_text(hjust = 0.5,vjust = 1),plot.subtitle = element_text(hjust = 0.5,size=6))+
	expand_limits(x = 0)
	
# Section B: For WCG
# read data from file
mydata.WCG_2<-read.table("chrY.WCG.bed",sep="\t")
# set column names
mydata.WCG_2<-cbind(mydata.WCG_2[,1:3],log(mydata.WCG_2[,4]),mydata.WCG_2[,5])
colnames(mydata.WCG_2)<-c("chrom", "start", "end", "Methyl_transform", "Repli")
# remove x-axis cutoff < -10
mydata.WCG_2<-mydata.WCG_2[with(mydata.WCG_2,mydata.WCG_2$Methyl_transform>-10),]
# Estimation of regression parameters
# fit linear models and carry out regression
model_2<-lm(Methyl_transform~Repli,data=mydata.WCG_2)
# Obtain predicted and residual values and draw figures
mydata.WCG_2$predicted <- predict(model_2)
mydata.WCG_2$residuals <- residuals(model_2)
suma<-summary(model_2)
# R2 represents the proportion of variance, in the outcome variable y, that may be predicted by knowing the value of the x variables. 
# An R2 value close to 1 indicates that the model explains a large portion of the variance in the outcome variable.
result_r.squared<-round(suma$r.squared,3)
# the coefficient of RT
a_2<-round(model_2$coefficients[2],3)
a_3<-round(model_2$coefficients[1],3)
n=dim(mydata.WCG_2)[1]
# combine
if(a_3<0){result2<-paste0("Methyl_transform=",a_2,"XRepli",a_3,"; R^2=",result_r.squared,"; n=",n)}else{result2<-paste0("Methyl_transform=",a_2,"XRepli+",a_3,"; R^2=",result_r.squared,"; n=",n)}

p2=ggplot()+
	geom_hex(mapping=aes(x=mydata.WCG_2$Methyl,y=mydata.WCG_2$predicted),color="black")+
	labs(x = "Observed", y = "Model-predicted", title = "chrY: The Scatter Plot of mCpG(JY608-bg)",subtitle=result2)+
	theme(plot.title = element_text(hjust = 0.5,vjust = 1),plot.subtitle = element_text(hjust = 0.5,size=6))+
	expand_limits(x = 0)
	
# Section A: For GCH
# read data from file
mydata.GCH_3<-read.table("chrX.GCH.bed",sep="\t")
# set column names
mydata.GCH_3<-cbind(mydata.GCH_3[,1:3],log(mydata.GCH_3[,4]),mydata.GCH_3[,5])
colnames(mydata.GCH_3)<-c("chrom", "start", "end", "Methyl_transform", "Repli")
# remove x-axis cutoff < -10
mydata.GCH_3<-mydata.GCH_3[with(mydata.GCH_3,mydata.GCH_3$Methyl_transform>-10),]
# Estimation of regression parameters
# fit linear models and carry out regression
model_3<-lm(Methyl_transform~Repli,data=mydata.GCH_3)
# Obtain predicted and residual values and draw figures
mydata.GCH_3$predicted <- predict(model_3)
mydata.GCH_3$residuals <- residuals(model_3)
suma<-summary(model_3)
# R2 represents the proportion of variance, in the outcome variable y, that may be predicted by knowing the value of the x variables. 
# An R2 value close to 1 indicates that the model explains a large portion of the variance in the outcome variable.
result_r.squared<-round(suma$r.squared,3)
# the coefficient of RT
a_2<-round(model_3$coefficients[2],3)
a_3<-round(model_3$coefficients[1],3)
n=dim(mydata.GCH_3)[1]
# combine
if(a_3<0){result3<-paste0("Methyl_transform=",a_2,"XRepli",a_3,"; R^2=",result_r.squared,"; n=",n)}else{result3<-paste0("Methyl_transform=",a_2,"XRepli+",a_3,"; R^2=",result_r.squared,"; n=",n)}

p3=ggplot()+
	geom_hex(mapping=aes(x=mydata.GCH_3$Methyl,y=mydata.GCH_3$predicted),color="black")+
	labs(x = "Observed", y = "Model-predicted", title = "chrX: The Scatter Plot of mGpC(JY608-bg)",subtitle=result3)+
	theme(plot.title = element_text(hjust = 0.5,vjust = 1),plot.subtitle = element_text(hjust = 0.5,size=6))+
	expand_limits(x = 0)
	
# Section B: For WCG
# read data from file
mydata.WCG_4<-read.table("chrX.WCG.bed",sep="\t")
# set column names
mydata.WCG_4<-cbind(mydata.WCG_4[,1:3],log(mydata.WCG_4[,4]),mydata.WCG_4[,5])
colnames(mydata.WCG_4)<-c("chrom", "start", "end", "Methyl_transform", "Repli")
# remove x-axis cutoff < -10
mydata.WCG_4<-mydata.WCG_4[with(mydata.WCG_4,mydata.WCG_4$Methyl_transform>-10),]
# Estimation of regression parameters
# fit linear models and carry out regression
model_4<-lm(Methyl_transform~Repli,data=mydata.WCG_4)
# Obtain predicted and residual values and draw figures
mydata.WCG_4$predicted <- predict(model_4)
mydata.WCG_4$residuals <- residuals(model_4)
suma<-summary(model_4)
# R2 represents the proportion of variance, in the outcome variable y, that may be predicted by knowing the value of the x variables. 
# An R2 value close to 1 indicates that the model explains a large portion of the variance in the outcome variable.
result_r.squared<-round(suma$r.squared,3)
# the coefficient of RT
a_2<-round(model_4$coefficients[2],3)
a_3<-round(model_4$coefficients[1],3)
n=dim(mydata.WCG_4)[1]
# combine
if(a_3<0){result4<-paste0("Methyl_transform=",a_2,"XRepli",a_3,"; R^2=",result_r.squared,"; n=",n)}else{result4<-paste0("Methyl_transform=",a_2,"XRepli+",a_3,"; R^2=",result_r.squared,"; n=",n)}

p4=ggplot()+
	geom_hex(mapping=aes(x=mydata.WCG_4$Methyl,y=mydata.WCG_4$predicted),color="black")+
	labs(x = "Observed", y = "Model-predicted", title = "chrX: The Scatter Plot of mCpG(JY608-bg)",subtitle=result4)+
	theme(plot.title = element_text(hjust = 0.5,vjust = 1),plot.subtitle = element_text(hjust = 0.5,size=6))+
	expand_limits(x = 0)
	
ggplot2.multiplot(p1,p2,p3,p4,cols=2)		





dev.off()

#remove(additionalx)
