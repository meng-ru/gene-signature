## this code analyse Celebi's biomarker data for melanoma staging 
library(survival)
library(rms)
library(survplot)
library(CPE)

melanoma <- read.table("data/melanoma_clin.csv", header=TRUE,sep=",")
#remove epi2 and those without survival data
mela1<-melanoma[melanoma$subtype!="Epi2",] # 6 removed left 45
mela2<-mela1[mela1$os!="n.a.",]# 5 removed left 40
mela2$death<-ifelse(mela2$status=="Deceased",1,0)
mela2$os <- as.numeric(as.character(mela2$os))
mela2$stage[mela2$ajcc=="Ia"|mela2$ajcc=="Ib"]<-1 
mela2$stage[mela2$ajcc=="IIa"|mela2$ajcc=="IIb"|mela2$ajcc=="IIc"]<-2
mela2$stage[mela2$ajcc=="IIIa"|mela2$ajcc=="IIIb"|mela2$ajcc=="IIIc"]<-3 
mela2$tcat[mela2$t=="T1"|mela2$t=="T1a"]<-1
mela2$tcat[mela2$t=="T2"|mela2$t=="T2a"|mela2$t=="T2b"]<-2
mela2$tcat[mela2$t=="T3"|mela2$t=="T3a"|mela2$t=="T3b"]<-3
mela2$tcat[mela2$t=="T4a"|mela2$t=="T4b"]<-4

mela2$rcat[mela2$subtype=="Epi1"]<-0
mela2$rcat[mela2$subtype=="Epi3"]<-1

mela2$ulcercat[mela2$ulcer=="Yes"]<-1
mela2$ulcercat[mela2$ulcer=="No"]<-0






table(mela2$death)
summary(mela2$os[mela2$death==1])

km<- survfit(Surv(os, death) ~ 1, data = mela2, conf.type = "log-log")
plot(km, conf.int=FALSE, mark.time=TRUE,
     mark=3, col=1, lty=1, lwd=1, cex=1, 
     log=FALSE, xscale=1, yscale=1, 
     firstx=0, firsty=1, ymin=0,
     main="Overall survival",
     xlab="Months from diagnosis", ylab="Survial probability", xaxs="S")
survdiff(Surv(os, death) ~ subtype, data=mela2)
legend(100, 0.6, c("15/40 died"), bty="n")
legend(100,0.5,c("Median survival 219 months, 95% CI 104-NA"),bty="n")
legend(100,0.4,c("Log rank test p value <0.0001"),bty="n")


km.by.type <- survfit(Surv(os, death == 1) ~ subtype, data = mela2, conf.type = "log-log")
plot(km.by.type, conf.int=FALSE, mark.time=TRUE,
     mark=3, col=c(1,2), lty=1, lwd=1, cex=1, 
     log=FALSE, xscale=1, yscale=1, 
     firstx=0, firsty=1, ymin=0,
     main="Overall survival stratified by genetic risk group",
     xlab="Months from diagnosis", ylab="Survial probability", xaxs="S")
survdiff(Surv(os, death) ~ subtype, data=mela2)
legend(100, 0.7, c("Low risk 4/26 died", "High risk 11/14 died"), 
       col = c(1,2),text.col = "black", lty = c(1,1),bty="n")
legend(110,0.5,c("Log rank test p value <0.0001"),bty="n")

km.by.stage <- survfit(Surv(os, death == 1) ~ stage, data = mela2, conf.type = "log-log")
plot(km.by.stage, conf.int=FALSE, mark.time=TRUE,
     mark=3, col=c(1,2,3), lty=1, lwd=1, cex=1, 
     log=FALSE, xscale=1, yscale=1, 
     firstx=0, firsty=1, ymin=0,
     main="Overall survival stratified by AJCC stage",
     xlab="Months from diagnosis", ylab="Survial probability", xaxs="S")
survdiff(Surv(os, death) ~ stage, data=mela2)
legend(150, 0.6, c("AJCC stage I 0/12 died", "AJCC stage II 8/20 died",
                   "AJCC stage III 7/8 died"), 
       col = c(1,2,3),text.col = "black", lty = c(1,1),bty="n")
legend(110,0.35,c("Log rank test p value =0.0001"),bty="n")

km.by.t <- survfit(Surv(os, death == 1) ~ tcat, data = mela2, conf.type = "log-log")
plot(km.by.t, conf.int=FALSE, mark.time=TRUE,
     mark=3, col=c(1,2,3,4), lty=1, lwd=1, cex=1, 
     log=FALSE, xscale=1, yscale=1, 
     firstx=0, firsty=1, ymin=0,main="Overall Survival stratified by T staging",
     xlab="Months from diagnosis", ylab="Survial probability", xaxs="S")
survdiff(Surv(os, death) ~ tcat, data=mela2)
legend(150, 0.65, c("T1 0/9 died", "T2 1/7 died","T3 5/13 died","T4 9/11 died"), 
       col = c(1,2,3,4),
       text.col = "black", lty =1,bty="n")
legend(150,0.35,c("Log rank test p value =0.0008"),bty="n")

### KM plots by tcat and genetic risk groups
km.by.tr <- survfit(Surv(os, death == 1) ~ tcat+subtype, data = mela2, conf.type = "log-log")
plot(km.by.tr, conf.int=FALSE, mark.time=TRUE,
     mark=3, col=c(1,2,3,4,5,6,7), lty=1, lwd=1, cex=1, 
     log=FALSE, xscale=1, yscale=1, 
     firstx=0, firsty=1, ymin=0,main="Overall Survival stratified by T staging and genetic risk groups",
        xlab="Months from diagnosis", ylab="Survial probability", xaxs="S")
survdiff(Surv(os, death) ~tcat+subtype, data=mela2)
legend(150, 0.85, c("T1 low risk 0/9 died", "T2 low risk 1/6 died",
                    "T2 high risk 0/1 died", "T3 low risk 1/9 died",
                    "T3 high risk 4/4 died", "T4 low risk 2/2 died",
                    "T4 high risk 7/9 died"), 
       col = c(1,2,3,4,5,6,7),
       text.col = "black", lty =1,bty="n")
legend(150,0.35,c("Log rank test p value <0.0001"),bty="n")


### KM plots by AJCC stage and genetic risk groups
km.by.sr <- survfit(Surv(os, death == 1) ~ stage+subtype, data = mela2, conf.type = "log-log")
plot(km.by.sr, conf.int=FALSE, mark.time=TRUE,
     mark=3, col=c(1,2,3,4,5,6), lty=1, lwd=1, cex=1, 
     log=FALSE, xscale=1, yscale=1, 
     firstx=0, firsty=1, ymin=0,main="Overall Survival stratified by AJCC staging and genetic risk groups",
     xlab="Months from diagnosis", ylab="Survial probability", xaxs="S")
survdiff(Surv(os, death) ~stage+subtype, data=mela2)
legend(150, 0.85, c("AJCC I low risk 0/11 died", "AJCC I high risk 0/1 died",
                    "AJCC II low risk 3/13 died", "AJCC II high risk 5/7 died",
                    "AJCC III low risk 1/2 died", "AJCC III high risk 6/6 died"), 
       col = c(1,2,3,4,5,6),cex=0.9,
       text.col = "black", lty =1,bty="n")
legend(150,0.45,c("Log rank test p value <0.0001"),bty="n",cex=0.9)

#univariate analysis: KM for categorical variables and Coxph for continuous variables
survdiff(Surv(os, death) ~stage, data=mela2)
survdiff(Surv(os, death) ~tcat, data=mela2)
survdiff(Surv(os, death) ~sex, data=mela2) # not significant
survdiff(Surv(os, death) ~ulcercat, data=mela2)
coxph(Surv(os,death)~age,data=mela2)


rs2<-coxph(Surv(os, death) ~ rcat + age+ as.factor(stage), data=mela2) 
phcpe(rs2, CPE.SE=TRUE,out.ties=FALSE)
rs2_nor<-coxph(Surv(os, death) ~ age+ stage, data=mela2)
phcpe(rs2_nor, CPE.SE=TRUE,out.ties=FALSE)

rt2<-coxph(Surv(os, death) ~ rcat + age+tcat, data=mela2)
phcpe(rt2, CPE.SE=TRUE,out.ties=FALSE)
rt2_nor<-coxph(Surv(os, death) ~ age+ tcat, data=mela2)
phcpe(rt2_nor, CPE.SE=TRUE,out.ties=FALSE)

mela2$subtype1<-ifelse(mela2$subtype=="Epi1",1,0)
coxph(Surv(os,death)~subtype1+stage,data=mela2)

int<-coxph(Surv(os,death)~factor(rcat,levels=c(0,1))*factor(stage,levels=
                                                  c(1,2,3)),data=mela2)

## create the variables based on both tumor thickness and genetic risk group
# create new variable rs which takes value 1 for ajcc stage 1, =2 for ajcc stage II and 
# genetic risk group low,=3 for ajcc II and high risk, =4 for ajcc III

mela2$rs[mela2$stage==1]<-1
mela2$rs[mela2$stage==2 & mela2$subtype=="Epi1"]<-2
mela2$rs[mela2$stage==2 & mela2$subtype=="Epi3"]<-3
mela2$rs[mela2$stage==3]<-4

mela2$rs<-factor(mela2$rs,levels=c(2,1,3,4))
rsint<-coxph(Surv(os,death)~rs,data=mela2)
phcpe(rsint, CPE.SE=TRUE,out.ties=FALSE)
rsint

s<-coxph(Surv(os,death)~stage,data=mela2)
phcpe(s,CPE.SE=TRUE,out.ties=FALSE)


mela2$tr[mela2$tcat==1]<-1
mela2$tr[mela2$tcat==2]<-2
mela2$tr[mela2$tcat==3 & mela2$subtype=="Epi1"]<-3
mela2$tr[mela2$tcat==3 & mela2$subtype=="Epi3"]<-4
mela2$tr[mela2$tcat==4]<-5


mela2$tr<-factor(mela2$tr,levels=c(3,1,2,4,5))
trint<-coxph(Surv(os,death)~tr,data=mela2)
phcpe(trint,CPE.SE=TRUE,out.ties=FALSE)
trint

t<-coxph(Surv(os,death)~tcat,data=mela2)
phcpe(t,CPE.SE=TRUE,out.ties=FALSE)

km.by.tr <- survfit(Surv(os, death == 1) ~ tr, data = mela2, conf.type = "log-log")
