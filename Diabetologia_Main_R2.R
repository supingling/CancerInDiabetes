###Oct 2022 Diabetologia R2 == SL
###trends in all-cause and cancer mortality among diabetes 
#setwd("/rfs/LRWE_Proj46/Shared/SL/analysis/Mortality/R2")
library(Epi)
library(survival)
library(haven)
library(foreign)
library(broom)
library(popEpi)
library(dplyr)


rm(list=ls())
out_overall <-list()
out_sex     <-list()
out_sex1    <-list()
out_eth     <-list()
out_eth1    <-list()
out_imd     <-list()
out_imd1    <-list()
out_smk     <-list()
out_smk1    <-list()
out_bmi     <-list()
out_bmi1    <-list()
fittest     <-list()

age    <-seq(55,85,10)
period <-seq(1998.5,2018.5,1)
dur    <-c(8.4)
nd1    <-cbind(expand.grid(age, period, dur), lex.dur=1000)
colnames(nd1)<-c("age","period", "dur", "lex.dur")

age    <-c(72)
period <-seq(1998.5,2018.5,1)
dur    <-c(8.4)
nd2    <-cbind(expand.grid(age, period, dur), lex.dur=1000)
colnames(nd2)<-c("age","period", "dur", "lex.dur")

dmdata = read_dta("dmdb.dta")
attach(dmdata)

#####main analysis updated in Oct 2022
for (i in c(16:25)) {
  
  dml1 <-Lexis(entry = list(period = yearin,
                            age    = agein,
                            dur    = 0),
               exit = list(period  = out),
               exit.status = dmdata[[i]],
               id  = patid,
               data = dmdata)
  dms1  <-splitMulti(dml1, period= seq(1998,2020,1),
                     age = seq(35,111,1))
  
  a.kn <- with(subset(dms1, lex.Xst==1), quantile(age+lex.dur,(1:5-0.5)/5))
  p.kn <- with(subset(dms1, lex.Xst==1), quantile(period+lex.dur,(1:5-0.5)/5))
  d.kn <- with(subset(dms1, lex.Xst==1), quantile(dur+lex.dur,(1:5-0.5)/5))
  
  dmc1  <- glm((lex.Xst==1) ~ Ns(period, knots = p.kn) * Ns(age, knots = a.kn) + Ns(dur, knots = d.kn),
               family = poisson,
               offset = log(lex.dur),
               data   =dms1)
  fittest[[i]] <- as.data.frame(cbind(out = i, 
                                      deviance = dmc1$deviance, df = dmc1$df.residual, AIC = dmc1$aic, 
                                      p = pchisq(dmc1$deviance, dmc1$df.residual, lower.tail = F)))  
  
  p1     <-as.data.frame(ci.pred(dmc1, newdata = nd1, Exp = T))
  colnames(p1)<-c("rate","lb","ub")
  out_overall[[i]]    <-cbind(nd1, p1, out = i, gender = NA, ethnic = NA, imd = NA, smk_status = NA, bmi_grp = NA )
  
  
  if (i<22) {  
    for (j in c(1:2)) {
      a.kn <- with(subset(dms1, lex.Xst==1 & gender == j), quantile(age+lex.dur,(1:5-0.5)/5))
      p.kn <- with(subset(dms1, lex.Xst==1 & gender == j), quantile(period+lex.dur,(1:5-0.5)/5))
      
      dmc1  <- glm((lex.Xst==1) ~ Ns(period, knots = p.kn) * Ns(age, knots = a.kn) + Ns(dur, knots = d.kn),
                   family = poisson,
                   offset = log(lex.dur),
                   data   = subset(dms1,gender==j))
      
      p1     <-as.data.frame(ci.pred(dmc1, newdata = nd2, Exp = T))
      colnames(p1)<-c("rate","lb","ub")
      
      out_sex1[[j]] <-cbind(nd2, p1, out = i, gender = j, ethnic = NA, imd = NA, smk_status = NA, bmi_grp = NA)
    }
    out_sex[[i]] <-data.frame(Reduce(rbind,out_sex1))
  }
  
  if (i!=22 & i!=24) {
    for (j in c(1:2)) {
      a.kn <- with(subset(dms1, lex.Xst==1 & ethnic == j), quantile(age+lex.dur,(1:5-0.5)/5))
      p.kn <- with(subset(dms1, lex.Xst==1 & ethnic == j), quantile(period+lex.dur,(1:5-0.5)/5))
      dmc1  <- glm((lex.Xst==1) ~ Ns(period, knots = p.kn) * Ns(age, knots = a.kn) + Ns(dur, knots = d.kn),
                   family = poisson,
                   offset = log(lex.dur),
                   data   = subset(dms1,ethnic==j))
      
      p1     <-as.data.frame(ci.pred(dmc1, newdata = nd2, Exp = T))
      colnames(p1)<-c("rate","lb","ub")
      out_eth1[[j]]    <-cbind(nd2, p1, out = i, gender = NA, ethnic = j, imd = NA, smk_status = NA, bmi_grp = NA)
    }
    out_eth[[i]] <-data.frame(Reduce(rbind,out_eth1))
    
    
    for (j in c(1:5)) {
      a.kn <- with(subset(dms1, lex.Xst==1 & imd == j), quantile(age+lex.dur,(1:5-0.5)/5))
      p.kn <- with(subset(dms1, lex.Xst==1 & imd == j), quantile(period+lex.dur,(1:5-0.5)/5))
      dmc1  <- glm((lex.Xst==1) ~ Ns(period, knots = p.kn) * Ns(age, knots = a.kn) + Ns(dur, knots = d.kn),
                   family = poisson,
                   offset = log(lex.dur),
                   data   = subset(dms1,imd==j))
      
      p1     <-as.data.frame(ci.pred(dmc1, newdata = nd2, Exp = T))
      colnames(p1)<-c("rate","lb","ub")
      out_imd1[[j]]    <-cbind(nd2, p1,out=i, gender = NA, ethnic = NA, imd = j, smk_status = NA, bmi_grp = NA)
    }
    out_imd[[i]] <-data.frame(Reduce(rbind,out_imd1))
    
    
    for (j in c(1:3)) {
      a.kn <- with(subset(dms1, lex.Xst==1 & smk_status == j), quantile(age+lex.dur,(1:5-0.5)/5))
      p.kn <- with(subset(dms1, lex.Xst==1 & smk_status == j), quantile(period+lex.dur,(1:5-0.5)/5))
      dmc1  <- glm((lex.Xst==1) ~ Ns(period, knots = p.kn) * Ns(age, knots = a.kn) + Ns(dur, knots = d.kn),
                   family = poisson,
                   offset = log(lex.dur),
                   data   = subset(dms1,smk_status==j))
      
      p1     <-as.data.frame(ci.pred(dmc1, newdata = nd2, Exp = T))
      colnames(p1)<-c("rate","lb","ub")
      out_smk1[[j]]    <-cbind(nd2, p1,out=i, gender = NA, ethnic = NA, imd = NA, smk_status = j, bmi_grp = NA)
    }
    out_smk[[i]] <-data.frame(Reduce(rbind,out_smk1)) 
    
    
    for (j in c(2:5)) {
      a.kn <- with(subset(dms1, lex.Xst==1 & bmi_grp == j), quantile(age+lex.dur,(1:5-0.5)/5))
      p.kn <- with(subset(dms1, lex.Xst==1 & bmi_grp == j), quantile(period+lex.dur,(1:5-0.5)/5))
      dmc1  <- glm((lex.Xst==1) ~ Ns(period, knots = p.kn) * Ns(age, knots = a.kn) + Ns(dur, knots = d.kn),
                   family = poisson,
                   offset = log(lex.dur),
                   data   = subset(dms1,bmi_grp==j))
      
      p1     <-as.data.frame(ci.pred(dmc1, newdata = nd2, Exp = T))
      colnames(p1)<-c("rate","lb","ub")
      out_bmi1[[j]]    <-cbind(nd2, p1,out=i, gender = NA, ethnic = NA, imd = NA, smk_status = NA, bmi_grp = j)
    }
    out_bmi[[i]] <-data.frame(Reduce(rbind,out_bmi1))
  }
}
fitest <- data.frame(Reduce(rbind, fittest))
write.dta(as.data.frame(fitest), file = paste0("DM_fit_R2.dta"))

out_all <- rbind(do.call(rbind.data.frame, out_overall),do.call(rbind.data.frame, out_sex), 
                 do.call(rbind.data.frame, out_eth), do.call(rbind.data.frame, out_imd),
                 do.call(rbind.data.frame, out_smk), do.call(rbind.data.frame, out_bmi))
write.dta(as.data.frame(out_all), file = paste0("DM_res_R2.dta"))
