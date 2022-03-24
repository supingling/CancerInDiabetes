###trends in all-cause and cancer mortality among diabetes non-sex-specific cancers

library(Epi)
library(survival)
library(haven)
library(foreign)
library(broom)
library(popEpi)
library(dplyr)

rm(list=ls())
dmdata <-read_dta("dmdb.dta")
nrow(dmdata)
attach(dmdata)

outcome <-list()
#####mortality in DM overall model####
for (i in c(16:22)) {
  dml1 <-Lexis(entry = list(period = yearin,
                            age    = agein),
               exit = list(period  = out),
               exit.status = dmdata[[i]],
               id  = patid,
               data = dmdata)
  dms1  <-splitMulti(dml1, period= seq(1998,2020,1),
                     age = seq(35,111,1))
  #  db0 = select(dms1, "period","age", "gender", "lex.dur","lex.Xst","lex.id")
  #  write.dta(as.data.frame(db0), file = paste0("lexis", i, ".dta"))
  
  a.kn <- with(subset(dms1, lex.Xst==1), quantile(age+lex.dur,(1:5-0.5)/5))
  p.kn <- with(subset(dms1, lex.Xst==1), quantile(period+lex.dur,(1:5-0.5)/5))
  
  dmc1  <- glm((lex.Xst==1) ~ Ns(period, knots = p.kn) * Ns(age, knots = a.kn),
               family = poisson,
               offset = log(lex.dur),
               data   =dms1)
  
  age    <-c(45:85)
  period <-seq(1998.5,2018.5,1)
  nd1    <-cbind(expand.grid(age, period), lex.dur=1000)
  colnames(nd1)<-c("age","period", "lex.dur")
  p1     <-ci.pred(dmc1, newdata = nd1, Exp = T)
  colnames(p1)<-c("rate","lb","ub")
  outcome[[i]]    <-cbind(nd1, p1,out=i)
}
write.dta(do.call(rbind.data.frame, outcome), file = "DM_res_overall.dta")

#####non-sex-specific mortality in DM by sex model####
rm(outcome)
outcome <-list()
out     <-list()
for (i in c(16:22)) {
  dml1 <-Lexis(entry = list(period = yearin,
                            age    = agein),
               exit = list(period  = out),
               exit.status = dmdata[[i]],
               id  = patid,
               data = dmdata)
  dms1  <-splitMulti(dml1, period= seq(1998,2020,1),
                     age = seq(35,111,1))
  #  db0 = select(dms1, "period","age", "gender", "lex.dur","lex.Xst","lex.id")
  #  write.dta(as.data.frame(db0), file = paste0("lexis", i, ".dta"))
  for (j in c(1:2)) {
    a.kn <- with(subset(dms1, lex.Xst==1 & gender == j), quantile(age+lex.dur,(1:5-0.5)/5))
    p.kn <- with(subset(dms1, lex.Xst==1 & gender == j), quantile(period+lex.dur,(1:5-0.5)/5))
    dmc1  <- glm((lex.Xst==1) ~ Ns(period, knots = p.kn) * Ns(age, knots = a.kn),
                 family = poisson,
                 offset = log(lex.dur),
                 data   = subset(dms1,gender==j))
    
    p1     <-ci.pred(dmc1, newdata = nd1, Exp = T)
    colnames(p1)<-c("rate","lb","ub")
    out[[j]] <-cbind(nd1, p1,out=i, gender=j)
  }
  outcome[[i]] <-data.frame(Reduce(rbind,out))
}
write.dta(do.call(rbind.data.frame, outcome), file = "DM_res_sex.dta")

#####non-sex-specific cancer-mortality in DM by ethnicity model####
rm(outcome)
outcome <-list()
out     <-list()
for (i in c(16:22)) {
  dml1 <-Lexis(entry = list(period = yearin,
                            age    = agein),
               exit = list(period  = out),
               exit.status = dmdata[[i]],
               id  = patid,
               data = dmdata)
  dms1  <-splitMulti(dml1, period= seq(1998,2020,1),
                     age = seq(35,111,1))
  for (j in c(1:2)) {
    a.kn <- with(subset(dms1, lex.Xst==1 & ethnic == j), quantile(age+lex.dur,(1:5-0.5)/5))
    p.kn <- with(subset(dms1, lex.Xst==1 & ethnic == j), quantile(period+lex.dur,(1:5-0.5)/5))
    dmc1  <- glm((lex.Xst==1) ~ Ns(period, knots = p.kn) * Ns(age, knots = a.kn),
                 family = poisson,
                 offset = log(lex.dur),
                 data   = subset(dms1,ethnic==j))
    
    p1     <-ci.pred(dmc1, newdata = nd1, Exp = T)
    colnames(p1)<-c("rate","lb","ub")
    out[[j]]    <-cbind(nd1, p1,out=i, ethnic=j)
  }
  outcome[[i]] <-data.frame(Reduce(rbind,out))
}
write.dta(do.call(rbind.data.frame, outcome), file = "DM_res_ethnic.dta")

#####non-sex-specific cancer-mortality in DM by IMD model####
rm(outcome)
outcome <-list()
out     <-list()
for (i in c(16:22)) {
  dml1 <-Lexis(entry = list(period = yearin,
                            age    = agein),
               exit = list(period  = out),
               exit.status = dmdata[[i]],
               id  = patid,
               data = dmdata)
  dms1  <-splitMulti(dml1, period= seq(1998,2020,1),
                     age = seq(35,111,1))
  for (j in c(1:5)) {
    a.kn <- with(subset(dms1, lex.Xst==1 & imd == j), quantile(age+lex.dur,(1:5-0.5)/5))
    p.kn <- with(subset(dms1, lex.Xst==1 & imd == j), quantile(period+lex.dur,(1:5-0.5)/5))
    dmc1  <- glm((lex.Xst==1) ~ Ns(period, knots = p.kn) * Ns(age, knots = a.kn),
                 family = poisson,
                 offset = log(lex.dur),
                 data   = subset(dms1,imd==j))
    
    p1     <-ci.pred(dmc1, newdata = nd1, Exp = T)
    colnames(p1)<-c("rate","lb","ub")
    out[[j]]    <-cbind(nd1, p1,out=i, imd=j)
  }
  outcome[[i]] <-data.frame(Reduce(rbind,out))
}
write.dta(do.call(rbind.data.frame, outcome), file = "DM_res_imd.dta")

#####non-sex-specific cancer-mortality in DM by smoking model####
rm(outcome)
outcome <-list()
out     <-list()
for (i in c(16:22)) {
  dml1 <-Lexis(entry = list(period = yearin,
                            age    = agein),
               exit = list(period  = out),
               exit.status = dmdata[[i]],
               id  = patid,
               data = dmdata)
  dms1  <-splitMulti(dml1, period= seq(1998,2020,1),
                     age = seq(35,111,1))
  for (j in c(1:3)) {
    a.kn <- with(subset(dms1, lex.Xst==1 & smk_status == j), quantile(age+lex.dur,(1:5-0.5)/5))
    p.kn <- with(subset(dms1, lex.Xst==1 & smk_status == j), quantile(period+lex.dur,(1:5-0.5)/5))
    dmc1  <- glm((lex.Xst==1) ~ Ns(period, knots = p.kn) * Ns(age, knots = a.kn),
                 family = poisson,
                 offset = log(lex.dur),
                 data   = subset(dms1,smk_status==j))
    
    p1     <-ci.pred(dmc1, newdata = nd1, Exp = T)
    colnames(p1)<-c("rate","lb","ub")
    out[[j]]    <-cbind(nd1, p1,out=i, smk_status=j)
  }
  outcome[[i]] <-data.frame(Reduce(rbind,out))
}
write.dta(do.call(rbind.data.frame, outcome), file = "DM_res_smoking.dta")

#####non-sex-specific cancer-mortality in DM by BMI model####
rm(outcome)
outcome <-list()
out     <-list()
for (i in c(16:22)) {
  dml1 <-Lexis(entry = list(period = yearin,
                            age    = agein),
               exit = list(period  = out),
               exit.status = dmdata[[i]],
               id  = patid,
               data = dmdata)
  dms1  <-splitMulti(dml1, period= seq(1998,2020,1),
                     age = seq(35,111,1))
  for (j in c(2:5)) {
    a.kn <- with(subset(dms1, lex.Xst==1 & bmi_grp == j), quantile(age+lex.dur,(1:5-0.5)/5))
    p.kn <- with(subset(dms1, lex.Xst==1 & bmi_grp == j), quantile(period+lex.dur,(1:5-0.5)/5))
    dmc1  <- glm((lex.Xst==1) ~ Ns(period, knots = p.kn) * Ns(age, knots = a.kn),
                 family = poisson,
                 offset = log(lex.dur),
                 data   = subset(dms1,bmi_grp==j))
    
    p1     <-ci.pred(dmc1, newdata = nd1, Exp = T)
    colnames(p1)<-c("rate","lb","ub")
    out[[j]]    <-cbind(nd1, p1,out=i, bmi_grp=j)
  }
  outcome[[i]] <-data.frame(Reduce(rbind,out))
}
write.dta(do.call(rbind.data.frame, outcome), file = "DM_res_bmi.dta")

###trends in all-cause and cancer mortality among diabetes female specific cancer

library(Epi)
library(survival)
library(haven)
library(foreign)
library(broom)
library(popEpi)
library(dplyr)

rm(list=ls())
setwd("C:/Users/SupingLing/Desktop/UoL/KB/paper1_mortality/rerun-Mar2022")
dmdata <-read_dta("dmdb.dta")
nrow(dmdata)
dmdata = subset(dmdata, gender == 2)
attach(dmdata)

outcome <-list()
#####female specific cancer-mortality in DM overall model####
for (i in c(23:24)) {
  dml1 <-Lexis(entry = list(period = yearin,
                            age    = agein),
               exit = list(period  = out),
               exit.status = dmdata[[i]],
               id  = patid,
               data = dmdata)
  dms1  <-splitMulti(dml1, period= seq(1998,2020,1),
                     age = seq(35,111,1))
  #  db0 = select(dms1, "period","age", "gender", "lex.dur","lex.Xst","lex.id")
  #  write.dta(as.data.frame(db0), file = paste0("lexis", i, ".dta"))
  
  a.kn <- with(subset(dms1, lex.Xst==1), quantile(age+lex.dur,(1:5-0.5)/5))
  p.kn <- with(subset(dms1, lex.Xst==1), quantile(period+lex.dur,(1:5-0.5)/5))
  
  dmc1  <- glm((lex.Xst==1) ~ Ns(period, knots = p.kn) * Ns(age, knots = a.kn),
               family = poisson,
               offset = log(lex.dur),
               data   =dms1)
  
  age    <-c(45:85)
  period <-seq(1998.5,2018.5,1)
  nd1    <-cbind(expand.grid(age, period), lex.dur=1000)
  colnames(nd1)<-c("age","period", "lex.dur")
  p1     <-ci.pred(dmc1, newdata = nd1, Exp = T)
  colnames(p1)<-c("rate","lb","ub")
  outcome[[i]]    <-cbind(nd1, p1,out=i)
}
write.dta(do.call(rbind.data.frame, outcome), file = "DM_res_overall2.dta")

#####female specific cancer-mortality in DM by ethnicity model endometrial cannot converge####
rm(outcome)
outcome <-list()
out     <-list()
for (i in c(23:23)) {
  dml1 <-Lexis(entry = list(period = yearin,
                            age    = agein),
               exit = list(period  = out),
               exit.status = dmdata[[i]],
               id  = patid,
               data = dmdata)
  dms1  <-splitMulti(dml1, period= seq(1998,2020,1),
                     age = seq(35,111,1))
  for (j in c(1:2)) {
    a.kn <- with(subset(dms1, lex.Xst==1 & ethnic == j), quantile(age+lex.dur,(1:5-0.5)/5))
    p.kn <- with(subset(dms1, lex.Xst==1 & ethnic == j), quantile(period+lex.dur,(1:5-0.5)/5))
    dmc1  <- glm((lex.Xst==1) ~ Ns(period, knots = p.kn) * Ns(age, knots = a.kn),
                 family = poisson,
                 offset = log(lex.dur),
                 data   = subset(dms1,ethnic==j))
    
    p1     <-ci.pred(dmc1, newdata = nd1, Exp = T)
    colnames(p1)<-c("rate","lb","ub")
    out[[j]]    <-cbind(nd1, p1,out=i, ethnic=j)
  }
  outcome[[i]] <-data.frame(Reduce(rbind,out))
}
write.dta(do.call(rbind.data.frame, outcome), file = "DM_res_ethnic2.dta")

#####female specific cancer cancer-mortality in DM by IMD model####
rm(outcome)
outcome <-list()
out     <-list()
for (i in c(23:24)) {
  dml1 <-Lexis(entry = list(period = yearin,
                            age    = agein),
               exit = list(period  = out),
               exit.status = dmdata[[i]],
               id  = patid,
               data = dmdata)
  dms1  <-splitMulti(dml1, period= seq(1998,2020,1),
                     age = seq(35,111,1))
  for (j in c(1:5)) {
    a.kn <- with(subset(dms1, lex.Xst==1 & imd == j), quantile(age+lex.dur,(1:5-0.5)/5))
    p.kn <- with(subset(dms1, lex.Xst==1 & imd == j), quantile(period+lex.dur,(1:5-0.5)/5))
    dmc1  <- glm((lex.Xst==1) ~ Ns(period, knots = p.kn) * Ns(age, knots = a.kn),
                 family = poisson,
                 offset = log(lex.dur),
                 data   = subset(dms1,imd==j))
    
    p1     <-ci.pred(dmc1, newdata = nd1, Exp = T)
    colnames(p1)<-c("rate","lb","ub")
    out[[j]]    <-cbind(nd1, p1,out=i, imd=j)
  }
  outcome[[i]] <-data.frame(Reduce(rbind,out))
}
write.dta(do.call(rbind.data.frame, outcome), file = "DM_res_imd2.dta")

#####female specific cancer-mortality in DM by smoking model####
rm(outcome)
outcome <-list()
out     <-list()
for (i in c(23:24)) {
  dml1 <-Lexis(entry = list(period = yearin,
                            age    = agein),
               exit = list(period  = out),
               exit.status = dmdata[[i]],
               id  = patid,
               data = dmdata)
  dms1  <-splitMulti(dml1, period= seq(1998,2020,1),
                     age = seq(35,111,1))
  for (j in c(1:3)) {
    a.kn <- with(subset(dms1, lex.Xst==1 & smk_status == j), quantile(age+lex.dur,(1:5-0.5)/5))
    p.kn <- with(subset(dms1, lex.Xst==1 & smk_status == j), quantile(period+lex.dur,(1:5-0.5)/5))
    dmc1  <- glm((lex.Xst==1) ~ Ns(period, knots = p.kn) * Ns(age, knots = a.kn),
                 family = poisson,
                 offset = log(lex.dur),
                 data   = subset(dms1,smk_status==j))
    
    p1     <-ci.pred(dmc1, newdata = nd1, Exp = T)
    colnames(p1)<-c("rate","lb","ub")
    out[[j]]    <-cbind(nd1, p1,out=i, smk_status=j)
  }
  outcome[[i]] <-data.frame(Reduce(rbind,out))
}
write.dta(do.call(rbind.data.frame, outcome), file = "DM_res_smoking2.dta")

#####female specific cancer-mortality in DM by BMI model####
rm(outcome)
outcome <-list()
out     <-list()
for (i in c(23:24)) {
  dml1 <-Lexis(entry = list(period = yearin,
                            age    = agein),
               exit = list(period  = out),
               exit.status = dmdata[[i]],
               id  = patid,
               data = dmdata)
  dms1  <-splitMulti(dml1, period= seq(1998,2020,1),
                     age = seq(35,111,1))
  for (j in c(2:5)) {
    a.kn <- with(subset(dms1, lex.Xst==1 & bmi_grp == j), quantile(age+lex.dur,(1:5-0.5)/5))
    p.kn <- with(subset(dms1, lex.Xst==1 & bmi_grp == j), quantile(period+lex.dur,(1:5-0.5)/5))
    dmc1  <- glm((lex.Xst==1) ~ Ns(period, knots = p.kn) * Ns(age, knots = a.kn),
                 family = poisson,
                 offset = log(lex.dur),
                 data   = subset(dms1,bmi_grp==j))
    
    p1     <-ci.pred(dmc1, newdata = nd1, Exp = T)
    colnames(p1)<-c("rate","lb","ub")
    out[[j]]    <-cbind(nd1, p1,out=i, bmi_grp=j)
  }
  outcome[[i]] <-data.frame(Reduce(rbind,out))
}
write.dta(do.call(rbind.data.frame, outcome), file = "DM_res_bmi2.dta")

###trends in all-cause and cancer mortality among diabetes male specific cancer

library(Epi)
library(survival)
library(haven)
library(foreign)
library(broom)
library(popEpi)
library(dplyr)

rm(list=ls())
setwd("C:/Users/SupingLing/Desktop/UoL/KB/paper1_mortality/rerun-Mar2022")
dmdata <-read_dta("dmdb.dta")
nrow(dmdata)
dmdata = subset(dmdata, gender == 1)
attach(dmdata)

outcome <-list()
##### male specific cancer-mortality in DM overall model####
for (i in c(25:25)) {
  dml1 <-Lexis(entry = list(period = yearin,
                            age    = agein),
               exit = list(period  = out),
               exit.status = dmdata[[i]],
               id  = patid,
               data = dmdata)
  dms1  <-splitMulti(dml1, period= seq(1998,2020,1),
                     age = seq(35,111,1))
  #  db0 = select(dms1, "period","age", "gender", "lex.dur","lex.Xst","lex.id")
  #  write.dta(as.data.frame(db0), file = paste0("lexis", i, ".dta"))
  
  a.kn <- with(subset(dms1, lex.Xst==1), quantile(age+lex.dur,(1:5-0.5)/5))
  p.kn <- with(subset(dms1, lex.Xst==1), quantile(period+lex.dur,(1:5-0.5)/5))
  
  dmc1  <- glm((lex.Xst==1) ~ Ns(period, knots = p.kn) * Ns(age, knots = a.kn),
               family = poisson,
               offset = log(lex.dur),
               data   =dms1)
  
  age    <-c(45:85)
  period <-seq(1998.5,2018.5,1)
  nd1    <-cbind(expand.grid(age, period), lex.dur=1000)
  colnames(nd1)<-c("age","period", "lex.dur")
  p1     <-ci.pred(dmc1, newdata = nd1, Exp = T)
  colnames(p1)<-c("rate","lb","ub")
  outcome[[i]]    <-cbind(nd1, p1,out=i)
}
write.dta(do.call(rbind.data.frame, outcome), file = "DM_res_overall1.dta")

##### male specific cancer-mortality in DM by ethnicity model####
rm(outcome)
outcome <-list()
out     <-list()
for (i in c(25:25)) {
  dml1 <-Lexis(entry = list(period = yearin,
                            age    = agein),
               exit = list(period  = out),
               exit.status = dmdata[[i]],
               id  = patid,
               data = dmdata)
  dms1  <-splitMulti(dml1, period= seq(1998,2020,1),
                     age = seq(35,111,1))
  for (j in c(1:2)) {
    a.kn <- with(subset(dms1, lex.Xst==1 & ethnic == j), quantile(age+lex.dur,(1:5-0.5)/5))
    p.kn <- with(subset(dms1, lex.Xst==1 & ethnic == j), quantile(period+lex.dur,(1:5-0.5)/5))
    dmc1  <- glm((lex.Xst==1) ~ Ns(period, knots = p.kn) * Ns(age, knots = a.kn),
                 family = poisson,
                 offset = log(lex.dur),
                 data   = subset(dms1,ethnic==j))
    
    p1     <-ci.pred(dmc1, newdata = nd1, Exp = T)
    colnames(p1)<-c("rate","lb","ub")
    out[[j]]    <-cbind(nd1, p1,out=i, ethnic=j)
  }
  outcome[[i]] <-data.frame(Reduce(rbind,out))
}
write.dta(do.call(rbind.data.frame, outcome), file = "DM_res_ethnic1.dta")

##### male specific cancer-mortality in DM by IMD model####
rm(outcome)
outcome <-list()
out     <-list()
for (i in c(25:25)) {
  dml1 <-Lexis(entry = list(period = yearin,
                            age    = agein),
               exit = list(period  = out),
               exit.status = dmdata[[i]],
               id  = patid,
               data = dmdata)
  dms1  <-splitMulti(dml1, period= seq(1998,2020,1),
                     age = seq(35,111,1))
  for (j in c(1:5)) {
    a.kn <- with(subset(dms1, lex.Xst==1 & imd == j), quantile(age+lex.dur,(1:5-0.5)/5))
    p.kn <- with(subset(dms1, lex.Xst==1 & imd == j), quantile(period+lex.dur,(1:5-0.5)/5))
    dmc1  <- glm((lex.Xst==1) ~ Ns(period, knots = p.kn) * Ns(age, knots = a.kn),
                 family = poisson,
                 offset = log(lex.dur),
                 data   = subset(dms1,imd==j))
    
    p1     <-ci.pred(dmc1, newdata = nd1, Exp = T)
    colnames(p1)<-c("rate","lb","ub")
    out[[j]]    <-cbind(nd1, p1,out=i, imd=j)
  }
  outcome[[i]] <-data.frame(Reduce(rbind,out))
}
write.dta(do.call(rbind.data.frame, outcome), file = "DM_res_imd1.dta")

##### male specific cancer-mortality in DM by smoking model####
rm(outcome)
outcome <-list()
out     <-list()
for (i in c(25:25)) {
  dml1 <-Lexis(entry = list(period = yearin,
                            age    = agein),
               exit = list(period  = out),
               exit.status = dmdata[[i]],
               id  = patid,
               data = dmdata)
  dms1  <-splitMulti(dml1, period= seq(1998,2020,1),
                     age = seq(35,111,1))
  for (j in c(1:3)) {
    a.kn <- with(subset(dms1, lex.Xst==1 & smk_status == j), quantile(age+lex.dur,(1:5-0.5)/5))
    p.kn <- with(subset(dms1, lex.Xst==1 & smk_status == j), quantile(period+lex.dur,(1:5-0.5)/5))
    dmc1  <- glm((lex.Xst==1) ~ Ns(period, knots = p.kn) * Ns(age, knots = a.kn),
                 family = poisson,
                 offset = log(lex.dur),
                 data   = subset(dms1,smk_status==j))
    
    p1     <-ci.pred(dmc1, newdata = nd1, Exp = T)
    colnames(p1)<-c("rate","lb","ub")
    out[[j]]    <-cbind(nd1, p1,out=i, smk_status=j)
  }
  outcome[[i]] <-data.frame(Reduce(rbind,out))
}
write.dta(do.call(rbind.data.frame, outcome), file = "DM_res_smoking1.dta")

##### male specific cancer-mortality in DM by BMI model####
rm(outcome)
outcome <-list()
out     <-list()
for (i in c(25:25)) {
  dml1 <-Lexis(entry = list(period = yearin,
                            age    = agein),
               exit = list(period  = out),
               exit.status = dmdata[[i]],
               id  = patid,
               data = dmdata)
  dms1  <-splitMulti(dml1, period= seq(1998,2020,1),
                     age = seq(35,111,1))
  for (j in c(2:5)) {
    a.kn <- with(subset(dms1, lex.Xst==1 & bmi_grp == j), quantile(age+lex.dur,(1:5-0.5)/5))
    p.kn <- with(subset(dms1, lex.Xst==1 & bmi_grp == j), quantile(period+lex.dur,(1:5-0.5)/5))
    dmc1  <- glm((lex.Xst==1) ~ Ns(period, knots = p.kn) * Ns(age, knots = a.kn),
                 family = poisson,
                 offset = log(lex.dur),
                 data   = subset(dms1,bmi_grp==j))
    
    p1     <-ci.pred(dmc1, newdata = nd1, Exp = T)
    colnames(p1)<-c("rate","lb","ub")
    out[[j]]    <-cbind(nd1, p1,out=i, bmi_grp=j)
  }
  outcome[[i]] <-data.frame(Reduce(rbind,out))
}
write.dta(do.call(rbind.data.frame, outcome), file = "DM_res_bmi1.dta")


