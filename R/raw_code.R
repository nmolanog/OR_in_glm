rm(list=ls())
library(pacman)
library(pacman)
p_load(here)
p_load(tidyverse)
p_load(reshape2)
p_load(boot)
p_load(epitools)
p_load(kableExtra)
p_load(knitr)
p_load(ggmosaic)
set.seed(150)

# Introduction
ca_ctr_r<-.3
n<-250
nCA<-round(n*ca_ctr_r)
z0<-data.frame(status=c(rep("flu",nCA),rep("no_flu",n-nCA)))
z0$vac<-NA
exp_CA<-.35
exp_CTR<-.75
z0[z0$status %in% "flu","vac"]<-ifelse(runif(nCA)<exp_CA,"yes","no")
z0[z0$status %in% "no_flu","vac"]<-ifelse(runif(n-nCA)<exp_CTR,"yes","no")
z0$vac<-factor(z0$vac,levels = c("yes","no"))

#contingency table of flu vs vac
t(table(z0)%>%addmargins)
#row proportions of flu vs vac  
t(table(z0))%>%prop.table(1)

or1<-t(table(z0))%>%oddsratio
or2<-table(z0)%>%oddsratio

#relevel
z0$status<-factor(z0$status,levels = c("no_flu","flu"))
or3<-t(table(z0))%>%oddsratio

or1$measure
or2$measure
or3$measure

# Odds ratios for $K \times 2$ tables
N<-1000

treatment<-sample(LETTERS[1:3],N,replace =T)

cons<-NA
cons[treatment %in% "A"]<-rnorm(sum(treatment %in% "A"),100,10)
cons[treatment %in% "B"]<-rnorm(sum(treatment %in% "B"),100,10)
cons[treatment %in% "C"]<-rnorm(sum(treatment %in% "C"),100,10)

sim_probs<-rep(NA,N)

eqmA<-matrix(c(1,min(cons),
               1,max(cons)),byrow = T,nrow=2)
eqbA<-qlogis(c(.01,.6))
Abetas<-solve(eqmA,eqbA)

eqmB<-matrix(c(1,min(cons),
               1,max(cons)),byrow = T,nrow=2)
eqbB<-qlogis(c(.3,.9))
Bbetas<-solve(eqmB,eqbB)

eqmC<-matrix(c(1,min(cons),
               1,max(cons)),byrow = T,nrow=2)
eqbC<-qlogis(c(.99,.4))
Cbetas<-solve(eqmC,eqbC)

sim_probs[treatment %in% "A"]<-inv.logit(Abetas[1]+Abetas[2]*cons[treatment %in% "A"])
sim_probs[treatment %in% "B"]<-inv.logit(Bbetas[1]+Bbetas[2]*cons[treatment %in% "B"])
sim_probs[treatment %in% "C"]<-inv.logit(Cbetas[1]+Cbetas[2]*cons[treatment %in% "C"])

z1<-data.frame(treatment,cons,sim_probs,stringsAsFactors = T)
z1$BM<-map_dbl(z1$sim_probs,~rbinom(1,1,.))
z1$BM<-factor(z1$BM,level=0:1,labels = c("No","Yes"))

z1$BM<-relevel(z1$BM,ref="Yes")

#contingency table of treatment vs behavior modifications
table(z1[,c("treatment","BM")])%>%addmargins
#OR's for reference treatment A
oddsratio(table(z1$treatment,z1$BM))

ggplot(data = z1) +
  geom_mosaic(aes(x = product(BM,treatment), fill=BM), na.rm=TRUE)+
  labs(x = "treatment",y="")+theme_bw()

##logistic regresion
#case 1
z0$vac<-relevel(z0$vac,ref="no")
summary(z0)

m0<-glm(status~vac,data = z0,family = binomial)
summary(m0)
#compare OR
coef(m0) %>% exp
or1

predict(m0,data.frame(vac="yes"),type ="response")
predict(m0,data.frame(vac="no"),type ="response")
t(table(z0))%>%prop.table(1)

ggplot(data = z0) +
  geom_mosaic(aes(x = product(status,vac), fill=status), na.rm=TRUE)+
  labs(x = "treatment",y="")+theme_bw()

#case 2
z1$BM<-relevel(z1$BM,ref="No")
summary(z1)

m1<-glm(BM~treatment,data = z1,family = binomial)
summary(m1)
#inverted OR
coef(m1) %>% exp
# OR as in oddsratio function
coef(m1) %>% {exp(-.)}

oddsratio(table(z1$treatment,z1$BM),rev ="columns")
