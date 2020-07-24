#######################
###load data
#######################
rm(list=ls())
options(max.print=999999)
library(pacman)
p_load(here)
p_load(tidyverse)
p_load(reshape2)
p_load(boot)
p_load(epitools)
library(bueri)
set.seed(150)
N<-1000

treatment<-sample(LETTERS[1:3],N,replace =T)

cons<-NA
cons[treatment %in% "A"]<-rnorm(sum(treatment %in% "A"),100,10)
cons[treatment %in% "B"]<-rnorm(sum(treatment %in% "B"),100,10)
cons[treatment %in% "C"]<-rnorm(sum(treatment %in% "C"),100,10)

summary(cons)
  

sim_probs<-rep(NA,N)

eqmA<-matrix(c(1,min(cons),
               1,max(cons)),byrow = T,nrow=2)
eqbA<-qlogis(c(.01,.6))
(Abetas<-solve(eqmA,eqbA))

eqmB<-matrix(c(1,min(cons),
               1,max(cons)),byrow = T,nrow=2)
eqbB<-qlogis(c(.3,.9))
(Bbetas<-solve(eqmB,eqbB))

eqmC<-matrix(c(1,min(cons),
               1,max(cons)),byrow = T,nrow=2)
eqbC<-qlogis(c(.99,.4))
(Cbetas<-solve(eqmC,eqbC))

sim_probs[treatment %in% "A"]<-inv.logit(Abetas[1]+Abetas[2]*cons[treatment %in% "A"])
sim_probs[treatment %in% "B"]<-inv.logit(Bbetas[1]+Bbetas[2]*cons[treatment %in% "B"])
sim_probs[treatment %in% "C"]<-inv.logit(Cbetas[1]+Cbetas[2]*cons[treatment %in% "C"])

z0<-data.frame(treatment,cons,sim_probs,stringsAsFactors = T)
z0$outcome<-map_dbl(z0$sim_probs,~rbinom(1,1,.))
z0$outcome<-factor(z0$outcome)
summary(z0)
table(z0$treatment,z0$outcome)%>%prop.table(1)
oddsratio(table(z0$treatment,z0$outcome),rev="columns")

m0<-glm(outcome~treatment,data=z0,family = "binomial")
summary(m0)

m01<-glm(outcome~treatment*cons,data=z0,family = "binomial")
summary(m01)



