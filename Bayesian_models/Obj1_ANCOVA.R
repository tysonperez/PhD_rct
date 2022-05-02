#Objective 1: 
#To compare the effectiveness of genuine versus sham ISF-NFB for treating IDs in a female population. 
#We will assess differences between genuine 1-region, genuine 2-region, & sham ISF-NFB in PROs & neurophysiological measures after 6 sessions.  
#We hypothesise that all groups will show clinical improvements via non-specific treatment effects (e.g. placebo) but genuine groups 
#will demonstrate additional improvements due to specific effects (i.e. effects specific to the modulation of the targeted EEG variable(s)).

#Load libraries
library(rjags)  
library(readxl)   
library(tidyverse)
library(bayesplot)
library(patchwork)

#Load & visualise data
clinical_data <- read_excel("Clinical.xlsx", sheet="Obj1-5")
clinical_data

#Subset data = Participant, Treatment, baseline 2=t1 & midtreatment=t2 columns 
get_data_obj1 <- function(data){
  #update last two column numbers in data[c(...)] to reflect outcome of interest
  tmp = data[c(1,2,12,13)]
  tmp$tx =  
    1*(tmp$Treatment=="control") +
    2*(tmp$Treatment=="isf1") +
    3*(tmp$Treatment=="isf2")
  return(tmp)
}
data_obj1 <- get_data_obj1(clinical_data)
data_obj1

#Define the model

  m1 <- "model{
          for(i in 1:n.p){                        #n.p=number of participants
            y6[i] ~ dnorm(mu.y6[i],tau[tx[i]])                       
            mu.y6[i] <-   mu +
                          treat.effect[tx[i]] + 
                          beta.baseline*y0[i]
        }
        
        mu <- 0                                   #prevents overparameterisation & infinite number of solutions
        
        #vague priors
        beta.baseline ~ dnorm(0,0.0001)
        for(i in 1:3){
          treat.effect[i] ~ dnorm(0,0.0001)
          tau[i] <- 1/(sd[i]*sd[i])
          sd[i] ~ dt(0,0.033,3)T(0,)              #half t distribution 
          }
        
        #nodes to monitor
        isf1_v_sham <- treat.effect[2]- treat.effect[1]
        isf2_v_sham <- treat.effect[3]- treat.effect[1]
        #isf2_v_isf1 <- treat.effect[3]- treat.effect[2] 
        
  }"

#Generate initial values
inits_func <- function(){
  sd <- runif(3,0,10)
  treat.effect <- runif(3,-10,10)
  beta.baseline <- runif(1,-5,5)
  return(list(sd=sd,treat.effect=treat.effect,beta.baseline=beta.baseline))
}

#Compile the model

  m1_jags <- jags.model(textConnection(m1),
                        data=list(y6=data_obj1[,4][[1]],
                                  y0=data_obj1[,3][[1]],
                                  n.p=length(data_obj1$Participant),
                                  tx=data_obj1$tx),
                        inits=inits_func,
                        n.chains = 3,
                        n.adapt=5000)
#Burn in
update(m1_jags,10000)

#Simulate the posterior

m1_out <- coda.samples(model=m1_jags,
                       variable.names=c("isf1_v_sham","isf2_v_sham"),
                       n.iter=25000)

######Posteriors######

#numerical summary
summary(m1_out)

#plot posteriors
post_plot <- function(data,parameter){
  color_scheme_set("brightblue")
  mcmc_areas(data, 
             prob=0.95, 
             point_est="mean",
             pars=parameter) +
    labs(title="Posterior distributions", 
         subtitle="with means & 95% credible intervals")
}
post_plot(m1_out,c("isf1_v_sham","isf2_v_sham"))



