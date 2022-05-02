#Objective 2: 
#To assess whether greater changes in EEG metrics 
#are associated with greater changes in PROs.  In particular, we are interested
#in the correlation between PRO score improvements and modulation of EEG variables 
#(i.e. activity & connectivity) in the targeted ROIs (MCC & PCC) 

#Load libraries
library(readxl)   
library(tidyverse)
library(bayesplot)
library(rstanarm)
library(bayestestR)
library(see)

#Load & visualize data
clinical_data <- read_excel("Clinical_PN.xlsx")#,sheet="Obj1-5")
clinical_data

#Gather relevant data & add delta variables 
get_data_obj2 <- function(data,changecolumns){
  #changecolumns are the columns used to calculate the change scores
  tmp = data[c(1,2,changecolumns)]
  tmp$delta.score = (tmp[,4]-tmp[,3])[[1]]
  tmp$delta.R1 = (tmp[,6]-tmp[,5])[[1]]
  return(tmp)
}

#update column #s in c(...) according to data of interest
data_obj2 <- get_data_obj2(clinical_data,c(7,8,35,36))
data_obj2

#remove rows with NAs
data_obj2 <- data_obj2[!is.na(data_obj2[,3])&!is.na(data_obj2[,4])&
                       !is.na(data_obj2[,5])&!is.na(data_obj2[,6]),]
data_obj2

######Regression plot of 500 random draws from posterior distribution using non-standardized variables######

#regress model
reg2 <- stan_glm(delta.score ~ delta.R1, data=data_obj2,
                     chains=3, iter=35000, warmup=10000)
summary(reg2)

#coerce model to data-frame & remove sigma
reg2_df <- reg2 %>% as_tibble() %>% rename(intercept=`(Intercept)`) %>% select(-sigma)
head(reg2_df)

#aesthetics
n_draws <- 500
alpha_level <- .15
color_draw <- "grey60"
color_mean <- "#3366FF"  

#plot
ggplot(data_obj2) +
  aes(x=delta.R1,y=delta.score) +
  coord_cartesian(ylim=c(-12,6)) +                       #restrict y-axis to focus on centre of data
  geom_abline(aes(intercept=intercept,slope=delta.R1),  #plot random sample rows from sim df as grey,semi-transparent lines
              data=sample_n(reg2_df,n_draws),
              color=color_draw,
              alpha=alpha_level) +
  geom_abline(intercept=mean(reg2_df$intercept),        #plot mean of intercept & slope coefficients in blue
              slope=mean(reg2_df$delta.R1),
              size=1,
              color=color_mean) +
  geom_point(aes(color=Treatment)) + 
  labs(x="?? PN's sgACC slow (0.2-1.5 Hz) log-CSD",
       y="?? HADS-D",
       title="Visualization of 500 Regression Lines from the Posterior Distribution")

######Standardized regression coefficient######

#standardize the variables
data_obj2_std <- data.frame(scale(data_obj2[,7:8]))
data_obj2_std

#regress model
reg2_std <- stan_glm(delta.score ~ delta.R1, data=data_obj2_std,
                     chains=3, iter=35000, warmup=10000)
summary(reg2_std)

#extract std reg coefficient from model object
b1 <- as.matrix(reg2_std)
head(b1)

#plot posterior distribution of std slope
post_plot <- function(data,parameter){
  color_scheme_set("brightblue")
  mcmc_areas(data, 
             prob=0.95, 
             point_est="mean",
             pars=parameter) +
    labs(title="Posterior distributions", 
         subtitle="with means & 95% HDI")
}
post_plot(b1,c("delta.R1"))

#describe the posterior of std slope
describe_posterior(reg2_std,
                   centrality = "mean",
                   ci_method="hdi",
                   rope_ci = 1,
                   rope_range = c(-0.05,0.05))

#plot rope of std slope
plot(rope(reg2_std, ci=1))

# #####Correlation coefficient####
# 
# #standardize the variables
# data_obj2_std <- data.frame(scale(data_obj2[,7:8]))
# data_obj2_std
# 
# cor <- correlationBF(data_obj2_std$delta.score, data_obj2_std$delta.R1)
# describe_posterior(cor,
#                    centrality = "mean",
#                    ci_method="hdi",
#                    rope_ci = 1,
#                    rope_range = c(-0.05,0.05))
# plot(rope(cor, ci=1))
# plot(p_direction(cor, ci=1)) 

#####Diagnostics#######

summary(reg2)           #Rhat & effective sample size
pp_check(reg2, "stat")  #posterior predictive distribution


