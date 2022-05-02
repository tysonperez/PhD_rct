#Objective 1: 
#To compare the effectiveness of genuine versus sham ISF-NFB for treating IDs in a female population. 
#We will assess differences between genuine 1-region, genuine 2-region, & sham ISF-NFB in PROs & neurophysiological measures after 6 sessions.  
#We hypothesize that all groups will show clinical improvements via non-specific treatment effects (e.g. placebo) but genuine groups 
#will demonstrate additional improvements due to specific effects (i.e. effects specific to the modulation of the targeted EEG variable(s)).

#Load libraries
  library(rjags)  
  library(readxl)   
  library(tidyverse)
  library(bayesplot)
  library(bayestestR)
  library(see)
  library(patchwork)

#Load & visualise data
  clinical_data <- read_excel("Clinical.xlsx", sheet="Obj1-5")
  clinical_data 

#Subset data = Participant, Treatment, baseline #2 = t1 & post-6 sessions = t2  
  get_data_obj1 <- function(data){
    #update last two column numbers in data[c(...)] to reflect outcome of interest
    tmp = data[c(1,2,35,36)]
    tmp$tx =  
      1*(tmp$Treatment=="control") +
      2*(tmp$Treatment=="isf1") +
      3*(tmp$Treatment=="isf2")
    return(tmp)
  }
  data_obj1 <- get_data_obj1(clinical_data) 
  data_obj1
  
  #remove rows with NAs
  data_obj1 <- data_obj1[!is.na(data_obj1[,3])&!is.na(data_obj1[,4]),]
  data_obj1
  
#Define the model (multivariate normal)
  m1 <- "model{
      for(i in 1:n){                            # n = participant number
          y[i,1:2] ~ dmnorm(mu[i,],omega[,])    # 1:2 = time1:time2 
          mu[i,1] <- mu1                        # time1 is baseline 
          mu[i,2] <- mu1 + cum.effect[tx[i]]    # time2 cumulative effects = nonspecific + tx effects
          
          #predicted data 
          y.pred[i,1:2] ~ dmnorm(mu.pred[i,],omega[,])
          mu.pred[i,1] <- mu1
          mu.pred[i,2] <- mu1 + cum.effect[tx[i]]
          
          #calculate Pearson residuals^2 for observed & predicted data 
          pearson[i,1] <- (y[i,1] - mu[i,1])/sqrt(sigma[1,1])
          pearson[i,2] <- (y[i,2] - mu[i,2])/sqrt(sigma[2,2])
          pearson.pred[i,1] <- (y.pred[i,1] - mu.pred[i,1])/sqrt(sigma[1,1])
          pearson.pred[i,2] <- (y.pred[i,2] - mu.pred[i,2])/sqrt(sigma[2,2])
          D1[i] <- pow(pearson[i,1],2)
          D2[i] <- pow(pearson[i,2],2)
          D1.pred[i] <- pow(pearson.pred[i,1],2)
          D2.pred[i] <- pow(pearson.pred[i,2],2)
          }
      
    #vague priors 
    mu1 ~ dnorm(10,0.0001)
    
    for(i in 1:3){
        cum.effect[i] ~ dnorm(0,0.0001)
    }
        
    omega[1:2,1:2] ~ dwish(R[,],2)
		    sigma[1:2,1:2] <- inverse(omega[,])
		    R[1,1] <- 0.1
		    R[2,2] <- 0.1
		    R[1,2] <- 0.005
		    R[2,1] <- 0.005
	      rho <- sigma[1,2]/sqrt(sigma[1,1]*sigma[2,2])
	
    #nodes to monitor
	  baseline <- mu1
	  post_sham <- mu1 + cum.effect[1]
	  post_isf1 <- mu1 + cum.effect[2]
	  post_isf2 <- mu1 + cum.effect[3]
	  delta_sham <- cum.effect[1]
	  delta_isf1 <- cum.effect[2]
	  delta_isf2 <- cum.effect[3]
	  isf1_v_sham <- cum.effect[2] - cum.effect[1]
	  isf2_v_sham <- cum.effect[3] - cum.effect[1]
	  isf2_v_isf1 <- cum.effect[3] - cum.effect[2]
	  std_sham <- cum.effect[1]/sqrt(sigma[1,1]+sigma[2,2]-2*sigma[1,2])
	  std_isf1 <- cum.effect[2]/sqrt(sigma[1,1]+sigma[2,2]-2*sigma[1,2])
  	std_isf2 <- cum.effect[3]/sqrt(sigma[1,1]+sigma[2,2]-2*sigma[1,2])
  	std_isf1_v_sham <- (cum.effect[2]-cum.effect[1])/sqrt(sigma[1,1]+sigma[2,2]-2*sigma[1,2])
  	std_isf2_v_sham <- (cum.effect[3]-cum.effect[1])/sqrt(sigma[1,1]+sigma[2,2]-2*sigma[1,2])
  	mcid <- -1.5/sqrt(sigma[1,1]+sigma[2,2]-2*sigma[1,2])
  	
  	#summary stat (T) = sum the Pearson res^2 for observed & predicted data
  	T1 <- sum(D1[])
  	T1.pred <- sum(D1.pred[])
  	T2 <- sum(D2[])
  	T2.pred <- sum(D2.pred[])
  	
  	#calculate probability T.pred >= T (Bayesian p-value)
  	Bayes.p1 <- step(T1.pred - T1)
  	Bayes.p2 <- step(T2.pred - T2)
    }"

#Generate initial values
  inits_func <- function(){
    mu1 <- runif(1,5,15)
    return(list(mu1=mu1))
  }

#Compile the model
  m1_jags <- jags.model(textConnection(m1),
                        data=list(y=cbind(data_obj1[,3],data_obj1[,4]),
                                  n=length(data_obj1$Participant),
                                  tx=data_obj1$tx),
                        inits=inits_func,
                        n.chains = 3,
                        n.adapt=5000)

#Burn in
  update(m1_jags,10000)

#Simulate the posterior
  m1_out <- coda.samples(model=m1_jags,
                         variable.names=c("rho","sigma[1,1]","sigma[1,2]","sigma[2,2]",
                                          "baseline",
                                          "post_sham","post_isf1","post_isf2",
                                          "delta_sham","delta_isf1","delta_isf2",
                                          "isf1_v_sham","isf2_v_sham","isf2_v_isf1",
                                          "std_sham","std_isf1","std_isf2",
                                          "std_isf1_v_sham","std_isf2_v_sham",
                                          "mcid",
                                          "T1","T1.pred","Bayes.p1",
                                          "T2","T2.pred","Bayes.p2"),
                         n.iter=25000)

#####Posteriors######

  #plot posteriors
  post_plot <- function(data,parameter){
    color_scheme_set("brightblue")
    mcmc_areas(data, 
               prob=0.95, 
               point_est="mean",
               pars=parameter) +
      labs(title="Posterior distributions", 
           subtitle="with means & 95% HDI")
  }

   # post_plot(m1_out,c("delta_sham","delta_isf1","isf1_v_sham"))
   # post_plot(m1_out,c("delta_sham","delta_isf2","isf2_v_sham"))
    post_plot(m1_out,c("std_sham","std_isf1","std_isf1_v_sham"))
    post_plot(m1_out,c("std_sham","std_isf2","std_isf2_v_sham"))
   
  #summary of posterior
  summary(m1_out)   

  
#####Standardized Effect size#####
  
  #describe the posterior
  describe_posterior(m1_out[,c(22:26)],         # [,c(#,#)] parameters of interest
                     centrality = "mean",
                     ci_method="hdi",           
                     rope_ci = 1,
                     rope_range = c(-0.1,0.1))
  
  #convert model object to matrix
  m1_out_matrix <- as.matrix(m1_out[,c(22:26)]) # [,c(#,#)] parameters of interest
  head(m1_out_matrix)
  
  
  #plot probability of direction
  plot(p_direction(m1_out_matrix[,4], ci=1))      # [,#] parameter of interest
  
  #plot % in rope
  plot(rope(m1_out_matrix[,4],                    # [,#] parameter of interest
            ci=1,
            ci_method="hdi",
            range=c(-0.1,0.1)))             
    
####Unstandardized Effect Size####
  
  # #describe the posterior
  # describe_posterior(m1_out[,c(8:10,16,18)],         # [,c(#,#)] parameters of interest
  #                    centrality = "mean",
  #                    ci_method="hdi",           
  #                    rope_ci = 1,
  #                    rope_range = c(-1.5,0))
  # 
  # #convert model object to matrix
  # m1_out_matrix <- as.matrix(m1_out[,c(8:10,16,18)]) # [,c(#,#)] parameters of interest
  # head(m1_out_matrix)
  # 
  # #describe the posterior
  # describe_posterior(m1_out[,c(8:10,16,18)],         # [,c(#,#)] parameters of interest
  #                    centrality = "mean",
  #                    ci_method="hdi",           
  #                    rope_ci = 1,
  #                    rope_range = c(-1.5,0))
  
  #######################DIAGNOSTICS#######################
  
  ######Convergence check##########
  
  #trace plots (look for fuzzy caterpillar)
  color_scheme_set("mix-blue-red")
  mcmc_trace(m1_out, 
             pars=c("std_isf1_v_sham","std_isf2_v_sham"),
             facet_args=list(ncol=1,strip.position="left"))
  
  #####Posterior predictive check of model fit#########
  
  #want mean Bayesian p-value around 0.5
  summary(m1_out)  
  
  #want points to be distributed equally on both sides of line
  m1_out_matrix <- as.matrix(m1_out)
  m1_out_matrix
  plot(m1_out_matrix[,3],m1_out_matrix[,4],
       xlim=c(0,100),ylim=c(0,100),
       xlab="T1",ylab="T1.pred",pch=".")
  lines(c(0,100),c(0,100))
  
  plot(m1_out_matrix[,5],m1_out_matrix[,6],
       xlim=c(0,100),ylim=c(0,100),
       xlab="T2",ylab="T2.pred",pch=".")
  lines(c(0,100),c(0,100))
  
  #####Residuals#######
  
  #Baseline (time 1): look for anything bizarre e.g. multi-modal but wouldn't expect it
  mu1 = summary(m1_out)[[1]][7,1]  #[[1]] = 1st of 2 items in the m1_out list
  mu1
  s11 = summary(m1_out)[[1]][18,1]  
  s11
  
  #Raw & standardised residuals at baseline
  res1 = data_obj1[,3][[1]] - mu1 
  res1_std = (data_obj1[,3][[1]] - mu1)/sqrt(s11)  
  
  #Residuals plot: look for outliers, patterns, changes in spread
  plot(seq(length(res1)),res1, ylim=c(-20, 20))  
  plot(seq(length(res1_std)),res1_std, ylim=c(-3,3))
  
  #QQ-plot: normality doesn't matter much, but look for extreme departures
  qqnorm(res1_std, main="QQ Plot: Baseline"); qqline(res1_std)
  
  #Post 6 sessions (time 2)
  
  #Subset data by group
  df_group <- function(data,treatment){
    df <- data %>% filter(Treatment==treatment)
    return(df)
  }
  df_sham <- df_group(data_obj1, "control")
  df_isf1 <- df_group(data_obj1, "isf1")
  df_isf2 <- df_group(data_obj1, "isf2")
  
  #select mean inferred value
  mu2_sham <- summary(m1_out)[[1]][16,1]
  mu2_sham
  mu2_isf1 <- summary(m1_out)[[1]][14,1]
  mu2_isf1
  mu2_isf2 <- summary(m1_out)[[1]][15,1]
  mu2_isf2
  s22 <- summary(m1_out)[[1]][20,1]
  s22
  
  #Check assumption of constant variance of residuals 
  res2 <- function(df, mu2, group){
    res <- df[,4][[1]] - mu2
    plot <- plot(seq(length(res)),
                 res, 
                 ylim=c(-40, 40),
                 main=group,
                 ylab="residuals",
                 xlab="n") + abline(h=0,col="red")
    return(res)           
  }
  res2_sham <- res2(df_sham, mu2_sham, "SHAM")
  res2_isf1 <- res2(df_isf1, mu2_isf1, "ISF1")
  res2_isf2 <- res2(df_isf2, mu2_isf2, "ISF2")
  
  #Check assumption of constant variance of standardized residuals 
  res2_std <- function(df, mu2, group){
    res_std = (df[,4][[1]] - mu2)/sqrt(s22)
    plot <- plot(seq(length(res_std)),
                 res_std, ylim=c(-4,4),
                 main=group,
                 ylab="standardized residuals",
                 xlab="n") + abline(h=0,col="red")
    return(res_std)
  }
  res2_std_sham <- res2_std(df_sham,mu2_sham,"Sham")
  res2_std_isf1 <- res2_std(df_isf1,mu2_isf1,"ISF1")
  res2_std_isf2 <- res2_std(df_isf2,mu2_isf2,"ISF2")
  
  #Check assumption of normally distributed residuals
  qqnorm(res2_std_sham,main="QQ Plot: Sham")
         #ylab="standardized residuals",
         #xlab="Inverse normal"); 
  qqline(res2_std_sham)
  qqnorm(res2_std_isf1,main="QQ Plot: ISF1")
         #ylab="standardized residuals",
         #xlab="Inverse normal"); 
  qqline(res2_std_isf1)
  qqnorm(res2_std_isf2,main="QQ Plot: ISF2")
         #ylab="standardized residuals",
         #xlab="Inverse normal"); 
  qqline(res2_std_isf2)
  
  #Check assumption of equal variances met
  res_vert <- function(df,data){
    res <- ggplot(df,aes(x=Treatment,y=data)) + 
      geom_point(alpha=.5) + 
      ylab("standardized residuals") + ylim(-4,4)
  }
  res_sham <- res_vert(df_sham,res2_std_sham)
  res_isf1 <- res_vert(df_isf1,res2_std_isf1)
  res_isf2 <- res_vert(df_isf2,res2_std_isf2)
  
  compare.plot <- res_sham + res_isf1 + res_isf2
  print(compare.plot)
  
