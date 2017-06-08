

######################################################################
#                                                                    #
#                           Load packages                            #
#                                                                    #
######################################################################
library(BayesFactor)
library(ggplot2)

######################################################################
#                                                                    #
#                            Functions                               #
#                                                                    #
######################################################################

######################################################################
#                Morey functions to allow for custom prior           #
######################################################################

# These functions are required to run the bayes factor analysis with custom priors
# Richard Morey provided these functions to make this possible

#Functions source: https://gist.github.com/richarddmorey/e49c345fcd6cfe8f32535eda4bd8c29b
#also availible with:
#devtools::source_gist("e49c345fcd6cfe8f32535eda4bd8c29b", filename='utility_functions.R')

fullAlt_norm = Vectorize(function(theta, y, N, shift, scale){
  p = plogis(theta)
  exp(dbinom(y, N, p, log = TRUE) + dnorm(theta, shift, scale * pi/sqrt(3), log = TRUE))
},"theta")

normalize_norm = function(shift, scale, interval){
  diff(pnorm(interval, shift, scale * pi/sqrt(3)))
}

restrictedAlt_norm = function(theta,y,N,shift,scale,interval){
  fullAlt_norm(theta,y,N,shift,scale) / normalize_norm(shift, scale, interval) * (theta>interval[1] & theta<interval[2])
}

margLike_norm = function(y, N, shift, scale, interval){
  theta_interval = qlogis(sort(interval))
  integrate(restrictedAlt_norm, theta_interval[1], theta_interval[2], 
            y = y, N = N, shift = shift, scale = scale, interval = theta_interval)[[1]]
}

BF01_norm = Vectorize(function(y, N, shift, scale, interval, null.p){
  dbinom(y,N,null.p) / margLike_norm(y, N, shift, scale, interval)
},"y")



######################################################################
#                         Simulation function                        #
######################################################################

# function to simulate data and data analysis
# the function is used for power analysis

simul <- function(H0_prob, H1_prob, min_num_trials, max_num_trials,
                  rho_scale_proptest, one_sided, prior_distribution, 
                  stop_when_BF_low, stop_when_BF_high){
  
  # one or two sided tests
  # if one_sided == T a one-sided test will be performed 
  if(one_sided == T){
    proptest_alternative = "less" # alternative hypothesis is "less" because the default level is the first cell of the proportion table, which is the 0, in the alternative hypothesis we expect that there will be less 0-s (representing misses) then 1-s (representing hits) 
    alter_Interval_proptest = c(0.5,1)
  } else {
    proptest_alternative = "two.sided"
    alter_Interval_proptest = c(0,1)
  }
  
  # determine if there is optional stopping based on BF or not
  # if no criteria are set for optional stopping, the minimum and maximum number of trials are set to be equal, making the simulation run until reaching the maximum number of trials
  if(is.na(stop_when_BF_low) & is.na(stop_when_BF_high)){min_num_trials = max_num_trials}
  
  ### Simulation assuming H1 is true
  
  # generate data asuming H1 is true, the simulated data asumes each trial to be independent
  # see reasoning on this decision here: https://osf.io/rh5hb/wiki/Reasoning%20behind%20data%20analysis%20decisions/
  data_all_H1 <- rbinom(max_num_trials, size = 1, prob=H1_prob)
  
  # compute bayes factor on the simulated data where H1 is true
  # bayes factor is computed starting with a part of the dataset including the minimum sample size set in the function parameters and this is repeated in a loop increasing the dataset analyzed by 100 trials each time
  # this number is set to 100 instead of 1 in order to preserve processing power and make the simulation run in a reasonable time scale
  # if any if the optional stopping BF thresholds are reached set in the parameters stop_when_BF_low and stop_when_BF_high, the loop stops
  # computed bayes factors are collected in the vector BF_vect for visualization purposes, the counter vector ensures that the computed bayes factor values are stored consecutively in BF_vect  
  
  counter = 0
  BF_vect = NA
  for(i in seq(min_num_trials, max_num_trials, 100)){
    
    # dataset used in the bayes factor calculation, it includes data from the beggining the data_all_H1 full simulated dataset, the size of the smaple it contains increases in the loop until the BF threshold or the maxumum sample size is reached 
    data_H1 = data_all_H1[1:i]
    
    # probability table from the data
    prob_table <- table(data_H1)
    
    # Bayes Factor calculation from the proportion test
    # runs the test depending on the prior distribution to be used
    # if prior_distribution == "BUJ", use the "knowledge based" published  in Bem, D. J., Utts, J., & Johnson, W. O. (2011). Must psychologists change the way they analyze their data? Journal of Personality and Social Psychology, 101(4), 716-719.
    # if prior_distribution == "default" use the Cauchy prior from the Bayesfactor package with the width set in the parameter rho_scale_proptest
    # the BUJ prior is calculated from Bem's paper where the prior distribution is defined so that the mean is at 0 and the 90th percintele is at d = 0.5. Asuming a one-tailed normal distribution, this means that the sd is 0.5/1.28 = 0.3906. Source: Bem, D. J., Utts, J., & Johnson, W. O. (2011). Must psychologists change the way they analyze their data? Journal of Personality and Social Psychology, 101(4), 716-719.
    # the sd of the log odds ratio is computed from the sd of cohen's d shown above, source for the conversion formula is found in the compute.es R package https://cran.r-project.org/web/packages/compute.es/compute.es.pdf
    # BF01_norm runs the BF calculation with a proportion test using the normal distribution prior converted to log odds ratio (function by Richard Morey)
    # proportionBF is the standard function to do bayesian proportion test in the BayesFactor package, using the Cauchy prior
    
    if(prior_distribution == "BUJ"){
      p0 = H0_prob 
      shift = qlogis(p0)
      interval = alter_Interval_proptest
      r.scale = sqrt((pi^2*0.3906^2)/3)
      BF_proptest_H1_true <- BF01_norm(prob_table[2],length(data_H1),shift,r.scale,interval,p0) #nem kell a reciproka, a magasabb szam a 0-t tamogatja
    }
    if(prior_distribution == "default"){
      bf_proptest_rev = proportionBF(prob_table[1], length(data_H1), p = H0_prob, 
                                     rscale = rho_scale_proptest, nullInterval = alter_Interval_proptest)
      BF_proptest_H1_true <- as.numeric(matrix(1/bf_proptest_rev[2]))
    }
    
    counter = counter+1
    BF_vect[counter] = BF_proptest_H1_true
    if(!is.na(stop_when_BF_low)){if(stop_when_BF_low > BF_proptest_H1_true) break}
    if(!is.na(stop_when_BF_high)){if(stop_when_BF_high < BF_proptest_H1_true) break}
  }
  
  # plots the final BF curve for the simulation asuming H1 is true
  # this was used to help determin the optimal BF threshold, and continues to be included in the function to be able to track progress of the simulation 
  plot(log(BF_vect), type = "l", ylim = c(log(stop_when_BF_low), log(stop_when_BF_high)), yaxt = "n", xaxt = "n")
  abline(h = log(stop_when_BF_low))
  abline(h = log(stop_when_BF_high))
  
  
  
  ### Simulation assuming H0 is true  
  # the above data generation and analysis loop repeated with the assumption that H0 is true
  data_all_H0 <- rbinom(max_num_trials, size = 1, prob=H0_prob)
  
  for(i in seq(min_num_trials, max_num_trials, 100)){
    
    data_H0 = data_all_H0[1:i]
    
    prob_table <- table(data_H0)

    if(prior_distribution == "BUJ"){
      p0 = H0_prob
      shift = qlogis(p0)
      interval = alter_Interval_proptest
      r.scale = sqrt((pi^2*0.3906^2)/3)
      BF_proptest_H0_true <- BF01_norm(prob_table[2],length(data_H0),shift,r.scale,interval,p0) #nem kell a reciproka, a magasabb szam a 0-t tamogatja
    }
    if(prior_distribution == "default"){
      bf_proptest_rev = proportionBF(prob_table[1], length(data_H0), p = H0_prob, 
                                     rscale = rho_scale_proptest, nullInterval = alter_Interval_proptest)
      BF_proptest_H0_true <- as.numeric(matrix(1/bf_proptest_rev[2]))
    }
    
    if(!is.na(stop_when_BF_low)){if(stop_when_BF_low > BF_proptest_H0_true) break}
    if(!is.na(stop_when_BF_high)){if(stop_when_BF_high < BF_proptest_H0_true) break}
  }
  
  # output of the simulation function is the following: 
  # final BF from the dataset where H1 was true,
  # the number of trials where the simulated study ended in which H1 was true
  # the final BF from the dataset where H0 was true
  # the number of trials where the simulated study ended in which H0 was true
  output <- c(BF_proptest_H1_true, length(data_H1), BF_proptest_H0_true, length(data_H0))
  return(output)
}


######################################################################
#                                                                    #
#                     Running the simulaton                          #
#                                                                    #
######################################################################

# set the number of times the simulation needs to be run
# note that each simulation simulates two datasets and analyses, one where H1 is simulated to be true and one where H0 is simulated to be true
iterations = 5000 #5000 iterations run for about 4 hours i7 processor

# set the poparties of the plot area
# this is used mostly for tracking progress of this long simulation
# but the progress idnicator plots are the actual BF curves produced in each simulation, so it can also be used to fine tune BF thresholds for optional stopping and minimum and maximum number of trials 
par(mfrow=c(round(sqrt(iterations),0),round(sqrt(iterations),0)))
par(mar=c(0.1,0.1,0.1,0.1))


# set BF optional stopping thresholds which will be also used in the power analysis
BF_threshold_null = 60
BF_threshold = 1/BF_threshold_null

# runs the simulations iterations number of times
out <- replicate(iterations, simul(H0_prob = 0.5,
                                   H1_prob = 0.51,
                                   min_num_trials = 14000,
                                   max_num_trials = 60000,
                                   rho_scale_proptest = 1/2,
                                   one_sided = T,
                                   prior_distribution = "BUJ",
                                   stop_when_BF_low = BF_threshold,
                                   stop_when_BF_high = BF_threshold_null))

# output of the simulations are stored
sim_output <- as.data.frame(t(out))
names(sim_output) = c("BF_proptest_H1_true", "num_trials_H1_true", "BF_proptest_H0_true", "num_trials_H0_true")


# setting up a dataframe which will summarize the results of the power analysis simulation
power_results <- as.data.frame(matrix(NA, nrow = 1, ncol = 8))
names(power_results) = c("True_H1_when_H1_true", "Insensitive_H1_true", "False_H0_when_H1_true", "mean_num_trials_H1_true", "True_H0_when_H0_true", "Insensitive_H0_true", "False_H1_when_H0_true", "mean_num_trials_H0_true")


power_results[1,"True_H1_when_H1_true"] <- mean(sim_output[,"BF_proptest_H1_true"] < BF_threshold)
power_results[1,"False_H0_when_H1_true"] <- mean(sim_output[,"BF_proptest_H1_true"] > BF_threshold_null)
power_results[1,"Insensitive_H1_true"] <- 1-power_results[1,"True_H1_when_H1_true"]-power_results[1,"False_H0_when_H1_true"]
power_results[1,"mean_num_trials_H1_true"] <- mean(sim_output[,"num_trials_H1_true"])
power_results[1,"True_H0_when_H0_true"] <- mean(sim_output[,"BF_proptest_H0_true"] > BF_threshold_null)
power_results[1,"False_H1_when_H0_true"] <- mean(sim_output[,"BF_proptest_H0_true"] < BF_threshold)
power_results[1,"Insensitive_H0_true"] <- 1-power_results[1,"True_H0_when_H0_true"]-power_results[1,"False_H1_when_H0_true"]
power_results[1,"mean_num_trials_H0_true"] <- mean(sim_output[,"num_trials_H0_true"])

# prints out result summary of the power analysis
# indicates the poportions of alpha and beta errors for both when H1 and when H0 was simulated to be true, and the proportion of insensitive tests in which the macimum number of trials were reached without reaching the optional stopping BF thresholds
power_results

