
######################################################################
#                                                                    #
#                           Load packages                            #
#                                                                    #
######################################################################
library(progress)
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



# function to find mode
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


######################################################################
#                         Simulation function                        #
######################################################################

# function to simulate data and data analysis
# the function is used for power analysis

simul <- function(H0_prob, H1_prob, max_num_trials,
                  rho_scale_proptest, one_sided, prior_distribution, 
                  stop_when_BF_low, stop_when_BF_high, when_to_check_BF, 
                  NHST_or_BF, NHST_minimum_effect_threshold,
                  stop_when_p_H1, stop_when_p_H0, when_to_check_NHST,
                  robustness_test_with){
  
  pb$tick()
  
  # one or two sided tests
  # if one_sided == T a one-sided test will be performed 
  if(one_sided == T){
    proptest_alternative = "less" # alternative hypothesis is "less" because the default level is the first cell of the proportion table, which is the 0, in the alternative hypothesis we expect that there will be less 0-s (representing misses) then 1-s (representing hits) 
    alter_Interval_proptest = c(0.5,1)
  } else {
    proptest_alternative = "two.sided"
    alter_Interval_proptest = c(0,1)
  }
  
  # output of the simulation function is the following: 
  # final BF from the dataset where H1 was true,
  # the number of trials where the simulated study ended in which H1 was true
  # the final BF from the dataset where H0 was true
  # the number of trials where the simulated study ended in which H0 was true
  BF_proptest_H1_true = NA
  BF_proptest_H0_true = NA
  data_length_at_stop_BF_H1_true = NA
  data_length_at_stop_BF_H0_true = NA
  p_greater_then_nullprob_H1_true = NA
  p_less_than_minimum_effect_threshold_H1_true = NA
  p_greater_then_nullprob_H0_true = NA
  p_less_than_minimum_effect_threshold_H0_true = NA
  data_length_at_stop_NHST_H1_true = NA
  data_length_at_stop_NHST_H0_true = NA
  robustness_NHST_probtest_p_greater_then_nullprob_H1_true = NA
  robustness_NHST_probtest_p_less_than_minimum_effect_threshold_H1_true = NA
  robustness_NHST_probtest_p_greater_then_nullprob_H0_true = NA
  robustness_NHST_probtest_p_less_than_minimum_effect_threshold_H0_true = NA
  

  # generate data asuming H1 is true, the simulated data asumes each trial to be independent
  # see reasoning on this decision here: https://osf.io/rh5hb/wiki/Reasoning%20behind%20data%20analysis%20decisions/
  data_all_H1 <- rbinom(max_num_trials, size = 1, prob=H1_prob)
  
  # generate data asuming H0 is true
  data_all_H0 <- rbinom(max_num_trials, size = 1, prob=H0_prob)
  
  if(NHST_or_BF == "BF" | NHST_or_BF == "both"){
    ### Analysis on the dataset in which H1 is true
    # compute bayes factor on the simulated data where H1 is true
    # bayes factor is computed starting with a part of the dataset including the minimum sample size set in the function parameters and this is repeated in a loop increasing the dataset analyzed by 100 trials each time
    # this number is set to 100 instead of 1 in order to preserve processing power and make the simulation run in a reasonable time scale
    # if any if the optional stopping BF thresholds are reached set in the parameters stop_when_BF_low and stop_when_BF_high, the loop stops
    # computed bayes factors are collected in the vector BF_vect for visualization purposes, the counter vector ensures that the computed bayes factor values are stored consecutively in BF_vect  
    
    counter = 0
    BF_vect = NA
    for(i in when_to_check_BF){
      
      # dataset used in the bayes factor calculation, it includes data from the beggining the data_all_H1 full simulated dataset, the size of the smaple it contains increases in the loop until the BF threshold or the maxumum sample size is reached 
      data_H1_BF = data_all_H1[1:i]
      
      # probability table from the data
      prob_table <- table(data_H1_BF)
      
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
        BF_proptest_H1_true <- BF01_norm(prob_table[2],length(data_H1_BF),shift,r.scale,interval,p0) #nem kell a reciproka, a magasabb szam a 0-t tamogatja
      }
      if(prior_distribution == "default"){
        bf_proptest_rev = proportionBF(prob_table[1], length(data_H1_BF), p = H0_prob, 
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
#    plot(log(BF_vect), type = "l", ylim = c(log(stop_when_BF_low), log(stop_when_BF_high)), yaxt = "n", xaxt = "n")
#    abline(h = log(stop_when_BF_low))
#    abline(h = log(stop_when_BF_high))
    
    
    
    ### Analysis on the dataset in which H0 is true
    # the above analysis loop repeated on the dataset where H0 was simulated to be true
    for(i in when_to_check_BF){
      
      data_H0_BF = data_all_H0[1:i]
      
      prob_table <- table(data_H0_BF)
      
      if(prior_distribution == "BUJ"){
        p0 = H0_prob
        shift = qlogis(p0)
        interval = alter_Interval_proptest
        r.scale = sqrt((pi^2*0.3906^2)/3)
        BF_proptest_H0_true <- BF01_norm(prob_table[2],length(data_H0_BF),shift,r.scale,interval,p0) #nem kell a reciproka, a magasabb szam a 0-t tamogatja
      }
      if(prior_distribution == "default"){
        bf_proptest_rev = proportionBF(prob_table[1], length(data_H0_BF), p = H0_prob, 
                                       rscale = rho_scale_proptest, nullInterval = alter_Interval_proptest)
        BF_proptest_H0_true <- as.numeric(matrix(1/bf_proptest_rev[2]))
      }
      
      if(!is.na(stop_when_BF_low)){if(stop_when_BF_low > BF_proptest_H0_true) break}
      if(!is.na(stop_when_BF_high)){if(stop_when_BF_high < BF_proptest_H0_true) break}
    }
    
    data_length_at_stop_BF_H1_true = length(data_H1_BF)
    data_length_at_stop_BF_H0_true = length(data_H0_BF)
    
    
   
  }

  
  # robustness of BF results is tested with other statistical tests
  if(robustness_test_with == "NHST_probtest"){
    robustness_NHST_probtest_p_greater_then_nullprob_H1_true = prop.test(x = sum(data_H1_BF), n = length(data_H1_BF), p = H0_prob, alternative = "greater")$p.value
    robustness_NHST_probtest_p_less_than_minimum_effect_threshold_H1_true = prop.test(x = sum(data_H1_BF), n = length(data_H1_BF), p = NHST_minimum_effect_threshold, alternative = "less")$p.value
    robustness_NHST_probtest_p_greater_then_nullprob_H0_true = prop.test(x = sum(data_H0_BF), n = length(data_H0_BF), p = H0_prob, alternative = "greater")$p.value
    robustness_NHST_probtest_p_less_than_minimum_effect_threshold_H0_true = prop.test(x = sum(data_H0_BF), n = length(data_H0_BF), p = NHST_minimum_effect_threshold, alternative = "less")$p.value
  }
  
  
  
  
  
  
  # NHST test 
  if(NHST_or_BF == "NHST" | NHST_or_BF == "both"){
    
    # Analysis on the dataset where H1 is simulated to be true

    comparisons = 0
    for(i in when_to_check_NHST){
      comparisons = comparisons + 2 #because we do two tests at each stop point
      data_H1_NHST = data_all_H1[1:i]
      p_greater_then_nullprob_H1_true = prop.test(x = sum(data_H1_NHST), n = length(data_H1_NHST), p = H0_prob, alternative = "greater")$p.value
      p_greater_then_nullprob_H1_true = p.adjust(p = p_greater_then_nullprob_H1_true, method = "bonferroni", n = comparisons)
      p_less_than_minimum_effect_threshold_H1_true = prop.test(x = sum(data_H1_NHST), n = length(data_H1_NHST), p = NHST_minimum_effect_threshold, alternative = "less")$p.value
      p_less_than_minimum_effect_threshold_H1_true = p.adjust(p = p_less_than_minimum_effect_threshold_H1_true, method = "bonferroni", n = comparisons)
      
      if(!is.na(stop_when_BF_low)){if(stop_when_p_H1 > p_greater_then_nullprob_H1_true) break}
      if(!is.na(stop_when_BF_high)){if(stop_when_p_H0 > p_less_than_minimum_effect_threshold_H1_true) break}
    }
    
    # Analysis on the dataset where H0 is simulated to be true
    # Number of stops is counted in the vector "comparisons" for p value adjustment
    comparisons = 0
    for(i in when_to_check_NHST){
      comparisons = comparisons + 2 #because we do two tests at each stop point
      data_H0_NHST = data_all_H0[1:i]
      p_greater_then_nullprob_H0_true = prop.test(x = sum(data_H0_NHST), n = length(data_H0_NHST), p = H0_prob, alternative = "greater")$p.value
      p_greater_then_nullprob_H0_true = p.adjust(p = p_greater_then_nullprob_H0_true, method = "bonferroni", n = comparisons)
      p_less_than_minimum_effect_threshold_H0_true = prop.test(x = sum(data_H0_NHST), n = length(data_H0_NHST), p = NHST_minimum_effect_threshold, alternative = "less")$p.value
      p_less_than_minimum_effect_threshold_H0_true = p.adjust(p = p_less_than_minimum_effect_threshold_H0_true, method = "bonferroni", n = comparisons)
      
      if(!is.na(stop_when_BF_low)){if(stop_when_p_H1 > p_greater_then_nullprob_H0_true) break}
      if(!is.na(stop_when_BF_high)){if(stop_when_p_H0 > p_less_than_minimum_effect_threshold_H0_true) break}
    }
    
    data_length_at_stop_NHST_H1_true = length(data_H1_NHST)
    data_length_at_stop_NHST_H0_true = length(data_H0_NHST)
  }

  output <- c(  BF_proptest_H1_true,
                BF_proptest_H0_true,
                data_length_at_stop_BF_H1_true,
                data_length_at_stop_BF_H0_true,
                p_greater_then_nullprob_H1_true,
                p_less_than_minimum_effect_threshold_H1_true,
                p_greater_then_nullprob_H0_true,
                p_less_than_minimum_effect_threshold_H0_true,
                data_length_at_stop_NHST_H1_true,
                data_length_at_stop_NHST_H0_true,
                robustness_NHST_probtest_p_greater_then_nullprob_H1_true,
                robustness_NHST_probtest_p_less_than_minimum_effect_threshold_H1_true,
                robustness_NHST_probtest_p_greater_then_nullprob_H0_true,
                robustness_NHST_probtest_p_less_than_minimum_effect_threshold_H0_true)
  
  return(output)
}




######################################################################
#                                                                    #
#             Running the simulaton for multiple parameters          #
#                                                                    #
######################################################################

######################################################################
#                          Changable parameters                      #
######################################################################

# set the number of times the simulation needs to be run
# note that each simulation simulates two datasets and analyses, one where H1 is simulated to be true and one where H0 is simulated to be true
iterations = 100000 #10000 iterations with 3 stops per simulation runs for 8 minutes on i7 6600U 2.6GHz processor

# Determine whether to do classical NHST or BF analysis. This can be set to "BF", "NHST" or "both"
NHST_or_BF = "BF"

# if doing NHST as a main analysis what is the p-threshold of accepting the hypotheses (without correction for multiple comparisons)
# correction for multiple comparisons is built into the simulation function, it uses Bonferroni correction
p_threshold_to_try = c(0.001)

# If there is a bayes factor analysis, determine whether to run a NHST after stopping, to see if the results are robust to statistical approach
# THIS IS NOT refering to the NHST test as a main/primary hypothesis testing method. This is using the NHST to "confirm" the results of the BF analysis
# For testing H0 this analysis uses the NHST_minimum_effect_threshold set in the parameters of the simulation function below
robustness_test_with = c("NHST_probtest") # if set to NA no robustness test is conducted after data collection stopped, if set to "NHST_probtest", the robustness test will be a NHST proportion test
robustness_test_p_threshold_for_H1 = c(0.001)
robustness_test_p_threshold_for_H0 = c(0.001)

# Which BF thresholds to try in the simulation
BFs_to_try = c(70)

# Which prior distributuins to try
prior_distribution_to_try = c("BUJ") # BUJ or default. If prior distribution is default, it uses an "uninformed" Cauchy Prior. 
rho_to_try_if_default_prior = c(1/2) # 1/2 is the default setting for the width of the default Cauchy Prior in the BayesFactor package. 1/2 is less powerful at detecting H0 than BUJ prior, but more powerful at detecting H1 if effect size is p = 0.51. If rho is set to 1, "Ultrawide prior" produces similar results as BUJ in terms of power. 

# Which minimal, maximal, and middle stopping point sample sizes to try
min_trialnumber_to_try = c(16000)
max_trialnumber_to_try = c(80000)
mid_trialnumber_to_try = 40000 # if set to NA, a middle stopping point will be automatically assigned 1/3 of the way between min and max trialnum










######################################################################
#                   Simulation and output                            #
######################################################################

# create a dataframe which will house the final results of the simulation series
output_frame <- expand.grid(iterations, robustness_test_with, robustness_test_p_threshold_for_H1, robustness_test_p_threshold_for_H0, BFs_to_try, prior_distribution_to_try, rho_to_try_if_default_prior, min_trialnumber_to_try, mid_trialnumber_to_try, max_trialnumber_to_try, p_threshold_to_try)
names(output_frame) = c("Iterations", "Robustness_test_with", "robustness_test_p_threshold_for_H1", "robustness_test_p_threshold_for_H0", "BF_threshold", "Prior_distribution", "Prior_rho", "Min_trial_n", "Mid_trial_n", "Max_trial_n","p_threshold")
output_frame[,c("True_H1_when_H1_true_BF", "Insensitive_H1_true_BF", "False_H0_when_H1_true_BF", "True_H0_when_H0_true_BF", "Insensitive_H0_true_BF", "False_H1_when_H0_true_BF", "Mode_num_trials_H1_true_BF", "Mode_num_trials_H0_true_BF", 
                "True_H1_when_H1_true_NHST", "Insensitive_H1_true_NHST", "False_H0_when_H1_true_NHST", "True_H0_when_H0_true_NHST", "Insensitive_H0_true_NHST", "False_H1_when_H0_true_NHST", "Mode_num_trials_H1_true_NHST", "Mode_num_trials_H0_true_NHST",
                "ROBUST_True_H1_when_H1_true_BF", "ROBUST_Insensitive_H1_true_BF", "ROBUST_False_H0_when_H1_true_BF", "ROBUST_True_H0_when_H0_true_BF", "ROBUST_Insensitive_H0_true_BF", "ROBUST_False_H1_when_H0_true_BF",
                "Stopped_at_Min_trial_n_BF_H0_true", "Stopped_at_Mid_trial_n_BF_H0_true", "Stopped_at_Max_trial_n_BF_H0_true", "Stopped_at_Min_trial_n_BF_H1_true", "Stopped_at_Mid_trial_n_BF_H1_true", "Stopped_at_Max_trial_n_BF_H1_true", 
                "Stopped_at_Min_trial_n_NHST_H0_true", "Stopped_at_Mid_trial_n_NHST_H0_true", "Stopped_at_Max_trial_n_NHST_H0_true", "Stopped_at_Min_trial_n_NHST_H1_true", "Stopped_at_Mid_trial_n_NHST_H1_true", "Stopped_at_Max_trial_n_NHST_H1_true")] <- NA

# set up a progress bar
pb <- progress_bar$new(
  format = " simulation progress [:bar] :percent eta: :eta",
  total = iterations*nrow(output_frame), clear = FALSE, width= 60)

# This for loop runs the simulation series with all the possible combinations of the parameters to try given above
for(i in 1:nrow(output_frame)){
  
  # set BF optional stopping thresholds which will be also used in the power analysis
  BF_threshold_null = output_frame[i,"BF_threshold"]
  BF_threshold = 1/BF_threshold_null

  # set prior distribution and rho if applicable
  prior_distribution = output_frame[i,"Prior_distribution"]
  prior_rho = output_frame[i,"Prior_rho"]
  
  #set p_threhold for the main NHST analyses if any
  p_threshold_H1 = output_frame[i,"p_threshold"]
  p_threshold_H0 = output_frame[i,"p_threshold"]
  
  # set minimum and maximum sample size, and middle stopping point
  min_trialnumber = output_frame[i,"Min_trial_n"]
  max_trialnumber = output_frame[i,"Max_trial_n"]
  # set the midpoint stopping point
  if(is.na(mid_trialnumber_to_try)){mid_trialnumber <- min_trialnumber+((max_trialnumber - min_trialnumber)/3)} else {
    mid_trialnumber = mid_trialnumber_to_try
  }
  output_frame[i,"Mid_trial_n"] = mid_trialnumber 
     
      
  # runs the simulations iterations number of times
  out <- replicate(iterations, simul(H0_prob = 0.5,
                                     H1_prob = 0.51,
                                     max_num_trials = max_trialnumber,
                                     rho_scale_proptest = prior_rho,
                                     one_sided = T,
                                     prior_distribution = prior_distribution,
                                     stop_when_BF_low = BF_threshold,
                                     stop_when_BF_high = BF_threshold_null,
                                     when_to_check_BF = c(min_trialnumber, mid_trialnumber, max_trialnumber),
                                     NHST_or_BF = NHST_or_BF,
                                     NHST_minimum_effect_threshold = 0.51,
                                     stop_when_p_H1 = p_threshold_H1,
                                     stop_when_p_H0 = p_threshold_H0,
                                     when_to_check_NHST = c(min_trialnumber, mid_trialnumber, max_trialnumber),
                                     robustness_test_with = robustness_test_with))
      
  # output of the simulations are stored
  sim_output <- as.data.frame(t(out))
  names(sim_output) = c("BF_proptest_H1_true",
                        "BF_proptest_H0_true",
                        "data_length_at_stop_BF_H1_true",
                        "data_length_at_stop_BF_H0_true",
                        "p_greater_then_nullprob_H1_true",
                        "p_less_than_minimum_effect_threshold_H1_true",
                        "p_greater_then_nullprob_H0_true",
                        "p_less_than_minimum_effect_threshold_H0_true",
                        "data_length_at_stop_NHST_H1_true",
                        "data_length_at_stop_NHST_H0_true",
                        "robustness_NHST_probtest_p_greater_then_nullprob_H1_true",
                        "robustness_NHST_probtest_p_less_than_minimum_effect_threshold_H1_true",
                        "robustness_NHST_probtest_p_greater_then_nullprob_H0_true",
                        "robustness_NHST_probtest_p_less_than_minimum_effect_threshold_H0_true")
  
  # adds new columns to the sim_output frame to record whether main and robustness analysis decisions match
  sim_output[,"robustness_result_H1_true"] <- NA
  sim_output[,"robustness_result_H0_true"] <- NA

  # set p-values to 0.999 where both greater_then_nullprob and less_than_minimum_effect_threshold was true, this way these runs will be considered insensitive
  if(NHST_or_BF=="NHST" | NHST_or_BF=="both"){
    for(j in 1:nrow(sim_output)){
      if((sim_output[j,"p_greater_then_nullprob_H1_true"] < p_threshold_H1) & (sim_output[j,"p_less_than_minimum_effect_threshold_H1_true"] < p_threshold_H0)){
        sim_output[j,"p_greater_then_nullprob_H1_true"] = 0.999
        sim_output[j,"p_less_than_minimum_effect_threshold_H1_true"] = 0.999
      }
      if((sim_output[j,"p_greater_then_nullprob_H0_true"] < p_threshold_H1) & (sim_output[j,"p_less_than_minimum_effect_threshold_H0_true"] < p_threshold_H0)){
        sim_output[j,"p_greater_then_nullprob_H0_true"] = 0.999
        sim_output[j,"p_less_than_minimum_effect_threshold_H0_true"] = 0.999
      }
    }
  }
  
  # record whether main and robustness analysis decisions match
  if((NHST_or_BF=="BF" | NHST_or_BF=="both") & !is.na(output_frame[i,"Robustness_test_with"])){
    for(j in 1:nrow(sim_output)){
      # set p-values to 0.999 where both greater_then_nullprob and less_than_minimum_effect_threshold was true, this way these runs will be considered insensitive
      
      if((sim_output[j,"robustness_NHST_probtest_p_greater_then_nullprob_H1_true"] < p_threshold_H1) & (sim_output[j,"robustness_NHST_probtest_p_less_than_minimum_effect_threshold_H1_true"] < p_threshold_H0)){
        sim_output[j,"robustness_NHST_probtest_p_greater_then_nullprob_H1_true"] = 0.999
        sim_output[j,"robustness_NHST_probtest_p_less_than_minimum_effect_threshold_H1_true"] = 0.999
      }
      if((sim_output[j,"robustness_NHST_probtest_p_greater_then_nullprob_H0_true"] < p_threshold_H1) & (sim_output[j,"robustness_NHST_probtest_p_less_than_minimum_effect_threshold_H0_true"] < p_threshold_H0)){
        sim_output[j,"robustness_NHST_probtest_p_greater_then_nullprob_H0_true"] = 0.999
        sim_output[j,"robustness_NHST_probtest_p_less_than_minimum_effect_threshold_H0_true"] = 0.999
      }

      # determin outcome result if we take into account both main and robustness analyses
      if((sim_output[j,"BF_proptest_H1_true"] < BF_threshold) & (sim_output[j,"robustness_NHST_probtest_p_greater_then_nullprob_H1_true"] < output_frame[i,"robustness_test_p_threshold_for_H1"])) {sim_output[j,"robustness_result_H1_true"] = "True_H1"}
      if((sim_output[j,"BF_proptest_H1_true"] > BF_threshold_null) & (sim_output[j,"robustness_NHST_probtest_p_less_than_minimum_effect_threshold_H1_true"] < output_frame[i,"robustness_test_p_threshold_for_H0"])) {sim_output[j,"robustness_result_H1_true"] = "False_H0"}
      if(is.na(sim_output[j,"robustness_result_H1_true"])){sim_output[j,"robustness_result_H1_true"] = "Insensitive"}
         
      if((sim_output[j,"BF_proptest_H0_true"] > BF_threshold_null) & (sim_output[j,"robustness_NHST_probtest_p_less_than_minimum_effect_threshold_H0_true"] < output_frame[i,"robustness_test_p_threshold_for_H0"])) {sim_output[j,"robustness_result_H0_true"] = "True_H0"}
      if((sim_output[j,"BF_proptest_H0_true"] < BF_threshold) & (sim_output[j,"robustness_NHST_probtest_p_greater_then_nullprob_H0_true"] < output_frame[i,"robustness_test_p_threshold_for_H1"])) {sim_output[j,"robustness_result_H0_true"] = "False_H1"}
      if(is.na(sim_output[j,"robustness_result_H0_true"])){sim_output[j,"robustness_result_H0_true"] = "Insensitive"}
      
    }
  }
  


  # Fills up the final output dataframe with the results of the sequential simulations
  # Main BF analysis results if any
  if(NHST_or_BF=="BF" | NHST_or_BF=="both"){
    output_frame[i,"True_H1_when_H1_true_BF"] <- mean(sim_output[,"BF_proptest_H1_true"] < BF_threshold)
    output_frame[i,"False_H0_when_H1_true_BF"] <- mean(sim_output[,"BF_proptest_H1_true"] > BF_threshold_null)
    output_frame[i,"Insensitive_H1_true_BF"] <- 1-output_frame[i,"True_H1_when_H1_true_BF"]-output_frame[i,"False_H0_when_H1_true_BF"]
    output_frame[i,"Mode_num_trials_H1_true_BF"] <- Mode(sim_output[,"data_length_at_stop_BF_H1_true"])
    output_frame[i,"True_H0_when_H0_true_BF"] <- mean(sim_output[,"BF_proptest_H0_true"] > BF_threshold_null)
    output_frame[i,"False_H1_when_H0_true_BF"] <- mean(sim_output[,"BF_proptest_H0_true"] < BF_threshold)
    output_frame[i,"Insensitive_H0_true_BF"] <- 1-output_frame[i,"True_H0_when_H0_true_BF"]-output_frame[i,"False_H1_when_H0_true_BF"]
    output_frame[i,"Mode_num_trials_H0_true_BF"] <- Mode(sim_output[,"data_length_at_stop_BF_H0_true"])
    output_frame[i,c("Stopped_at_Min_trial_n_BF_H0_true", "Stopped_at_Mid_trial_n_BF_H0_true", "Stopped_at_Max_trial_n_BF_H0_true")] <- table(sim_output[,"data_length_at_stop_BF_H0_true"])/iterations
    output_frame[i,c("Stopped_at_Min_trial_n_BF_H1_true", "Stopped_at_Mid_trial_n_BF_H1_true", "Stopped_at_Max_trial_n_BF_H1_true")] <- table(sim_output[,"data_length_at_stop_BF_H1_true"])/iterations
  }
  # Main NHST analysis results if any
  if(NHST_or_BF=="NHST" | NHST_or_BF=="both"){
    output_frame[i,"True_H1_when_H1_true_NHST"] <- mean(sim_output[,"p_greater_then_nullprob_H1_true"] < p_threshold_H1)
    output_frame[i,"False_H0_when_H1_true_NHST"] <- mean(sim_output[,"p_less_than_minimum_effect_threshold_H1_true"] < p_threshold_H0)
    output_frame[i,"Insensitive_H1_true_NHST"] <- 1-output_frame[i,"True_H1_when_H1_true_NHST"]-output_frame[i,"False_H0_when_H1_true_NHST"]
    output_frame[i,"Mode_num_trials_H1_true_NHST"] <- Mode(sim_output[,"data_length_at_stop_NHST_H1_true"])
    output_frame[i,"True_H0_when_H0_true_NHST"] <- mean(sim_output[,"p_less_than_minimum_effect_threshold_H0_true"] < p_threshold_H0)
    output_frame[i,"False_H1_when_H0_true_NHST"] <- mean(sim_output[,"p_greater_then_nullprob_H0_true"] < p_threshold_H0)
    output_frame[i,"Insensitive_H0_true_NHST"] <- 1-output_frame[i,"True_H0_when_H0_true_NHST"]-output_frame[i,"False_H1_when_H0_true_NHST"]
    output_frame[i,"Mode_num_trials_H0_true_NHST"] <- Mode(sim_output[,"data_length_at_stop_NHST_H0_true"])
    output_frame[i,c("Stopped_at_Min_trial_n_NHST_H0_true", "Stopped_at_Mid_trial_n_NHST_H0_true", "Stopped_at_Max_trial_n_NHST_H0_true")] <- table(sim_output[,"data_length_at_stop_NHST_H0_true"])/iterations
    output_frame[i,c("Stopped_at_Min_trial_n_NHST_H1_true", "Stopped_at_Mid_trial_n_NHST_H1_true", "Stopped_at_Max_trial_n_NHST_H1_true")] <- table(sim_output[,"data_length_at_stop_NHST_H1_true"])/iterations
  }
  # Robustness analysis results if any
  if((NHST_or_BF=="BF" | NHST_or_BF=="both") & !is.na(output_frame[i,"Robustness_test_with"])){
    output_frame[i,"ROBUST_True_H1_when_H1_true_BF"] <- mean(sim_output[,"robustness_result_H1_true"] == "True_H1")
    output_frame[i,"ROBUST_False_H0_when_H1_true_BF"] <- mean(sim_output[,"robustness_result_H1_true"] == "False_H0")
    output_frame[i,"ROBUST_Insensitive_H1_true_BF"] <- mean(sim_output[,"robustness_result_H1_true"] == "Insensitive")

    output_frame[i,"ROBUST_True_H0_when_H0_true_BF"] <- mean(sim_output[,"robustness_result_H0_true"] == "True_H0")
    output_frame[i,"ROBUST_False_H1_when_H0_true_BF"] <- mean(sim_output[,"robustness_result_H0_true"] == "False_H1")
    output_frame[i,"ROBUST_Insensitive_H0_true_BF"] <- mean(sim_output[,"robustness_result_H0_true"] == "Insensitive")
  }
}

output_frame


# you may want to save results into a csv file so they are easier to look through
# setwd("C:\\Users\\zo0052ke\\Desktop\\Temp\\Simulation") # if so, set working directory
# write.csv(output_frame, file = "output_frame_csv.csv")



# just one way to visualize the results and compare different parameters
# plot_H1_Power_BF <- ggplot(data=output_frame, aes(x=BF_threshold, y=True_H1_when_H1_true_BF, group1=factor(Min_trial_n), group2=factor(Max_trial_n))) +
#              geom_line(aes(color=factor(Min_trial_n), linetype=factor(Max_trial_n)), size=1.2)+
#              geom_point(size=3)
#plot_H1_Power_BF


