
# This script contains the statistical analysis that will be used for the confirmatory
# analysis in the Transparent Psi Project
# Please read the Statistical analysis section to make sense of the code

######################################################################
#                                                                    #
#                           Load packages                            #
#                                                                    #
######################################################################
library(HDInterval) # needed to calcluate HDI credible intervals in the Bayesian parameter estimation robustness test
library(ggplot2) # for plotting
library(emdist) # to calcluate earth mover's distance (EMD)


######################################################################
#                                                                    #
#               Set parameters for example data                      #
#                                                                    #
######################################################################

### Set parameters of generating example dataset

# number of erotic trials performed per participant
trial_size_per_participant = 18

# probability of successful guess among non-ESP users 
# (this is also used as the null hypothesis in statistical tests)
M0_prob = 0.5

# ESP_user_percentage sets the percentage of ESP-users, or ESP-capable individuals 
# in the population 
# set this to 1 to simulate that everyone has the same level of ESP ability
# set this to 0 to simulate that noone has ESP ability
# set this to something in between to simulate that only a portion of the population can
# use ESP, and the others are just guessing randomly, 
# For example ESP_user_percentage = 0.03 and M1_prob = 0.65 simulates that only 3% 
# of the total poupulation can use ESP, and prediction success rate in this subgroup is 65%.
ESP_user_percentage = 1 

# probability of successful guess among ESP-users 
M1_prob = 0.51

# chance of stopping prematurely in each trial, to simulate some missing data due to unexpected events
# setting this to 0.002 would generate roughly 2% missing data
# note that premature stopping cannot introduce bias in our confirmatory analyses
chance_for_stopping_session = 0.002



######################################################################
#                                                                    #
#                            Functions                               #
#                                                                    #
######################################################################

######################################################################
#                  Bayes factor calculation functions                #
######################################################################


### Functions for Bayes factor caclulation using beta prior
# These functions are required to run the Bayes factor analysis 
# we thank Richard Morey for his help in developing these functions


fullAlt_beta = Vectorize(function(p, y, N, alpha, beta){
  exp(dbinom(y, N, p, log = TRUE) + dbeta(p, alpha, beta, log = TRUE)) 
},"p")

normalize_beta = function(alpha, beta, interval){
  diff(pbeta(interval, alpha, beta))
}

restrictedAlt_beta = function(p,y,N,y_prior,N_prior,interval){
  alpha = y_prior + 1
  beta = N_prior - y_prior + 1
  fullAlt_beta(p, y, N, alpha, beta) / normalize_beta(alpha, beta, interval) * (p>interval[1] & p<interval[2])
}

margLike_beta = function(y, N, y_prior, N_prior, interval){
  integrate(restrictedAlt_beta, interval[1], interval[2], 
            y = y, N = N, y_prior = y_prior, N_prior = N_prior, interval = interval)[[1]]
}

BF01_beta = Vectorize(function(y, N, y_prior, N_prior, interval, null_prob){
  dbinom(y, N, null_prob) / margLike_beta(y = y, N = N, y_prior = y_prior, N_prior = N_prior, interval = interval)
},"y")


######################################################################
#                       Other supporting functions                   #
######################################################################

### Function calculating the highest density interval using sampling
# We use hdi() from the library(HDInterval)
# this function is needed for the Bayesian parameter estimation robustness test

mode_HDI <- function(scale, density, crit_width = 0.95, n_samples = 1e5){
  samp <- sample(x = scale, size = n_samples, replace = TRUE, prob = density)
  hdi_result = hdi(samp, credMass=crit_width)
  result = c(scale[which(density == max(density))], # mode
             hdi_result[1], # lower bound
             hdi_result[2]) # upper bound
  
  # only needed for the names of the result
  Crit_lb = (1-crit_width)/2
  Crit_ub = crit_width + (1-crit_width)/2
  
  names(result) = c("mode", paste(Crit_lb*100, "%", sep = ""), paste(Crit_ub*100, "%", sep = ""))
  return(result)
}






######################################################################
#                                                                    #
#                       Generate example data                        #
#                                                                    #
######################################################################


# simulate the performance of 10,000 potential participants
# In this code we only simulate erotic trials
# Non-erotic trials will be ommitted from the hypothesis testing analysis 
# during data management.

list = list(NA)

for(i in 1:10000){
  ESP_user = rbinom(1, 1, prob = ESP_user_percentage)
  prob = M0_prob + (M1_prob - M0_prob)*ESP_user
  subj_data = as.data.frame(matrix(NA, nrow = 1, ncol = 2))
  
  for(j in 1:trial_size_per_participant){
    subj_data[j,1] = i
    subj_data[j,2] = rbinom(1, 1, prob = prob)
    list[[i]] = subj_data
    if(runif(1, min = 0, max = 1) < chance_for_stopping_session) break
  }
}

data = do.call("rbind", list)
names(data) = c("participant_ID", "success")


######################################################################
#                                                                    #
#                          Data analysis                             #
#                                                                    #
######################################################################

### analysis parameters
# these are the analysis parameters currently specified in our protocol

# interim analysis points (in total number of erotic trials performed)
when_to_check = c(30060, 37836, 45612, 53406, 61182, 68958, 76734, 84528, 92304, 100080)

# thresholds to infer support for M0 (high) or M1 (low)
Inference_threshold_BF_high = 25
Inference_threshold_BF_low = 1/Inference_threshold_BF_high

# this information is used both for calculating replication Bayes factor, and the Bayesian parameter estimation robustness test. 
# Here we use data from Bem's experiment 1, 828 successes within 1560 erotic trials
y_prior = 828 #number of successes in erotic trials in Bem's experiment 1
N_prior = 1560 # number of erotic trials in Bem's experiment 1


# smallest effect size of interest in the NHST equivalence test 
minimum_effect_threshold_NHST = 0.01
# p threshold for the NHST proportion test robustness test
Inference_threshold_robustness_NHST = 0.005

# in the Bayesian parameter estimation robustness test this will determine the region of practical 
#equivalence (ROPE) interval. The ROPE is interpreted similarly to SESOI, but not entireli the same. 
# See Kruschke, J. K., & Liddell, T. M. (2017). The Bayesian New Statistics: Hypothesis testing, 
# estimation, meta-analysis, and power analysis from a Bayesian perspective. 
# Psychonomic Bulletin & Review, 1-29. 
minimum_effect_threshold_Bayes_Par_Est = 0.006
# this threshold is used to set the HDI width to check against the ROPE in the Bayesian parameter 
# estimation robustness test, if ths parameter is set to 0.05 for example, it means that we would 
# expect that 95% of the probability mass would be within the ROPE to accept a hypothesis
Inference_threshold_robustness_Bayes_Par_Est = 0.05 


######################################################################
#                       Primary confirmatory test                    #
######################################################################

for(i in when_to_check){
  
  # start sampling from the top of the full simulated dataset (from the first trial of the first participant) 
  # until reaching the next interim analysis point
  data_BF = data[1:i,]
  
  # number of successes and total N of trials
  successes = sum(data_BF[,"success"])
  total_N = i
  
  #================================================================#
  #        Calculating Bayes factors using different priors        #
  #================================================================#
  
  ### Replication Bayes factor, with the Bem 2011 experiment 1 results providing the prior information
  
  BF_replication <- BF01_beta(y = successes, N = total_N, y_prior = y_prior, N_prior = N_prior, interval = c(0.5,1), null_prob = M0_prob) #numbers higher than 1 support the null
  
  
  ### Bayes factor with uniform prior
  # using a non-informative flat prior distribution with alpha = 1 and beta = 1
  
  BF_uniform <- BF01_beta(y = successes, N = total_N, y_prior = 0, N_prior = 0, interval = c(0.5,1), null_prob = M0_prob) #numbers higher than 1 support the null
  
  
  ### Bayes factor with BUJ prior
  # the BUJ prior is calculated from Bem's paper where the prior distribution is defined as a
  # normal distribution with a mean at 0 and 90th percentele is at medium effect size d = 0.5 
  # (we asume that this is one-tailed). Source: Bem, D. J., Utts, J., & Johnson, W. O. (2011). 
  # Must psychologists change the way they analyze their data? Journal of Personality and Social Psychology, 101(4), 716-719.
  # We simulate this in this binomial framework with a one-tailed beta distribution with alpha = 7 and beta = 7.
  # This distribution has 90% of its probability mass under p = 0.712, which we determined 
  # to be equivalent to d = 0.5 medium effect size. We used the formula to convert d to log odds ratio logodds = d*pi/sqrt(3), 
  # found here: Borenstein, M., Hedges, L. V., Higgins, J. P. T., & Rothstein, H. R. (2009). 
  # Converting Among Effect Sizes. In Introduction to Meta-Analysis (pp. 45-49): John Wiley & Sons, Ltd.
  # Then, log odds ratio vas converted to probability using the formula: p = exp(x)/(1+exp(x))
  # The final equation: exp(d*pi/sqrt(3))/(1+exp(d*pi/sqrt(3)))
  
  BF_BUJ <- BF01_beta(y = successes, N = total_N, y_prior = 6, N_prior = 12, interval = c(0.5,1), null_prob = M0_prob) #numbers higher than 1 support the null
  
  
  #================================================================#
  #                    Main analysis inference                     #
  #================================================================#
  
  # determine inference (supported model) based on the Bayes factors calculated above  
  if(Inference_threshold_BF_low >= max(c(BF_replication, BF_uniform, BF_BUJ))) {
    inference_BF = "M1"
    break} else if(Inference_threshold_BF_high <= min(c(BF_replication, BF_uniform, BF_BUJ))) {
      inference_BF = "M0"
      break} else {inference_BF = "Inconclusive"}
  
}





##################################################
#                 Robustness tests               #
##################################################

#================================================================#
#               Robustness test of BF results with NHST          #
#================================================================#

# robustness of BF results is tested with NHST proportion tests
# here we perform both an equality test and an equivalence test to draw statistical inference


# equality test
equality_test_p = prop.test(x = successes, n = total_N, p = M0_prob, alternative = "greater")$p.value

# equivalence test
equivalence_test_p = prop.test(x = successes, n = total_N,
                                                                          p = M0_prob+minimum_effect_threshold_NHST, 
                                                                          alternative = "less")$p.value

# descritives of the frequentist estimate
pbar = successes/total_N
SE = sqrt(pbar * (1-pbar)/total_N)
E = qnorm(.9975) * SE
proportion_995CI = round(pbar + c(-E , E), 3)



# making inference decision
if((Inference_threshold_robustness_NHST > equality_test_p) & 
   (Inference_threshold_robustness_NHST <= equivalence_test_p)){
  inference_robustness_NHST = "M1"} else if((Inference_threshold_robustness_NHST > equivalence_test_p) & 
                                            (Inference_threshold_robustness_NHST <= equality_test_p)){
    inference_robustness_NHST = "M0" 
  } else {inference_robustness_NHST = "Inconclusive"}


#=======================================================================#
#   Robustness test of BF results with Bayesian parameter estimation    #
#=======================================================================#

# robustness of BF results is tested by calculating HDI of the posteriro distribution and checking its relation to
# the region of practical equivalence (ROPE), promoted in Kruschke, J. K., & Liddell, T. M. 
# (2017). The Bayesian New Statistics: Hypothesis testing, estimation, meta-analysis, and power 
# analysis from a Bayesian perspective. Psychonomic Bulletin & Review, 1-29. 

# calculate posterior distribution using beta distribution
prior_alpha = y_prior + 1
prior_beta = N_prior-y_prior+1

posterior_alpha = prior_alpha + successes
posterior_beta = prior_beta + total_N - successes

scale = seq(0, 1, length = 1001)
posterior_density = dbeta(scale, posterior_alpha, posterior_beta)

# calculate HDI for the posterior distribution
# (here we calculate the upper and lower bound of the 90% of the probability mass
# because we use a one-tailed test. This means that the total probability mass below
# the upper bound of the 90% HDI will be 95%)
hdi_result = mode_HDI(scale = scale, density = posterior_density, crit_width = 1-Inference_threshold_robustness_Bayes_Par_Est*2, n_samples = 1e6)

# parameters for decision making
HDI_lb = hdi_result[2]
HDI_ub = hdi_result[3]

# making inference decision
ROPE = M0_prob+minimum_effect_threshold_Bayes_Par_Est
if(HDI_lb >= ROPE){inference_robustness_Bayes_Par_Est = "M1"
} else if(HDI_ub <= ROPE){inference_robustness_Bayes_Par_Est = "M0"
} else {inference_robustness_Bayes_Par_Est = "Inconclusive"}


#=======================================================================#
#      Determine final inference of all robustness tests combined       #
#=======================================================================#

# the main analysis inference is only robust if all robustness tests came to the same inference as the main test
inferences = c(inference_robustness_NHST, inference_robustness_Bayes_Par_Est)
inference_robustness = if(all(inferences == inferences[1])){inferences[1]} else {"mixed"}

Robust = if(inference_BF == inference_robustness){"robust"} else {"not robust"}




##################################################
#       Report of confirmatory analysis          #
##################################################

# number of erotic trials performed (without the last participant, whose data might
# be cut in half because of the interim stopping point)
numtrials_without_last_part = nrow(data_BF[data_BF[,"participant_ID"] != max(data_BF[,"participant_ID"]),])
# number of missing trials/data points
missing_data_points = (length(unique(data_BF[,"participant_ID"])) * trial_size_per_participant) - numtrials_without_last_part 


general_text = paste("The study stopped at ", total_N, " erotic trials gathered from a total of ",
      max(data_BF[,"participant_ID"]), " participants.", " There has been ", missing_data_points, " (", round(missing_data_points/total_N*100,2),"%) missing data points due to unfinished sessions.", 
      " We observed a total of ", round(mean(data_BF[,"success"]), 4)*100, "% successful guesses within ",
      total_N, " erotic trials (99.5% CI = ", proportion_995CI[1]*100, "%", ", ", proportion_995CI[2]*100, "%", 
      "; posterior mode = ", hdi_result[1]*100, "%", ", posterior 90% HDI = ", HDI_lb*100, "%", ", ",
      HDI_ub*100, "%", ").", sep = "")

robustness_text = if(Robust == "robust"){paste(" The results proved to be robust to different statistical approaches, increasing our confidence in our inference.")} else if(Robust == "not robust"){
  paste(" However, the results did not prove to be robust to different statistical approaches.")}


if(inference_BF == "M1"){
  paste(general_text, " Observing this success rate is ", round(1/max(c(BF_replication, BF_uniform, BF_BUJ)),0),
        " times more likely if humans can guess future randomly determined events than if they are guessing randomly. Taken at face value, the data provide strong evidence that the probability of successfully guessing later computer-generated random events is higher than chance level as previously reported by Bem (2011) and others (Bem, Tressoldi, Rabeyron, & Duggan, 2015).",
        robustness_text, sep = "")
} else if(inference_BF == "M0"){
  paste(general_text, " Observing this success rate is ", round(min(c(BF_replication, BF_uniform, BF_BUJ)),0),
        " times more likely if humans are guessing randomly than if they can guess future randomly determined events. Taken at face value, the data provide strong evidence that the probability of successfully guessing later computer-generated random events is not higher than chance level as previously reported by Bem (2011) and others (Bem, Tressoldi, Rabeyron, & Duggan, 2015).",
        robustness_text, sep = "")
  } else if(inference_BF == "Inconclusive" & max(c(BF_replication, BF_uniform, BF_BUJ)) < 1){
  paste(general_text, " Observing this success rate is ", round(1/max(c(BF_replication, BF_uniform, BF_BUJ)),0),
        " times more likely if humans can guess future randomly determined events than if they are guessing randomly. However, this study outcome did not reach the pre-specified criteria of strong support for either model.",
        robustness_text, sep = "")
} else if(inference_BF == "Inconclusive" & min(c(BF_replication, BF_uniform, BF_BUJ)) > 1){
  paste(general_text, " Observing this success rate is ", round(min(c(BF_replication, BF_uniform, BF_BUJ)),0),
        " times more likely if humans are guessing randomly than if they can guess future randomly determined events. However, this study outcome did not reach the pre-specified criteria of strong support for either model.",
        robustness_text, sep = "")
} else {paste(general_text, "However, this study outcome did not reach the pre-specified criteria of strong support for either model.",
              robustness_text, sep = "")}











######################################################################
#                         Exploratory analysis                       #
######################################################################
# The existence of individual differences between participants will be evaluated 
# using two statistical approaches, one using a frequentist approach, the other Bayesian 
# inference. In both approaches, we assume that it is possible that only a small 
# subgroup of the population we sample from is capable of predicting future events
# with better than chance accuracy.


# Set parameters ESP_user_percentage and M1_prob to simulate a samll subgroup of
# ESP users. For example set ESP_user_percentage = 0.03 and M1_prob = 0.65 to simulate
# that only 3% of the total poupulation can use ESP, and prediction success rate
# in this subgroup is 65%.

# EXPLORATORY ANALYSIS RESULTS WILL NOT AFFECT THE CONCLUSIONS OF OUR STUDY

#=======================================================================#
#           Comparison of expected and observed distributions           #
#=======================================================================#

# calculate proportion of successful guesses for each participant in the observed data
data_BF_split = split(data_BF, f = data_BF[,"participant_ID"])
success_proportions_empirical = sapply(data_BF_split, function(x) mean(x[,"success"]))

# samples 1,000,000 participants from a population with H0 success rate
# this is used for the stochastic dominance test as the null model
# we call this the theoretical sample, because it approximates the theoretical null model
sim_null_participant_num = 1000000
success_proportions_theoretical <- rbinom(sim_null_participant_num, size = trial_size_per_participant, prob=M0_prob)/trial_size_per_participant

# determine possible values of success rates
possible_success_rates = 0
for(i in 1:trial_size_per_participant){
  possible_success_rates[i+1] = round(1/(trial_size_per_participant/i), 2)
}
possible_success_rates_char = as.character(possible_success_rates)
success_proportions_theoretical_char_rounded = as.character(round(success_proportions_theoretical, 2))
success_proportions_empirical_char_rounded = as.character(round(success_proportions_empirical, 2))

success_rates_theoretical = NA
for(i in 1:length(possible_success_rates)){
  success_rates_theoretical[i] = sum(success_proportions_theoretical_char_rounded == possible_success_rates_char[i])
}
success_rates_theoretical_prop = matrix(success_rates_theoretical/sum(success_rates_theoretical))


success_rates_empirical = NA
for(i in 1:length(possible_success_rates)){
  success_rates_empirical[i] = sum(success_proportions_empirical_char_rounded == possible_success_rates_char[i])
}
success_rates_empirical_prop = matrix(success_rates_empirical/sum(success_rates_empirical))



# plot 
histogram_plot_data = as.data.frame(c(success_rates_theoretical_prop, success_rates_empirical_prop))
histogram_plot_data = cbind(histogram_plot_data, factor(c(rep("Expected if M0 is true", length(success_rates_theoretical_prop)), rep("Observed", length(success_rates_empirical_prop)))))
histogram_plot_data = cbind(histogram_plot_data, factor(rep(possible_success_rates_char, 2)))
names(histogram_plot_data) = c("proportion", "group", "success")

figure_1 =  ggplot(histogram_plot_data, aes(y = proportion, x = success, group = group))+
  geom_bar(aes(fill = group), alpha = 0.5, stat = "identity", position = "identity")+
  scale_fill_manual(values = c("darkgrey", "black")) +
  xlab("Successful guess rate") +
  ylab("Proportion") +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        axis.line = element_line(colour = "black", size = 1.2),
        axis.text.x = element_text(angle = 90, color = "black", face = "bold", size = 12, margin = margin(1,0,0,0,"mm"), vjust = 0.5),
        axis.text.y = element_text(face = "bold", color = "black", size = 12),
        axis.title = element_text(size = 16))
figure_1

# earth mover's distance
emd2d(success_rates_theoretical_prop,success_rates_empirical_prop)



#=======================================================================#
#                             Bayesian aproach                          #
#=======================================================================#

# The computer code of the Bayesian exploratory analysis is still under construction

