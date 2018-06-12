######################################################################
#                                                                    #
#                           Load packages                            #
#                                                                    #
######################################################################

library(HDInterval) # needed to calcluate HDI credible intervals in the Bayesian parameter estimation robustness test
library(ggplot2) # for visualization


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
#                         Data management                            #
#                                                                    #
######################################################################

# get up to date date from github
pilot_data_pre <- read.csv("https://raw.githubusercontent.com/gy0p4k/transparent-psi-results/master/results.csv")

# only read pilot sessions
# pilot sessions are marked with session_type = "pilot" and they were all collected in the lab with the laboratory_ID_code = "lab_ELTE_01"
lab_ID <- "lab_ELTE_01"
pilot_data <- pilot_data_pre[pilot_data_pre[,"session_type"] == "pilot" & pilot_data_pre[,"laboratory_ID_code"] == "lab_ELTE_01", ]


# Number of participants tested in the pilot test
sample_size_participants = length(unique(pilot_data[, "participant_ID"]))


pilot_data[, "trial_number"] = as.numeric(pilot_data[, "trial_number"])

data_BF = pilot_data[!is.na(pilot_data[, "trial_number"]), ]




######################################################################
#                                                                    #
#                          Data analysis                             #
#                                                                    #
######################################################################

### analysis parameters
# these are the analysis parameters currently specified in our protocol

# probability of success if M0 is true
M0_prob = 0.5

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

#================================================================#
#                Calculate number of successes                   #
#================================================================#

# number of successes and total N of trials
sides_match = data_BF[,"guessed_side"] == data_BF[,"target_side"]
successes = sum(sides_match)
total_N = nrow(data_BF)


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














######################################################################
#                                                                    #
#                   Interpretation and visualization                 #
#                                                                    #
######################################################################


# visualize results
# ALL THREE BF values need to cross the threshold for the interpretation to be accepted

BF_results_for_plotting = cbind(as.data.frame(c(BF_replication, BF_uniform, BF_BUJ)), c("BF_replication", "BF_uniform", "BF_BUJ"))
names(BF_results_for_plotting) = c("Bayes_factor_01", "BF_type")

plot <- ggplot(BF_results_for_plotting, aes(y = Bayes_factor_01, x = BF_type))+
  geom_point()+
  geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=c(Inference_threshold_BF_high), ymax=c(Inf)), alpha = 0.2, fill=c("pink"))+
  geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=c(Inference_threshold_BF_low), ymax=c(Inference_threshold_BF_high)), alpha = 0.2, fill=c("grey80"))+
  geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=c(0), ymax=c(Inference_threshold_BF_low)), alpha = 0.2, fill=c("lightgreen"))+
  geom_point(size = 3.5, shape = 21, fill = "white")+
  scale_y_log10(limits = c(0.005,200), breaks=c(0.01, Inference_threshold_BF_low, 0.1, 0.33, 0, 3, 10, Inference_threshold_BF_high, 100))+
  geom_hline(yintercept = c(Inference_threshold_BF_low, Inference_threshold_BF_high), linetype = "dashed")+
  geom_text(aes(x=0.5, y=c(100, 1, 0.01), label=c("Supports M0", "Inconclusive", "Supports M1"), angle = 270))

# there is a warning because of the log scale. We suppres that and print the plot.
suppressWarnings(print(plot))


# plot results of the Bayesian parameter estimation used as a robustness test
plot(scale, posterior_density, type="l", lty=1, xlab="x value", xlim = c(0.45, 0.55),
     ylab="Density", main="Bayesian parameter esitmation updating prior data")
abline(v=c(M0_prob, ROPE), lty = c(1, 2))

density_table = as.data.frame(cbind(scale, posterior_density))

height_lim_lb = density_table[density_table[, "scale"] == HDI_lb, "posterior_density"]
height_lim_ub = density_table[density_table[, "scale"] == HDI_ub, "posterior_density"]
clip(0,1,-10,height_lim_lb)
abline(v=HDI_lb, lty = 3)
clip(0,1,-10,height_lim_ub)
abline(v=HDI_ub, lty = 3)





# Text summary of the results  
print(paste("So far ", sample_size_participants, " participants took part in the study providing data for a total of ", total_N, " trials. The Bayes factors are BF_replication = ", round(BF_replication, 2), ", BF_uniform = ", round(BF_uniform, 2), ", and BF_BUJ = ", round(BF_BUJ, 2), ". If the study would stop now, the interpretation would be the following: '", inference_BF, "'.", " According to the robustness analysis, the main analysis result is ", Robust, ".", sep = ""))

