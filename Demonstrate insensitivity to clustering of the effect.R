# This scipt demonstrates that individual differences in ESP ability
# does not affect the proportion test approach we propose.
# Whether the ESP effects are clustered within a subgroup of the whole sample
# or the ESP ability is evenly distributed in the sample,
# the proportion test approach, where we pool all trials irrespective
# of participants, will produce the same statistical result.


############################################
#              Load packages               #
############################################
# BayesFactor package for the bayesian analysis 
library(BayesFactor)


############################################
#              Set parameters              #
############################################

# number of trials to simulate
max_num_trials = 16000
# number of trials performed per participant
trial_size_per_participant = 20
# proportion of hits in the total sample if alternative hypothesis is true
H1_prob = 0.51
# proportion of hits in the total sample if the null hypothesis is true
H0_prob = 0.5
# percentage of ESP-users, or ESP-capable individuals in the population 
ESP_virtuoso_percentage = 0.021 # this proportion needs to be more than twice as big as the effect size, so for H1_prob = 0.51, this needs to be higher than 0.02


############################################
#              Run simulation              #
############################################

##### no individual differences in performance

# simulate data with the assumption that the null hypothesis is true (no individual differences are possible in this scenario, 
# as each trial s just a random guess and everyone's guess is just as good as the others')
data_all_H0_pre <- rbinom(max_num_trials, size = 1, prob=H0_prob)

##### asuming that only a group of people are able to use ESP to predict future events

average_effect_size = H1_prob-H0_prob # overall effect size averaged over the whole sample
H1_part_prob = H0_prob + average_effect_size*(1/ESP_virtuoso_percentage) # calclulate effect size within the group of ESP users
ESP_virtuoso_trial_size = max_num_trials*ESP_virtuoso_percentage # number of trials (sample size) of ESP users

# simulate data for the ESP users
data_virt <- rbinom(ESP_virtuoso_trial_size, size = 1, prob=H1_part_prob)

# add ESP-user data to above simulated null effect data to get the total sample
data_part_H1_pre = c(data_all_H0_pre[1:(max_num_trials-ESP_virtuoso_trial_size)],data_virt)
data_part_H1 = split(data_part_H1_pre, ceiling(seq_along(data_part_H1_pre)/20))


# simulate data with the assumption that the alternative hypothesis is true, but there are no individual differences in performance
# here we do this by simply shuffling randomly the simulated data where higher than chance hits were clustered within a few participants
# by the shuffling the higher than chance hit rate is now distributed evenly in the sample
data_all_H1_pre <- sample(data_part_H1_pre)
data_all_H1 = split(data_all_H1_pre, ceiling(seq_along(data_all_H1_pre)/20))




############################################
#           Statistical analysis           #
############################################
# we average successes within each participant. This data will be used in the classical one-sample t-test approach
success_proportions_all_H1 = sapply(data_all_H1, mean)
success_proportions_part_H1 = sapply(data_part_H1, mean)

# run the t-test and extract its p-value
# everyone the same
t.test(success_proportions_all_H1, mu = H0_prob, alternative = "greater")$p.value
# only a handful of ESP-users
t.test(success_proportions_part_H1, mu = H0_prob, alternative = "greater")$p.value
# the t-test is affected by the grouping of the effect within a few individuals, because SD is higher in this case, so the t-test becomes less powerful
# although, the effect clustering only has a minute effect even on the t-test results if sample size is high enough.

# the proportion test
# everyone the same
prop.test(x = sum(unlist(data_all_H1)), n = length(unlist(data_all_H1)), p = H0_prob, alternative = "greater")$p.value
# only a handful of ESP-users
prop.test(x = sum(unlist(data_part_H1)), n = length(unlist(data_part_H1)), p = H0_prob, alternative = "greater")$p.value
# the proportion test's effectiveness is not hindered by grouping of the effect within a few individuals

# the Bayesian proportion test
# everyone the same
bf_proptest_rev = proportionBF(sum(unlist(data_all_H1)), length(unlist(data_all_H1)), p = H0_prob, 
                               rscale = 1/2, nullInterval = c(0.5,1))
BF_proptest <- as.numeric(matrix(1/bf_proptest_rev[1]))
BF_proptest # higher number supports H0

# only a handful of ESP-users
bf_proptest_rev = proportionBF(sum(unlist(data_part_H1)), length(unlist(data_part_H1)), p = H0_prob, 
                               rscale = 1/2, nullInterval = c(0.5,1))
BF_proptest <- as.numeric(matrix(1/bf_proptest_rev[1]))
BF_proptest # higher number supports H0
# the bayesian varient of the proportion test is also not affected, the same BF is returned for both randomly distributed and clustered effects.
