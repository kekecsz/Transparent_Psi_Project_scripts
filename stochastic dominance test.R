
############################################
#              Set parameters              #
############################################
# number of participants
n = 800
# number of trials performed per participant
trial_size_per_participant = 20
# number of trials to simulate
max_num_trials = n * trial_size_per_participant
# proportion of hits in the total sample if alternative hypothesis is true
H1_prob = 0.51
# proportion of hits in the total sample if the null hypothesis is true
H0_prob = 0.5
# percentage of ESP-users, or ESP-capable individuals in the population 
ESP_virtuoso_percentage = 0.021 # this proportion needs to be more than twice as big as the effect size, so for H1_prob = 0.51, this needs to be higher than 0.02


############################################
#              Run simulation              #
############################################

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
data_part_H1 = split(data_part_H1_pre, ceiling(seq_along(data_part_H1_pre)/trial_size_per_participant))

# calculate proportion of correct guesses for each participant
success_proportions_part_H1 = sapply(data_part_H1, mean)

##### asuming that ESP ability is evenly distributed in the population

# simulate data with the assumption that the alternative hypothesis is true, but there are no individual differences in performance
# here we do this by simply shuffling randomly the simulated data where higher than chance hits were clustered within a few participants
# by the shuffling the higher than chance hit rate is now distributed evenly in the sample
data_all_H1_pre <- sample(data_part_H1_pre)
data_all_H1 = split(data_all_H1_pre, ceiling(seq_along(data_all_H1_pre)/trial_size_per_participant))

# calculate proportion of correct guesses for each participant
success_proportions_all_H1 = sapply(data_all_H1, mean)

############################################
#              Statistical test            #
############################################

##### if there are only a small number of ESP users
t.test(success_proportions_part_H1, mu = H0_prob, alternative = "greater")$p.value
prop.test(x = sum(data_part_H1_pre), n = max_num_trials, p = H0_prob, alternative = "greater")$p.value
# this is a K-S test to compare the sampling distribution under the H0 and the observer sampling distribution (which was simulated as if H1 was true)
# however the K-S test is not appropriate here because there are a a lot of "ties" (several people having the same observed proportion of successes) in the observed data
ks.test(success_proportions_part_H1, "pnorm", m = 0.5, sd = 0.112 , alternative = c("less"))$p.value

##### if everyone has the same ESP ability
t.test(success_proportions_all_H1, mu = H0_prob, alternative = "greater")$p.value
prop.test(x = sum(data_all_H1_pre), n = max_num_trials, p = H0_prob, alternative = "greater")$p.value
# this is a K-S test to compare the sampling distribution under the H0 and the observer sampling distribution (which was simulated as if H1 was true)
# however the K-S test is not appropriate here because there are a a lot of "ties" (several people having the same observed proportion of successes) in the observed data
ks.test(success_proportions_all_H1, "pnorm", m = 0.5, sd = 0.112 , alternative = c("less"))$p.value




