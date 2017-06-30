

############################################
#         Function for simulation          #
############################################

# function to simulate biased stopping
# this function generates a vector of random biomial numbers where success (1) has a probability defined in "hit_chance".
# if we simulate H0 being true (no precognition), hit_chance is 0.5
# we also have to set the maximum number of trials a person is allowed to carry out, and the minimum number of trials the person will carry out, and the biased stopping rule
# stop_if_result_over defines the stopping rule of the participant. 
# If stop_if_result_over is set to 0.5, it means that the participant will stop whenever his success rate across all of his own trials is becomes greater than 0.5

biased_sampler <- function(max_trials_per_person,
                           hit_chance,
                           min_trials_per_person,
                           stop_if_result_over){
  trial_results_vect = NA
  for(i in 1:max_trials_per_person){
    trial_results_vect[i] <- rbinom(1, size = 1, prob=hit_chance)
    if((length(trial_results_vect) > (min_trials_per_person-1)) & (mean(trial_results_vect) > stop_if_result_over)) break
  }
  return(trial_results_vect)
}



############################################
#              Run simulation              #
############################################

# set the number of trials we will stop the experiment at 
# (note that this is not the number of participants, which is irrelevant from our perspective,
# rather, the total number of trials we will want to analyze)
# we work with the number of trials because that is what we do in our experiment, because we cannot know 
# how many people will stop pre-maturely we rather specify how many trials we will analyze in the end. 
# This makes power analysis possible

trial_number = 1000

# simulate the study
# this loops generates data until we reach the total number of trials we pre-specified. 
# If the last participant had more trials then we need, those last few data points are deleted to so that we
# will only have as many trials as we pre-specified
trial_results_list = list(NA)
for(i in 1:trial_number){
  trial_results_list[[i]] = biased_sampler(max_trials_per_person = 20,
                                         hit_chance = 0.5,
                                         min_trials_per_person = 1, # to generate a more subtle bias where any individual will only be able to stop when he/she already completed a certain number of trials, increase this number
                                         stop_if_result_over = 0.5)

  # these ifs make sure that we stop when we reached the pre-specified number of trials and that we don't have more data than we pre-specified
  # the demonstration also works if we specify number of participants instead of number of trials
  # but it resembles our study design less
  # to see this, these two ifs needs to be deleted
  # in that case trial_number will mean number of participants
    if(length(unlist(trial_results_list)) > trial_number){
    trial_results_list[[i]] = trial_results_list[[i]][1:(length(trial_results_list[[i]])-(length(unlist(trial_results_list)) - trial_number))]
  }
  if(length(unlist(trial_results_list)) == trial_number) break
}

############################################
#           Statistical analysis           #
############################################

# set the null hypothesis success probability, which will be used in the statistical tests
H0_prob = 0.5

# we average successes within each participant. This data will be used in the classical one-sample t-test approach
success_proportions = sapply(trial_results_list, mean)
# run the t-test and extract its p-value
# note that with this extreme stopping rule by the participants
# (stop whenever overall hit rate is abolve p = 0.5 regardless of how many trials the participant did)
# the t-test will always be fooled (as long as the number of trials we analyze reaches around 100)
# this is solely because of the averaging of succeses within participants, which is the data used by the t-test
# so it is not the t-test which is biased, it is the data entered into the t-test that is already biased
t.test(success_proportions, mu = H0_prob, alternative = "greater")$p.value
# this is shown in the mean of the averaged success probabilities
mean(success_proportions)


# by contrast, the proportion test is not fooled by the biased stopping rule of the participants
# because it works with the raw data disregarding participants.
prop.test(x = sum(unlist(trial_results_list)), n = length(unlist(trial_results_list)), p = H0_prob, alternative = "greater")$p.value
# the raw data where we pool all trials is NOT biased
sum(unlist(trial_results_list))/length(unlist(trial_results_list))



