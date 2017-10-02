# This script demonstrates that if participants are instructed to stop when they
# have a few unsuccessful trials in the beggining, the t-test 
# will produce biassed results. On the other hand pooling all trial results
# irrespective of participants and doing a proportion test on this pooled
# results remains unbiased.

## The simulation runs for about 10 minutes with 1000 simulations on an i7 6600U 2.6 Ghz CPU

############################################
#             Load packages                #
############################################

library(progress) # this is necessary for the progress bar

############################################
#         Function for simulation          #
############################################

# function to simulate biased stopping
# this function generates a vector of random biomial numbers where 
# success (1) has a probability defined in "hit_chance".
# if we simulate H0 being true (no precognition), hit_chance is 0.5
# we also have to set the maximum number of trials a person is allowed
# to carry out, and the minimum number of trials the person will carry out, and the biased stopping rule
# stop_if_result_under defines the stopping rule of the participant. 
# If stop_if_result_under is set to 0.4, and min_trials_per_person is set to 3,
# it means that the participant will stop whenever his success rate across all
# of his own trials becomes lower than 0.4, but not until reaching 3 total completed trials

biased_sampler <- function(max_trials_per_person,
                           hit_chance,
                           min_trials_per_person,
                           stop_if_result_under){
  trial_results_vect = NA
  for(i in 1:max_trials_per_person){
    trial_results_vect[i] <- rbinom(1, size = 1, prob=hit_chance)
    if((length(trial_results_vect) > (min_trials_per_person-1)) & (mean(trial_results_vect) < stop_if_result_under)) break
  }
  return(trial_results_vect)
}


simul_stop_if_low = function(H0_prob, total_num_trials, max_trials_per_person,
                             hit_chance, min_trials_per_person,
                             stop_if_result_under){
  
  pb$tick() # adds a tick to the progress bar
  
  ############################################
  #              Run simulation              #
  ############################################
  
  # set the number of trials we will stop the experiment at 
  # (note that this is not the number of participants, which is irrelevant from our perspective,
  # rather, the total number of trials we will want to analyze)
  # we work with the number of trials because that is what we do in our experiment, because we cannot know 
  # how many people will stop pre-maturely we rather specify how many trials we will analyze in the end. 
  # This makes power analysis possible
  
  trial_number = total_num_trials
  
  # simulate the study
  # this loop generates data until we reach the total number of trials we pre-specified. 
  # If the last participant had more trials then we need, those last few data points are deleted so that we
  # will only have as many trials as we pre-specified
  trial_results_list = list(NA)
  for(i in 1:trial_number){
    trial_results_list[[i]] = biased_sampler(max_trials_per_person = max_trials_per_person,
                                             hit_chance = hit_chance,
                                             min_trials_per_person = min_trials_per_person, # to generate a more subtle bias where any individual will only be able to stop when he/she already completed a certain number of trials, increase this number
                                             stop_if_result_under = stop_if_result_under)
    
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
  H0_prob = H0_prob
  
  # we average successes within each participant. This data will be used in the classical one-sample t-test approach
  success_proportions = sapply(trial_results_list, mean)
  # run the t-test and extract its p-value and the mean of the proportion of successes
  t_p = t.test(success_proportions, mu = H0_prob)$p.value
  t_mean = mean(success_proportions)
  
  # run the proportion test (on the pooled data) and extract its p-value and the proportion of successes in the total dataset
  prop_p = prop.test(x = sum(unlist(trial_results_list)), n = length(unlist(trial_results_list)), p = H0_prob)$p.value
  prop_mean = sum(unlist(trial_results_list))/length(unlist(trial_results_list))
  
  return(c(t_p, t_mean, prop_p, prop_mean))
}



num_sim = 1000

# set up a progress bar
pb <- progress_bar$new(
  format = " simulation progress [:bar] :percent eta: :eta",
  total = num_sim, clear = FALSE, width= 60)

out = replicate(n = num_sim, simul_stop_if_low (H0_prob = 0.5,
                                   total_num_trials = 30060,
                                   max_trials_per_person = 18,
                                   hit_chance = 0.5,
                                   min_trials_per_person = 3,
                                   stop_if_result_under = 0.4))

output_frame = as.data.frame(t(out))
names(output_frame) = c("t_p", "t_mean", "prop_p", "prop_mean")

#set type one error probability
alpha = 0.005

# note that the proportion test is not fooled by the biased stopping rule of the participants
# because it works with the raw data disregarding participants.
# the percentage of simulations indicating significant effect is what we would expect
# based on the type one error rate we allo for
mean(output_frame[,"prop_p"] < alpha)
# and the raw data where we pool all trials is NOT biased
mean(output_frame[,"prop_mean"])

# on the other hand, the type-one error rate of the t-test is extremely high in this case
# in fact, the t-test will almost always be fooled as long as the number of trials we analyze reaches around 500
# this is because the t-test uses the data where we average succeses within participants, and 
# it is this average that is biased by the stopping rule of the participants
mean(output_frame[,"t_p"] < alpha)
mean(output_frame[,"t_mean"])


