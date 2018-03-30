##########################################
#          set working directory         #
##########################################
file_location <- "C:\\Users\\zo0052ke\\Dropbox\\Bayesian Statistics\\STAN\\Transparent Psi Project in rstan\\"
# file_location <- "C:\\Users\\kekec\\Dropbox\\Bayesian Statistics\\STAN\\Transparent Psi Project in rstan\\"

setwd(file_location)

##########################################
#              Load packages             #
##########################################
library(rstan)


############################################
#              Custom functions            #
############################################

# function for rounding to a specific base number
mround <- function(x,base){ 
  base*round(x/base) 
} 


############################################
#   Set parameters for demonstration data  #
############################################

# number of trials performed per participant
trial_size_per_participant = 20
# number of participants to simulate
participant_num = 2000
# proportion of hits in the total sample if alternative hypothesis is true
H1_prob = 0.55
# proportion of hits in the total sample if the null hypothesis is true
H0_prob = 0.50
# percentage of ESP-users, or ESP-capable individuals in the population 
ESP_virtuoso_percentage = 0.20 # this proportion needs to be more than twice as big as the effect size, so for H1_prob = 0.51, this needs to be higher than 0.02

# if exact is T it will generate exactly the above set proportion of successes
# if exact is F it will randomly generate the data from bernulli trials with the above set parameters
exact = T

##########################################
#              Generate data             #
##########################################

max_num_trials = participant_num*trial_size_per_participant # number of trials to simulate

average_effect_size = H1_prob-H0_prob # overall effect size averaged over the whole sample

H1_part_prob = H0_prob + average_effect_size*(1/ESP_virtuoso_percentage) # calclulate effect size within the group of ESP users
ESP_virtuoso_trial_size = mround(max_num_trials*ESP_virtuoso_percentage, trial_size_per_participant) # number of trials (sample size) of ESP users

##### no individual differences in performance

# simulate data with the assumption that the null hypothesis is true (no individual differences are possible in this scenario, 
# as each trial s just a random guess and everyone's guess is just as good as the others')
if(exact == T){
  data_all_H0_pre <- sample(c(rep(1, round(max_num_trials*H0_prob, 0)), rep(0, round(max_num_trials*(1-H0_prob), 0))))
  # simulate data for the ESP users
  data_virt <- sample(c(rep(1, round(ESP_virtuoso_trial_size*H1_part_prob, 0)), rep(0, round(ESP_virtuoso_trial_size*(1-H1_part_prob), 0))))
  data_null_to_add <- sample(c(rep(1, round((max_num_trials-ESP_virtuoso_trial_size)*H0_prob, 0)), rep(0, round((max_num_trials-ESP_virtuoso_trial_size)*(1-H0_prob), 0))))
} else {
  data_all_H0_pre <- rbinom(max_num_trials, size = 1, prob=H0_prob)
  # simulate data for the ESP users
  data_virt <- rbinom(ESP_virtuoso_trial_size, size = 1, prob=H1_part_prob)
  data_null_to_add <- data_all_H0_pre[1:(max_num_trials-ESP_virtuoso_trial_size)]
}

data_all_H0_pre2 <- as.data.frame(data_all_H0_pre)
data_all_H0 <- cbind(rep(1:(max_num_trials/trial_size_per_participant), each = trial_size_per_participant), data_all_H0_pre2)
names(data_all_H0) = c("participant_ID", "response")


# add ESP-user data to above simulated null effect data to get the total sample
data_part_H1_pre = as.data.frame(c(data_null_to_add,data_virt))
random_participant_ID = sample(1:(max_num_trials/trial_size_per_participant))
data_part_H1 = cbind(rep(random_participant_ID, each = trial_size_per_participant), data_part_H1_pre)
names(data_part_H1) = c("participant_ID", "response")
data_part_H1 = data_part_H1[order(data_part_H1[,"participant_ID"]),]



# get data with the assumption that the alternative hypothesis is true, but there are no individual differences in performance
# here we do this by simply shuffling randomly the simulated data where higher than chance hits were clustered within a few participants
# by the shuffling the higher than chance hit rate is now distributed evenly in the sample
data_all_H1_pre <- as.data.frame(sample(data_part_H1[,"response"]))
data_all_H1 <- cbind(rep(1:(max_num_trials/trial_size_per_participant), each = trial_size_per_participant), data_all_H1_pre)
names(data_all_H1) = c("participant_ID", "response")




# check that the effect is roughly what we want it to be
nonESP_participants = random_participant_ID[1:((max_num_trials-ESP_virtuoso_trial_size)/trial_size_per_participant)]
ESPvirt_participants = random_participant_ID[(((max_num_trials-ESP_virtuoso_trial_size)/trial_size_per_participant)+1):length(random_participant_ID)]

print(paste("ESP user success rate = ", round(mean(data_part_H1[data_part_H1[,"participant_ID"] %in% ESPvirt_participants,"response"]), 3)))
print(paste("ESP non-user success rate = ", round(mean(data_part_H1[data_part_H1[,"participant_ID"] %in% nonESP_participants,"response"]), 3)))
print(paste("total success rate in the clustered data = ", round(mean(data_part_H1[,"response"]), 3)))

print(paste("total success rate in the unclustered non-null data = ", round(mean(data_all_H1[,"response"]), 3)))

print(paste("total success rate in the unclustered null data = ", round(mean(data_all_H0[,"response"]), 3)))



##### write generated data files for reproducability
# write.csv(data_part_H1, file = paste(file_location, "data_part_H1_virt_perc_", ESP_virtuoso_percentage, "_virt_effect_", round(H1_part_prob, 3), ".csv", sep = ""), row.names = F)
# write.csv(data_all_H1, file = paste(file_location, "data_all_H1_total_effect_", round(H1_prob, 3), ".csv", sep = ""), row.names = F)
# write.csv(data_all_H0, file = paste(file_location, "data_all_H0_total_effect_", round(H0_prob, 3), ".csv", sep = ""), row.names = F)



##### Or load previously simulated data from file
# In this dataset 75% of the participants have no ESP ability and are just guessing randomly (responses generated with rbinom with p = 0.5), 
# and 25% of the participants have a small ESP ability (responses generated with rbinom with p = 0.54). 
# It has 700 simulated participants all of whom completed 20 trials.
# rsimdata <- read.csv("simdata.csv",header=TRUE)

##########################################
#             Data management            #
##########################################

# create data as list for Stan:
# three datasets are used:

# 1. where a small subset of participants are simulated to have ESP
stansimdata_part_H1 <- list(response = data_part_H1$response, subj = data_part_H1$participant_ID, N = nrow(data_part_H1), J = length(unique(data_part_H1$participant_ID)))
# 2. where ESP ability is the same for every individual (no personal differences)
stansimdata_all_H1 <- list(response = data_all_H1$response, subj = data_all_H1$participant_ID, N = nrow(data_all_H1), J = length(unique(data_all_H1$participant_ID)))
# 3. where there is no ESP ability
stansimdata_all_H0 <- list(response = data_all_H0$response, subj = data_all_H0$participant_ID, N = nrow(data_all_H0), J = length(unique(data_all_H0$participant_ID)))


##########################################
#               Model fit                #
##########################################

# provide starting parameter values
initpar <- function() {
  list(beta = c(0,0), alpha = 0.97)
} 

# STAN model
model <-  "data {
            int<lower=1> N;                  //number of data points
            int<lower=0, upper=1> response[N];        // response
            int<lower=1> J;                  //number of subjects
            int<lower=1, upper=J> subj[N];   //subject id
            }
            
            parameters {
            real<lower=0, upper=1> alpha;  
            vector[2] beta;            //fixed means beta=[beta1,beta2]
            vector[J] beta_rand;               //subject intercepts
            }
            
            model {
            real  mu;
            
            // alpha ~ uniform(0.95,1);  //  Optional unform prior.	WARNING!!! This is a tight support! Stan might fail if the initial values provided to the stan() function are outside this interval
            alpha ~ beta(30,1.5);
            beta[1] ~ normal(0,0.01); // prior of beta1
            beta[2] ~ normal(0.05,0.01); // prior of beta2
            
            target += log_mix(alpha, normal_lpdf(beta_rand|beta[1], 0.01), normal_lpdf(beta_rand|beta[2], 0.01));    //subj random effects
            
            // likelihood
            for (i in 1:N){
            mu = exp(beta_rand[subj[i]])/(1+exp(beta_rand[subj[i]])) ;
            response[i] ~ bernoulli(mu);
            }
            }

"

# this writes the STAN model out into a file in the working directory
# this is recommended instead of using a string to specify model for optimal performance
write(model, file = "model.stan",
      ncolumns = if(is.character(model)) 1 else 5,
      append = FALSE, sep = " ")

# if you are using rstan locally on a multicore machine and have plenty of RAM to estimate your model in parallel
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())  # detect number of cores on the current processor


# Fit model with Monte Carlo sampling
# With 40000 data points (2000 simulated participants), 3 chains and 3000 iterations, this 
# ran for almost 4 hours and 40 minutes on Intel i7 6600U 2.6 GHz CPU.
Fit1 <- stan(file = "model.stan", data = stansimdata_part_H1, init=initpar, iter = 3000, chains = 3, control = list(adapt_delta = 0.999,max_treedepth=15))  # chains can be 1, 2, 3, 4 (or even more, but not recommended unless you have >4 cores at your disposal)
Fit2 <- stan(file = "model.stan", data = stansimdata_all_H1, init=initpar, iter = 3000, chains = 3, control = list(adapt_delta = 0.999,max_treedepth=15))  # chains can be 1, 2, 3, 4 (or even more, but not recommended unless you have >4 cores at your disposal)
Fit3 <- stan(file = "model.stan", data = stansimdata_all_H0, init=initpar, iter = 3000, chains = 3, control = list(adapt_delta = 0.999,max_treedepth=15))  # chains can be 1, 2, 3, 4 (or even more, but not recommended unless you have >4 cores at your disposal)

# Or the fit of a previous run can be loaded from a file:
# load(paste(file_location, "output_fit_mixture_alpha_has_betaprior", sep = "")) # loads the Stan output in the workspace


# plot traceplots, with and without warm-up:
traceplot(Fit1, pars = c("beta","alpha"),inc_warmup = TRUE)  # plot traces with burnin
traceplot(Fit1, pars = c("beta", "alpha"),inc_warmup = FALSE)  # plot traces without burnin
# examine quantiles of posterior distributions:
print(Fit1, pars = c("beta", "alpha"),probs = c(0.025, 0.5, 0.975))
#marginal posteriors and correlations
pairs(Fit1,pars = c("beta", "alpha"))


# plot traceplots, with and without warm-up:
traceplot(Fit2, pars = c("beta","alpha"),inc_warmup = TRUE)  # plot traces with burnin
traceplot(Fit2, pars = c("beta", "alpha"),inc_warmup = FALSE)  # plot traces without burnin
# examine quantiles of posterior distributions:
print(Fit2, pars = c("beta", "alpha"),probs = c(0.025, 0.5, 0.975))
#marginal posteriors and correlations
pairs(Fit2,pars = c("beta", "alpha"))


# plot traceplots, with and without warm-up:
traceplot(Fit3, pars = c("beta","alpha"),inc_warmup = TRUE)  # plot traces with burnin
traceplot(Fit3, pars = c("beta", "alpha"),inc_warmup = FALSE)  # plot traces without burnin
# examine quantiles of posterior distributions:
print(Fit3, pars = c("beta", "alpha"),probs = c(0.025, 0.5, 0.975))
#marginal posteriors and correlations
pairs(Fit3,pars = c("beta", "alpha"))
