# https://cran.r-project.org/web/packages/brms/vignettes/brms_overview.pdf

library(brms)



### Function to round to nearest x

mround <- function(x,base){ 
  base*round(x/base) 
} 

### randomize psi talent for each person if personal differences are simulated
# True_prob is the population mean effect in probability of success
# SD_personal_differences is the standard deviation of the effect across individuals

talent_randomizer <- function(True_prob, SD_personal_differences){
  talent = True_prob + rnorm(mean = 0, sd = SD_personal_differences, n = 1)
  if(talent<0.001){talent = 0.001}else if(talent>0.999){talent = 0.999}
  return(talent)
}


### characteristics of the study
max_num_trials = 100080 # maximum number of trials that can be run in this study

True_prob = 0.5 # true mean probability of success in the population

SD_personal_differences = 0.15 # the standard deviation of the probability of success across individuals, set this to 0 to simulate no personal differences, set this to e.g. 0.15 to generate some personal differences in ability/talent

# randomly determin the precognition ability of each person, and simulate successes according to their respective ability
data_all_M1 <- as.vector(replicate(mround(max_num_trials, base = 18)/18, rbinom(18, size = 1, prob=talent_randomizer(True_prob = True_prob, SD_personal_differences = SD_personal_differences))))
# add ID for each participant
data_all_M1_dataframe = as.data.frame(cbind(data_all_M1, rep(1:(mround(max_num_trials, base = 18)/18), each = 18)))
data_all_M1_dataframe[,2] = as.factor(data_all_M1_dataframe[,2])
# name columns
names(data_all_M1_dataframe) = c("outcome","ID")

# Note that this code does not simulate the sequential analysis plan (here we only analyze once at the end of collecting the maximum amount of data)
# It also does not inclide the different priors or additional robustnes analyses planned.

# Fit brms model
fit2 <- brm(outcome ~ 0 + intercept + (0 + intercept|ID), 
            data = data_all_M1_dataframe, 
            family = bernoulli, 
            chains = 4, iter = 2000,
            control = list(adapt_delta = 0.90),
            prior = set_prior("normal(0,0.6)", class = "b"),
            sample_prior = TRUE, 
            save_all_pars = TRUE, 
            cores = 3)
pairs(fit2)
summary(fit2)
plot(fit2)
hypothesis(fit2, hypothesis = 'intercept = 0')


# or use this prior for a more informed prior about realistic personal differences 
# prior = c(set_prior("normal(0,0.6)", class = "b"), set_prior("normal(0,0.5)", class = "sd")),
# the default prior for sd in this model is set_prior("student_t(3, 0, 10)", class = "sd"))




