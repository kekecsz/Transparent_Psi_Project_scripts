
#####################################################################
#              Conversion of effect size using formula              #
#####################################################################

# these data come from a recent meta-analysis: https://f1000research.com/articles/4-1188/v2
# Bem, D., Tressoldi, P., Rabeyron, T., & Duggan, M. (2015). Feeling the future: A meta-analysis of 90 experiments on the anomalous anticipation of random future events. F1000Research, 4.

g = 0.14 # hedges g effect size calculated using t-values from one-sample t-tests performed in the studies using this paradigm. Here a fixed effect model was used to estimate the effect size.
CI_95_lower_bound = 0.08 #lower bound of the 95% CI of the meta-nalaytic effect size estimate
n = 863 # the total pooled smaple size in the 14 studies using this paradigm
CI_95 = g-CI_95_lower_bound # the length of one arm of the CI
SE = CI_95/1.96 # standard error of the mean calculated from the CI
SD = SE*sqrt(n) # standard deviation calculated from the SE and the sample size (I am not sure we need this)
Var = SD^2 # variance derived from the SD

# conversion formula from cohen's d to log odds ratio can be found here: https://www.meta-analysis.com/downloads/Meta-analysis%20Converting%20among%20effect%20sizes.pdf
# same formula is used here: https://cran.r-project.org/web/packages/compute.es/compute.es.pdf
logoddsratio = g*(pi/sqrt(3))
Varlogoddsratio = Var*((pi^2)/3)

# BUT it returns too high p (proportion of succesful hits, this should be much lower, around 51-53% hit rate)
p = exp(logoddsratio)/(1+exp(logoddsratio))
p



# so I tried to figure out the p, and logOR associated with the hedges g effect sizes using simulation
############################################
#              Set parameters              #
############################################


# number of participants (pooled sample size in the meta-analysis)
n = 863
# number of trials performed per participant
trial_size_per_participant = 20
# number of trials to simulate in total
max_num_trials = n*trial_size_per_participant
# proportion of hits in the total sample if the null hypothesis is true
H0_prob = 0.5


############################################
#              Run simulation              #
############################################


es_data = data.frame(ptotry = seq(0.508, 0.524, 0.001), g = NA)
es_data = cbind(sapply(es_data[,"ptotry"], function(x) log((x/(1-x)))), es_data)
names(es_data)[1] = "logOR"

for(i in 1:nrow(es_data)){
  print(paste("trying", es_data[i,"ptotry"]))
  H1_prob = es_data[i,"ptotry"]
  g = NA  
  for(j in 1:1000){
    ##### no individual differences in performance
    # simulate data with the assumption that the null hypothesis is true (no individual differences are possible in this scenario, 
    # as each trial s just a random guess and everyone's guess is just as good as the others')
    data_all_H1_pre <- sample(c(rep(1, round(H1_prob*max_num_trials,0)), rep(0, round((1-H1_prob)*max_num_trials,0))))
    data_all_H1 = split(data_all_H1_pre, ceiling(seq_along(data_all_H1_pre)/trial_size_per_participant))
    
    mean(unlist(data_all_H1))
    
    success_proportions_all_H1 = sapply(data_all_H1, mean)
    t = t.test(success_proportions_all_H1, mu = H0_prob, alternative = "greater")$statistic
    
    g[j] = t/sqrt(n-1) # the 
  }
  es_data[i,"g"] = mean(g)
}

############################################
#                  Results                 #
############################################

#This shows that the logOR associated with g = 0.14 (95%CI 0.08-0.21) is logOR = 0.064 (95%CI 0.036-0.094)
es_data
