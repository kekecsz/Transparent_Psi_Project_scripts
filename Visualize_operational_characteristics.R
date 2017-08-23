
## The code below is used to visualize operational characteristics
## of the Transparent Psi Project stage 2 replication study

########################################################
#   Visualize inference decisions for all effect sizes #
########################################################


# The file contains the results of 50000 simulations with 4 sequential
# stopping points (at reaching 26000, 44000, 62000, and 80000 
# completed trials). The Bayes factor (01) threshold is 45 and 1/45 
# for stopping with a decision to support the H0 and H1 respectively. 
# The following true probabilities were simulated, each 50000 times: 
# 0.45, 0.48, 0.49, 0.495, 0.498, 0.499, 0.5, 0.501, 0.502, 0.503, 
# 0.504, 0.505, 0.506, 0.507, 0.508, 0.509, 0.510, 0.511, 0.512, 
# 0.515, 0.520, 0.530, 0.560.
output_frame <- read.csv(url("https://github.com/kekecsz/Transparent_Psi_Project_scripts/blob/master/operational_characteristics_50000sims_beta_BF_45_prior_BEM_allESs.csv"))



### Data for visualizing the inference decisions derived from the simulation

# data to visualize
wide <- output_frame[,c("True_prob", "Supported_H1_BF", "Supported_H0_BF", "Inconclusive_BF")]

# convert to long format
long = reshape(wide, idvar = "True_prob", varying = c("Supported_H1_BF", "Supported_H0_BF", "Inconclusive_BF"),
               v.names = "conc", direction = "long")

# names of the inference decisions
long[long[,"time"] == 1,"time"] = "Support H1"
long[long[,"time"] == 2,"time"] = "Support H0"
long[long[,"time"] == 3,"time"] = "Inconclusive"
long[,"time"] = factor(long[,"time"])
# ordered so "Inconclusive" will be displayed between support H1 and support H0 in the chart
long[,"time"] = ordered( long[,"time"], levels = c("Support H0", "Inconclusive", "Support H1"))
# variable names
names(long) = c("Simulated_probability_of_success", "Inference", "Percentage")
# detect the plotted area (range of the simulated true probability values, the x axis)
plotted_area = range(long$Simulated_probability_of_success) 


### For visualizing the prior distribution

# prior density, based on data gathered in Bem 2011 experiment 1
# this density plot will be overlayed on the inference decisions plot
scale = seq(0,1, length = 100000)
alpha = 829
beta = 1560-828+1
interval = c(0.5, 1)
normalizing_constant = diff(pbeta(interval, alpha, beta))
range = scale>interval[1] & scale<interval[2]
# this is the final density vector
prior_density = range*dbeta(scale, alpha, beta)/normalizing_constant

# data for the density plot (scale and density vector) is put in a data frame and cut to fit the plotted area
plot_data = as.data.frame(cbind(scale, (prior_density/(max(prior_density)*(3/2)))))
names(plot_data)[2] = "density"
plot_data = plot_data[plot_data$scale >= plotted_area[1] & plot_data$scale <= plotted_area[2],]

# 90% HDI is calculated for the prior distribution
mode_and_HDI = mode_HDI(scale = scale, density = prior_density, crit_width = 0.90, n_samples = 1e5)
Bem_HDI_lb = mode_and_HDI[2]
Bem_HDI_ub = mode_and_HDI[3]

### Plotting
# this is just a first sketch of the figure I imagined
# the important thing is to visualize the how the probability
# of te different inference decisions change with the changes
# in the simulated true effect size (probability of success displayed
# on the x axis), and to put this in perspective we can also show
# the credible intervals of the previous study. 
# The lower bound of the 90% highest density credible interval is
# conveniently located at around p = 0.51, where we still have good
# operational characteristics (roughly 0.95 power and 0.005 alpha error for both H0 and H1)
figure <- ggplot(long, aes( Simulated_probability_of_success, Percentage))+
                 geom_area(aes(colour = Inference, fill= Inference), position = 'stack')+
                 geom_line(data = plot_data, aes(x=scale,y=density), linetype = 1)+
                 geom_vline(xintercept=interval[1], linetype = "solid")+
                 geom_vline(xintercept=Bem_HDI_lb, linetype = "dashed")+
                 geom_vline(xintercept=Bem_HDI_ub, linetype = "dashed")
figure

#############################################################
#   Visualize inference decisions for sampled effect sizes  #
#############################################################

## data frames containing BF and simulated true probability of success for each  of the simulations at stopping

# BFs and True_probs from simulations when True probability
# of success was sampled from the prior distribution for each study
BFs_sample <- read.csv(url("https://github.com/kekecsz/Transparent_Psi_Project_scripts/blob/master/..."))

# BFs and True_probs from simulations when H0 was simulated to be true
BFs_H0 <- read.csv(url("https://github.com/kekecsz/Transparent_Psi_Project_scripts/blob/master/..."))

# Data frame containing the operational characteristics, and the parameters of the two above mentioned simulations
# this is not necesseraly needed for the visualization
output_frame_sample_and_H0 <- read.csv(url("https://github.com/kekecsz/Transparent_Psi_Project_scripts/blob/master/..."))


#### possible plots for H1 (sampled effect sizes)

# checking that p-s were properly sampled from the prior beta distribution
# the histogram should be roughly following the curve
hist(BFs_sample[,"True_prob"], freq = F, breaks = 100)
curve(prior_density, add = T)

# very primitive visualization of the BF
plot(density(log(BFs_sample[,"BF"])))

#### possible plots for H0 (sampled True_prob always p = 0.5)

# very primitive visualization of the BF
plot(density(log(BFs_H0[,"BF"])))
