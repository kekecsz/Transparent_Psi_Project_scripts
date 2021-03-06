## The code below is used to visualize operational characteristics
## of the Transparent Psi Project stage 2 replication study

########################################################
#                     Load packages                    #
########################################################

library(HDInterval) # needed to calcluate HDI credible interval for visualization
library(reshape2)
library(tidyverse)

########################################################
#                    Custom Functions                  #
########################################################

# Function calculating the highest density interval using sampling
# We use hdi() from the library(HDInterval)
# needed for visualization of the power analysis results
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




# round to nearest x function
# this is to make sure that rounding does not result in a participant who has half of her trials with null effect and half with the virtuoso effect.
mround <- function(x,base){ 
  base*round(x/base) 
} 




###########################################################################
#   Visualize operational characteristics for all effect sizes (method 1) #
###########################################################################

# The file contains the results of 50000 simulations with 4 sequential
# The Bayes factor (01) threshold is 25 and 1/25 
# for stopping with a decision to support the M0 and M1 respectively. 
# The following true probabilities were simulated, each 50000 times: 
# 0.45, 0.48, 0.49, 0.495, 0.498, 0.499, 0.5, 0.501, 0.502, 0.503, 
# 0.504, 0.505, 0.506, 0.507, 0.508, 0.509, 0.510, 0.511, 0.512, 
# 0.515, 0.520, 0.530, 0.560.
output_frame <- read.csv(url("https://raw.githubusercontent.com/kekecsz/Transparent_Psi_Project_scripts/master/Operational_characteristics_all_ESs.csv"))
### Data for visualizing the inference decisions derived from the simulation
# data to visualize
wide <- output_frame[,c("True_prob", "Supported_H1_BF", "Supported_H0_BF", "Inconclusive_BF")]

# convert to long format
long <- reshape(wide, idvar = "True_prob", varying = c("Supported_H1_BF", "Supported_H0_BF", "Inconclusive_BF"),
                v.names = "conc", direction = "long")

# names of the inference decisions
long[long[,"time"] == 1,"time"] <- "evidence for M1"
long[long[,"time"] == 2,"time"] <- "evidence for M0"
long[long[,"time"] == 3,"time"] <- "inconclusive"
long[,"time"] <- factor(long[,"time"])

# ordered so "Inconclusive" will be displayed between support M1 and support M0 in the chart
long[,"time"] <- ordered( long[,"time"], levels = c("evidence for M0", "inconclusive", "evidence for M1"))

# variable names
names(long) <- c("Simulated_probability_of_success", "Inference", "Percentage")

# detect the plotted area (range of the simulated true probability values, the x axis)
plotted_area <- range(long$Simulated_probability_of_success) 

### For visualizing the prior distribution

# prior density, based on data gathered in Bem 2011 experiment 1
# this density plot will be overlayed on the inference decisions plot
scale <- seq(0,1, length = 100000)
alpha <- 829
beta <- 1560-828+1
interval <- c(0.5, 1)
normalizing_constant <- diff(pbeta(interval, alpha, beta))
range <- scale>interval[1] & scale<interval[2]
# this is the final density vector
prior_density <- range*dbeta(scale, alpha, beta)/normalizing_constant

# data for the density plot (scale and density vector) is put in a data frame and cut to fit the plotted area
plot_data <- as.data.frame(cbind(scale, (prior_density/(max(prior_density)*(3/2)))))
names(plot_data)[2] <- "density"
plot_data <- plot_data[plot_data$scale >= plotted_area[1] & plot_data$scale <= plotted_area[2],]

# 90% HDI is calculated for the prior distribution
mode_and_HDI <- mode_HDI(scale = scale, density = prior_density, crit_width = 0.90, n_samples = 1e5)
Bem_HDI_lb <- mode_and_HDI[2]
Bem_HDI_ub <- mode_and_HDI[3]


# visualization of how the probability
# of the different inference decisions change with the changes
# in the simulated true effect size (probability of success displayed
# on the x axis), and to put this in perspective we can also show
# the credible intervals of the previous study. 
# The lower bound of the 90% highest density credible interval is located at around p = 0.51

# Plots used as inspiration:
# Ly, A., Etz, A., & Wagenmakers, E. J. (2017). Replication Bayes Factors from Evidence Updating. Figure 1

# Creating a color palette for the plots
greypalette <- gray.colors(3, 0.5, 0.9, gamma = 2.2, alpha = NULL)

# Creating costum x axis labels and breaks
xlabs <- seq(head(long$Simulated_probability_of_success,1), tail(long$Simulated_probability_of_success,1), 0.025)

# Plot
figure_1 <- ggplot(long, aes(Simulated_probability_of_success, Percentage))+
  geom_area(aes(fill= Inference), position = 'stack')+
  geom_line(data = plot_data, aes(x=scale,y=density), linetype = 1, size = 1)+
  geom_vline(xintercept=interval[1], linetype = "longdash", size = 1) +
  scale_x_continuous(expand=c(0,0), breaks = xlabs, labels = xlabs) +
  scale_y_continuous(expand=c(0,0)) +
  geom_errorbarh(aes(xmax = as.numeric(Bem_HDI_ub), xmin = as.numeric(Bem_HDI_lb), y = 0.75, height = .05), size = 1.2) +
  annotate("text", x = 0.53, y = 0.79, label = "90% HDI: [0.51,0.551]", size = 4, fontface = 2) +
  geom_vline(xintercept=Bem_HDI_lb, linetype = "dashed", size = 1, alpha = 0.2)+
  geom_vline(xintercept=Bem_HDI_ub, linetype = "dashed", size = 1, alpha = 0.2) +
  scale_fill_manual(values = greypalette) +
  ylab("Proportion") +
  xlab("Simulated probability of success") +
  theme(axis.line = element_line(colour = "black", size = 1.2),
        axis.text.x = element_text(size = 12, face = "bold", margin = margin(4,0,0,0,"mm"), color = "black"),
        axis.text.y = element_text(size = 12, face = "bold", margin = margin(0,4,0,0,"mm"), color = "black"),
        axis.title = element_text(size = 16),
        axis.ticks = element_line(size = 1.2, linetype = 'solid', color = "black"),
        axis.ticks.length = unit(2, "mm"),
        legend.text = element_text(size = 10, face = "bold"),
        legend.title = element_blank(),
        legend.background = element_blank())

figure_1

# to save the plot as high res image
# tiff(file= "ES.tiff", width=10, height=6, units="in", pointsize="12", compression = "lzw", bg="white", res=600, antialias = "none")
# figure_1
# dev.off()
tiff(file= "C:\\Users\\kekec\\Dropbox\\_Active projects - Credibility Workgroup\\001 - TPP - Transparent Psi Project\\_Shared (Bem Replication)\\Stage 2 - Replication study\\R script\\Figures\\Supplement_Fig_1.tiff", width=10, height=6, units="in", pointsize="12", compression = "lzw", bg="white", res=600, antialias = "none")
figure_1
dev.off()

###############################################################################
#   Visualize operational characteristics for sampled effect sizes (method 2) #
###############################################################################

## data frames containing BFs and simulated true probability of success for each  of the simulations at stopping

# BFs and True_probs from simulations when True probability
# of success was sampled from the prior distribution for each study
BFs_sample <- read.csv(url("https://raw.githubusercontent.com/kekecsz/Transparent_Psi_Project_scripts/master/Simulation_Bayes_factors_sample.csv"))

# BFs and True_probs from simulations when M0 was simulated to be true
BFs_M0 <- read.csv(url("https://raw.githubusercontent.com/kekecsz/Transparent_Psi_Project_scripts/master/Simulation_Bayes_factors_050.csv"))

#### Plot for M1 (sampled effect sizes) with BF01

# Inspiration for the plot:
# Schonbrodt, F. D., & Wagenmakers, E. J. (2016). Bayes factor design analysis: Planning for compelling evidence. Psychonomic Bulletin & Review, 1-15. Figure 3a

# Creating thresholds for inference decisions
BF_threshold_high <- 25
BF_threshold_low <- 1/BF_threshold_high

# Creating a variable containing the inference of the given BF based on the thresholdss
BFs_sample$inference <- NA

BFs_sample[, "inference"] = apply(BFs_sample[,c("BF_replication", "BF_uniform", "BF_BUJ")], 1, function(x) if(min(x) >= BF_threshold_high){"evidence for M0"} else if(max(x) <= BF_threshold_low){"evidence for M1"} else {"Inconclusive"})

# Selecting the most conservative prior BF
BFs_sample$BF_most_conservative <- apply(BFs_sample[,c("BF_replication", "BF_uniform", "BF_BUJ")], 1, max)

# Creating axis labels, breaks and limits
BF_x_labels_M1 <- c("1/1000", paste("1/", BF_threshold_high, sep = ""), 1, BF_threshold_high, 1000)
BF_x_breaks_M1 <- log(c(1/1000, BF_threshold_low, 1, BF_threshold_high, 1000))

BF_y_labels_M1 <- c(0, 0.25, 0.5, 0.75, 0)

xlims_M1 <- c(-70,20)
ylims_M1 <- c(0, 0.01)

# Calculating density for the most conservative BF
dens_M1 <- density(log(BFs_sample[, "BF_most_conservative"])) # BFs are logged transformed for better visualization
plot_dataM1 <- data.frame(x = dens_M1$x, y = dens_M1$y)
plot_dataM1 <- plot_dataM1[plot_dataM1[, "x"] >= xlims_M1[1] & plot_dataM1[, "x"] <= xlims_M1[2], ]
thresholds <- c(log(BF_threshold_low), log(BF_threshold_high))
plot_dataM1$thresholds <- factor(findInterval(plot_dataM1$x, thresholds))
# Transforming thresholds to text for legend
plot_dataM1$legend <- ifelse(plot_dataM1$thresholds == 0, "true evidence for M1",
                             ifelse(plot_dataM1$thresholds == 1, "inconclusive evidence", "false evidence for M0"))

# Creating percentages of inference
false_M1 <- round(mean(BFs_sample[, "inference"] == "evidence for M0")*100, 2)
incon_M1 <- round(mean(BFs_sample[, "inference"] == "Inconclusive")*100, 2)
true_M1 <- round(mean(BFs_sample[, "inference"] == "evidence for M1")*100, 2)

false_M1+incon_M1+true_M1

# Plot
figure_2_b <- ggplot(plot_dataM1, aes(x,y)) +
  geom_line(size = 2) +
  geom_ribbon(aes(ymin = 0, ymax = y, fill = legend)) +
  geom_vline(xintercept = log(BF_threshold_low), linetype = "dashed", size = 1.2, color = "black") +
  geom_vline(xintercept = log(BF_threshold_high), linetype = "dashed", size = 1.2, color = "black") +
  annotate("text", x = (xlims_M1[1]-log(BF_threshold_low))/2, y = 0.005, label = paste("true evidence for M1", true_M1, "%", ""), angle = 90, fontface = 2) +
  annotate("text", x = 1, y = 0.0076, label = paste("inconclusive evidence ",incon_M1, "%", ""), angle = 90, fontface = 2) +
  annotate("text", x = (xlims_M1[2]-log(BF_threshold_high))/2, y = 0.005, label = paste("false evidence for M0", false_M1, "%", ""), angle = 90, fontface = 2) +
  scale_x_continuous(limits = xlims_M1, breaks = BF_x_breaks_M1, labels = BF_x_labels_M1, expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), limits = ylims_M1) +
  scale_fill_brewer(guide = "none") +
  xlab("Bayes Factor (BF01)") +
  ylab("Density") +
  ggtitle("Under M1") +
  labs(subtitle = "with Bayes Factor of M0 against M1") +
  scale_fill_manual(values = greypalette) +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.line = element_line(colour = "black", size = 1.2),
        axis.text.x = element_text(angle = 90, color = "black", face = "bold", size = 12, margin = margin(1,0,0,0,"mm"), vjust = 0.5),
        axis.text.y = element_text(face = "bold", color = "black", size = 12),
        axis.title = element_text(size = 16),
        axis.ticks = element_line(size = 1.2, linetype = 'solid', color = "black"),
        axis.ticks.length = unit(2, "mm"),
        plot.title = element_text(hjust = 0.5, size = 22),
        plot.subtitle = element_text(hjust = 0.5, size = 14))

figure_2_b

# to save the plot
# tiff(file= "M1.tiff", width=9, height=6.5, units="in", pointsize="12", compression = "lzw", bg="white", res=600, antialias = "none")
# jpeg(file= "M1.jpeg", width=9, height=6.5, units="in", pointsize="12", bg="white", res=600, antialias = "none")
# figure_3_b
# dev.off()

#### Plot for M0 (sampled True_prob always p = 0.5)

# Creating thresholds for inference decisions
BF_threshold_high = 25
BF_threshold_low = 1/BF_threshold_high

BFs_M0$inference <- NA
BFs_M0[, "inference"] = apply(BFs_M0[,c("BF_replication", "BF_uniform", "BF_BUJ")], 1, function(x) if(min(x) >= BF_threshold_high){"evidence for M0"} else if(max(x) <= BF_threshold_low){"evidence for M1"} else {"Inconclusive"})

# Selecting the most conservative prior BF
BFs_M0$BF_most_conservative <- apply(BFs_M0[,c("BF_replication", "BF_uniform", "BF_BUJ")], 1, min)

# Creating axis labels, breaks and limits
xlims_M0 <- c(-6, 6)

BF_x_labels_M0 <- c("1/1000", paste("1/", BF_threshold_high, sep = ""), 1, BF_threshold_high, 1000)
BF_x_breaks_M0 <- log(c(1/1000, BF_threshold_low, 1, BF_threshold_high, 1000))

# Calculating density for the least conservative BF
dens_M0 <- density(log(BFs_M0[, "BF_most_conservative"])) # BFs are logged transformed for better visualization
plot_dataM0 <- data.frame(x = dens_M0$x, y = dens_M0$y)
plot_dataM0 <- plot_dataM0[plot_dataM0[, "x"] >= xlims_M0[1] & plot_dataM0[, "x"] <= xlims_M0[2], ]
thresholds <- c(log(BF_threshold_low), log(BF_threshold_high))
plot_dataM0$thresholds <- factor(findInterval(plot_dataM0$x, thresholds))
plot_dataM0$legend <- ifelse(plot_dataM0$thresholds == 0, "false evidence for M1",
                             ifelse(plot_dataM0$thresholds == 1, "inconclusive evidence", "true evidence for M0"))

false_M0 <- round(mean(BFs_M0[, "inference"] == "evidence for M1")*100, 2)
incon_M0 <- round(mean(BFs_M0[, "inference"] == "Inconclusive")*100, 2)
true_M0 <- round(mean(BFs_M0[, "inference"] == "evidence for M0")*100, 2)

false_M0+incon_M0+true_M0

# Plot
figure_2_a <- ggplot(plot_dataM0, aes(x,y)) +
  geom_line(size = 2) +
  geom_ribbon(aes(ymin = 0, ymax = y, fill = legend)) +
  geom_vline(xintercept = log(BF_threshold_low), linetype = "dashed", size = 1.2, color = "black") +
  geom_vline(xintercept = log(BF_threshold_high), linetype = "dashed", size = 1.2, color = "black") +
  annotate("text", x = (xlims_M0[1]+log(BF_threshold_low))/2, y = 0.4, label = paste("false evidence for M1", false_M0, "%", ""), angle = 90, fontface = 2) +
  annotate("text", x = 0, y = 0.4, label = paste("inconclusive evidence", incon_M0, "%", ""), angle = 90, fontface = 2) +
  annotate("text", x = (xlims_M0[2]+log(BF_threshold_high)-1)/2, y = 0.4, label = paste("true evidence for M0", true_M0, "%", ""), angle = 90, fontface = 2) +
  scale_x_continuous(breaks = BF_x_breaks_M0, labels = BF_x_labels_M0,expand = c(0,0), limits = xlims_M0) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 0.75)) +
  scale_fill_brewer(guide = "none") +
  xlab("Bayes Factor (BF01)") +
  ylab("Density") +
  ggtitle("Under M0") +
  labs(subtitle = "with Bayes Factor of M0 against M1") +
  scale_fill_manual(values = greypalette) +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.line = element_line(colour = "black", size = 1.2),
        axis.text.x = element_text(angle = 90, color = "black", face = "bold", size = 12, margin = margin(1,0,0,0,"mm"), vjust = 0.5),
        axis.text.y = element_text(face = "bold", color = "black", size = 12),
        axis.title = element_text(size = 16),
        axis.ticks = element_line(size = 1.2, linetype = 'solid', color = "black"),
        axis.ticks.length = unit(2, "mm"),
        plot.title = element_text(hjust = 0.5, size = 22),
        plot.subtitle = element_text(hjust = 0.5, size = 14))

figure_2_a

# to save the plot
# tiff(file= "M0.tiff", width=9, height=6.5, units="in", pointsize="12", compression = "lzw", bg="white", res=600, antialias = "none")
# jpeg(file= "M0.jpeg", width=9, height=6.5, units="in", pointsize="12", bg="white", res=600, antialias = "none")
# figure_3_a
# dev.off()


###########################################################################
#        Visualize Comparison of expected and observed distributions      #
###########################################################################

#               Set parameters for example data                      


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
ESP_user_percentage = 0.1 

# probability of successful guess among ESP-users 
M1_prob = 0.6

# chance of stopping prematurely in each trial, to simulate some missing data due to unexpected events
# setting this to 0.002 would generate roughly 2% missing data
# note that premature stopping cannot introduce bias in our confirmatory analyses
chance_for_stopping_session = 0



#                       Generate example data                        


# simulate the performance of 3,000 potential participants
# In this code we only simulate erotic trials
# Non-erotic trials will be ommitted from the hypothesis testing analysis 
# during data management.

list = list(NA)

for(i in 1:3000){
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

data_BF = do.call("rbind", list)
names(data_BF) = c("participant_ID", "success")



#                       Data management for plot                     


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



#                        Plot 

histogram_plot_data = as.data.frame(c(success_rates_theoretical_prop, success_rates_empirical_prop))
histogram_plot_data = cbind(histogram_plot_data, factor(c(rep("Expected if M0 is true", length(success_rates_theoretical_prop)), rep("Observed", length(success_rates_empirical_prop)))))
histogram_plot_data = cbind(histogram_plot_data, factor(rep(possible_success_rates_char, 2)))
names(histogram_plot_data) = c("proportion", "group", "success")

figure_3 =  ggplot(histogram_plot_data, aes(y = proportion, x = success, group = group))+
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
figure_3




# to save the plot as high res image
# jpeg(file= "figure_3.jpeg", width=9, height=6.5, units="in", pointsize="12", bg="white", res=600, antialias = "none")
# figure_3
# dev.off()

