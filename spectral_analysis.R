# SEMELI FRANGOPOULOU - SNR 2020452

# PART 1: SPECTRAL ANALYSIS

#installing packages
library(dplyr)
library(readr)
install.packages("tidyverse")
library(tidyverse)

# Load AD and HC datasets with mean power per band
ad_mean <- read_csv("ad_mean.csv")
hc_mean <- read_csv("hc_mean.csv")

# AD band division
ad_delta <- select(ad_mean, Delta)
ad_theta <- select(ad_mean, Theta)
ad_alpha <- select(ad_mean, Alpha)
ad_beta <- select(ad_mean, Beta)
ad_gamma <- select(ad_mean, Gamma)

# HC band division
hc_delta <- select(hc_mean, Delta)
hc_theta <- select(hc_mean, Theta)
hc_alpha <- select(hc_mean, Alpha)
hc_beta <- select(hc_mean, Beta)
hc_gamma <- select(hc_mean, Gamma)

# Checking for normality in AD bands -> not normal
shapiro.test(ad_mean$Delta)
shapiro.test(ad_mean$Theta)
shapiro.test(ad_mean$Alpha)
shapiro.test(ad_mean$Beta)
shapiro.test(ad_mean$Gamma)

# Checking for normality in HC bands -> not normal
shapiro.test(hc_mean$Delta)
shapiro.test(hc_mean$Theta)
shapiro.test(hc_mean$Alpha)
shapiro.test(hc_mean$Beta)
shapiro.test(hc_mean$Gamma)

# Running a non-parametric Mann Whitney-U test (Wilcoxon rank-sum test)
# Independent comparison for every frequency band between AD and HC
wilcox.test(ad_mean$Delta, hc_mean$Delta)
wilcox.test(ad_mean$Theta, hc_mean$Theta)
wilcox.test(ad_mean$Alpha, hc_mean$Alpha)
wilcox.test(ad_mean$Beta, hc_mean$Beta)
wilcox.test(ad_mean$Gamma, hc_mean$Gamma)