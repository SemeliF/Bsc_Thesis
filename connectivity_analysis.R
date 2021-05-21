# SEMELI FRANGOPOULOU - SNR 2020452


# PART 2: CONNECTIVITY ANALYSIS (running statistical tests)
#installing packages
#installing packages

install.packages("tidyverse")
install.packages("devtools")
devtools::install_github("hrbrmstr/hrbrthemes")
install.packages("lsr")
library(tidyverse)
library(dplyr)
library(readr)
library(devtools)
library(hrbrthemes)
library(viridisLite)
library(lsr)
library(ggpubr)

## STEP 1: Load and read files
# importing all AD files and storing them into separate data frames
my_ad_files <- list.files("C:/Users/Semeli/Documents/R/AD", pattern = ".csv", full.names= TRUE)
my_ad_dfs <- lapply(my_ad_files, read.csv)
names(my_ad_dfs) <- my_ad_files

# importing all HC files and storing them into separate data frames
my_hc_files <- list.files("C:/Users/Semeli/Documents/R/HC", pattern = ".csv", full.names= TRUE)
my_hc_dfs <- lapply(my_hc_files, read.csv)
names(my_hc_dfs) <- my_hc_files

################################################################################

## STEP 2: Band mean calculation for every subject
# calculating band means for every AD subject and appending them into a list
ad_mean_per_band <- list()
for (i in my_ad_dfs){
  tmp <- rowMeans(i[,-1])
  ad_mean_per_band <- append(ad_mean_per_band, list(tmp))
}

# calculating band means for every HC subject and appending them into a list
hc_mean_per_band <- list()
for (i in my_hc_dfs){
  tmp <- rowMeans(i[,-1])
  hc_mean_per_band <- append(hc_mean_per_band, list(tmp))
}

## STEP 3: Average band mean calculation for AD vs HC subjects
#order: delta, theta, alpha, beta, gamma
# mean AD delta band
ad_delta_vals <- list()
for (subj in ad_mean_per_band) {
  tmp_ad_delta <- subj[1]
  ad_delta_vals <- append(ad_delta_vals, list(tmp_ad_delta))
}
ad_mean_delta <- mean(unlist(ad_delta_vals))

#mean AD theta band
ad_theta_vals <- list()
for (subj in ad_mean_per_band) {
  tmp_ad_theta <- subj[2]
  ad_theta_vals <- append(ad_theta_vals, list(tmp_ad_theta))
}
ad_mean_theta <- mean(unlist(ad_theta_vals))

#mean AD alpha band
ad_alpha_vals <- list()
for (subj in ad_mean_per_band) {
  tmp_ad_alpha <- subj[3]
  ad_alpha_vals <- append(ad_alpha_vals, list(tmp_ad_alpha))
}
ad_mean_alpha <- mean(unlist(ad_alpha_vals))

#mean AD beta band
ad_beta_vals <- list()
for (subj in ad_mean_per_band) {
  tmp_ad_beta <- subj[4]
  ad_beta_vals <- append(ad_beta_vals, list(tmp_ad_beta))
}
ad_mean_beta <- mean(unlist(ad_beta_vals))

#mean AD gamma band
ad_gamma_vals <- list()
for (subj in ad_mean_per_band) {
  tmp_ad_gamma <- subj[5]
  ad_gamma_vals <- append(ad_gamma_vals, list(tmp_ad_gamma))
}
ad_mean_gamma <- mean(unlist(ad_gamma_vals))


# mean HC delta band
hc_delta_vals <- list()
for (subj in hc_mean_per_band) {
  tmp_hc_delta <- subj[1]
  hc_delta_vals <- append(hc_delta_vals, list(tmp_hc_delta))
}
hc_mean_delta <- mean(unlist(hc_delta_vals))

#mean HC theta band
hc_theta_vals <- list()
for (subj in hc_mean_per_band) {
  tmp_hc_theta <- subj[2]
  hc_theta_vals <- append(hc_theta_vals, list(tmp_hc_theta))
}
hc_mean_theta <- mean(unlist(hc_theta_vals))

#mean HC alpha band
hc_alpha_vals <- list()
for (subj in hc_mean_per_band) {
  tmp_hc_alpha <- subj[3]
  hc_alpha_vals <- append(hc_alpha_vals, list(tmp_hc_alpha))
}
hc_mean_alpha <- mean(unlist(hc_alpha_vals))

#mean HC beta band
hc_beta_vals <- list()
for (subj in hc_mean_per_band) {
  tmp_hc_beta <- subj[4]
  hc_beta_vals <- append(hc_beta_vals, list(tmp_hc_beta))
}
hc_mean_beta <- mean(unlist(hc_beta_vals))

#mean HC gamma band
hc_gamma_vals <- list()
for (subj in hc_mean_per_band) {
  tmp_hc_gamma <- subj[5]
  hc_gamma_vals <- append(hc_gamma_vals, list(tmp_hc_gamma))
}
hc_mean_gamma <- mean(unlist(hc_gamma_vals))

###############################################################################

## STEP 4 : Box plots per band for comparison + Mann-Whitney U-test
# Delta Connectivity

delta_con <- data.frame(AD = unlist(ad_delta_vals),
                        HC = unlist(hc_delta_vals))
boxplot(delta_con,
        ylab = "Phase-Locking Value (PLV)",
        main = "Delta Connectivity",
        border = "blue")


# Theta Connectivity
theta_con <- data.frame(AD = unlist(ad_theta_vals),
                        HC = unlist(hc_theta_vals))
boxplot(theta_con,
        ylab = "Phase-Locking Value (PLV)",
        main = "Theta Connectivity",
        border = "red")


# Alpha Connectivity
alpha_con <- data.frame(AD = unlist(ad_alpha_vals),
                        HC = unlist(hc_alpha_vals))
boxplot(alpha_con,
        ylab = "Phase-Locking Value (PLV)",
        main = "Alpha Connectivity",
        border = "dark green")


# Beta Connectivity
beta_con <- data.frame(AD = unlist(ad_beta_vals),
                       HC = unlist(hc_beta_vals))
boxplot(beta_con,
        ylab = "Phase-Locking Value (PLV)",
        main = "Beta Connectivity",
        border = "magenta")


# Gamma Connectivity
gamma_con <- data.frame(AD = unlist(ad_gamma_vals),
                        HC = unlist(hc_gamma_vals))
boxplot(gamma_con,
        ylab = "Phase-Locking Value (PLV)",
        main = "Gamma Connectivity",
        border = "orange")


## Mann-Whitney U-test between the groups 5 times, each time on the AVG PLV in one freq. band
wilcox.test(unlist(ad_delta_vals), unlist(hc_delta_vals))
wilcox.test(unlist(ad_theta_vals), unlist(hc_theta_vals)) # -> SIGNIFICANCE
wilcox.test(unlist(ad_alpha_vals), unlist(hc_alpha_vals))
wilcox.test(unlist(ad_beta_vals), unlist(hc_beta_vals))
wilcox.test(unlist(ad_gamma_vals), unlist(hc_gamma_vals))

###############################################################################
## STEP 6: evaluating differences for c3-p3 & c4-p4 -> PAIR A
# loading data
ad_c3p3_c4p4 <- read_csv("ad_c3p3-c4p4_2.csv")
hc_c3p3_c4p4 <- read_csv("hc_c3p3-c4p4.csv")

# Wilcoxon rank-sum test per band between AD and HC
wilcox.test(ad_c3p3_c4p4$delta, hc_c3p3_c4p4$delta)
wilcox.test(ad_c3p3_c4p4$theta, hc_c3p3_c4p4$theta) #significance but AD > HC
wilcox.test(ad_c3p3_c4p4$beta, hc_c3p3_c4p4$beta)
wilcox.test(ad_c3p3_c4p4$alpha, hc_c3p3_c4p4$alpha)
wilcox.test(ad_c3p3_c4p4$gamma, hc_c3p3_c4p4$gamma)

# plotting the wilcoxon tests for every band
# delta
pairA_delta <- data.frame(AD = ad_c3p3_c4p4$delta,
                        HC = hc_c3p3_c4p4$delta)
boxplot(pairA_delta,
        ylab = "Phase-Locking Value (PLV)",
        main = "Delta",
        border = "blue")
# theta
pairA_theta <- data.frame(AD = ad_c3p3_c4p4$theta,
                          HC = hc_c3p3_c4p4$theta)
boxplot(pairA_theta,
        ylab = "Phase-Locking Value (PLV)",
        main = "Theta *",
        border = "red")
# alpha
pairA_alpha <- data.frame(AD = ad_c3p3_c4p4$alpha,
                          HC = hc_c3p3_c4p4$alpha)
boxplot(pairA_alpha,
        ylab = "Phase-Locking Value (PLV)",
        main = "Alpha",
        border = "dark green")
# beta
pairA_beta <- data.frame(AD = ad_c3p3_c4p4$beta,
                          HC = hc_c3p3_c4p4$beta)
boxplot(pairA_beta,
        ylab = "Phase-Locking Value (PLV)",
        main = "Beta",
        border = "magenta")
# gamma
pairA_gamma <- data.frame(AD = ad_c3p3_c4p4$gamma,
                          HC = hc_c3p3_c4p4$gamma)
boxplot(pairA_gamma,
        ylab = "Phase-Locking Value (PLV)",
        main = "Gamma",
        border = "orange")

###############################################################################
## STEP 7: evaluating differences for c3-f3 & c4-f4 -> PAIR B
# loading data
ad_f3c3_f4c4 <- read_csv("ad_f3c3_f4c4.csv")
hc_f3c3_f4c4 <- read_csv("hc_f3c3_f4c4.csv")

# Wilcoxon rank-sum test per band between AD and HC
wilcox.test(ad_f3c3_f4c4$delta, hc_f3c3_f4c4$delta) #significance but AD > HC
wilcox.test(ad_f3c3_f4c4$theta, hc_f3c3_f4c4$theta) #significance but AD > HC
wilcox.test(ad_f3c3_f4c4$beta, hc_f3c3_f4c4$beta)
wilcox.test(ad_f3c3_f4c4$alpha, hc_c3p3_c4p4$alpha)
wilcox.test(ad_c3p3_c4p4$gamma, hc_f3c3_f4c4$gamma)

# plotting the wilcoxon tests for every band
# delta
pairB_delta <- data.frame(AD = ad_f3c3_f4c4$delta,
                          HC = hc_f3c3_f4c4$delta)
boxplot(pairB_delta,
        ylab = "Phase-Locking Value (PLV)",
        main = "Delta *",
        border = "blue")
# theta
pairB_theta <- data.frame(AD = ad_f3c3_f4c4$theta,
                          HC = hc_f3c3_f4c4$theta)
boxplot(pairB_theta,
        ylab = "Phase-Locking Value (PLV)",
        main = "Theta *",
        border = "red")
# alpha
pairB_alpha <- data.frame(AD = ad_f3c3_f4c4$alpha,
                          HC = hc_f3c3_f4c4$alpha)
boxplot(pairB_alpha,
        ylab = "Phase-Locking Value (PLV)",
        main = "Alpha",
        border = "dark green")
# beta
pairB_beta <- data.frame(AD = ad_f3c3_f4c4$beta,
                          HC = hc_f3c3_f4c4$beta)
boxplot(pairB_beta,
        ylab = "Phase-Locking Value (PLV)",
        main = "Beta",
        border = "magenta")
# gamma
pairB_gamma <- data.frame(AD = ad_f3c3_f4c4$gamma,
                          HC = hc_f3c3_f4c4$gamma)
boxplot(pairB_gamma,
        ylab = "Phase-Locking Value (PLV)",
        main = "Gamma",
        border = "orange")

###############################################################################
## STEP 8: evaluating differences for p4-o2 & p3-o1 -> PAIR C
# loading data
ad_p4o2_p3o1 <- read_csv("ad_p4o2_p3o1.csv")
hc_p4o2_p3o1 <- read_csv("hc_p4o2_p3o1.csv")

wilcox.test(ad_p4o2_p3o1$delta, hc_p4o2_p3o1$delta) # -> SIGNIFICANCE
wilcox.test(ad_p4o2_p3o1$theta, hc_p4o2_p3o1$theta)
wilcox.test(ad_p4o2_p3o1$beta, hc_p4o2_p3o1$beta)
wilcox.test(ad_p4o2_p3o1$alpha, hc_p4o2_p3o1$alpha)
wilcox.test(ad_p4o2_p3o1$gamma, hc_p4o2_p3o1$gamma)

# plotting the wilcoxon tests for every band
# delta
pairC_delta <- data.frame(AD = ad_p4o2_p3o1$delta,
                          HC = hc_p4o2_p3o1$delta)
boxplot(pairC_delta,
        ylab = "Phase-Locking Value (PLV)",
        main = "Delta *",
        border = "blue")
# theta
pairC_theta <- data.frame(AD = ad_p4o2_p3o1$theta,
                          HC = hc_p4o2_p3o1$theta)
boxplot(pairC_theta,
        ylab = "Phase-Locking Value (PLV)",
        main = "Theta",
        border = "red")
# alpha
pairC_alpha <- data.frame(AD = ad_p4o2_p3o1$alpha,
                          HC = hc_p4o2_p3o1$alpha)
boxplot(pairC_alpha,
        ylab = "Phase-Locking Value (PLV)",
        main = "Alpha",
        border = "dark green")
# beta
pairC_beta <- data.frame(AD = ad_p4o2_p3o1$beta,
                          HC = hc_p4o2_p3o1$beta)
boxplot(pairC_beta,
        ylab = "Phase-Locking Value (PLV)",
        main = "Beta",
        border = "magenta")
# gamma
pairC_gamma <- data.frame(AD = ad_p4o2_p3o1$gamma,
                          HC = hc_p4o2_p3o1$gamma)
boxplot(pairC_gamma,
        ylab = "Phase-Locking Value (PLV)",
        main = "Gamma",
        border = "orange")

###############################################################################
## STEP 9: evaluating differences for T4-C4 & T3-C3 -> PAIR D
# loading data
hc_t4c4_t3c3 <- read_csv("hc_t4c4_t3c3.csv")
ad_t4c4_t3c3 <- read_csv("ad_t4c4_t3c3.csv")

wilcox.test(ad_t4c4_t3c3$delta, hc_t4c4_t3c3$delta)
wilcox.test(ad_t4c4_t3c3$theta, hc_t4c4_t3c3$theta) # -> SIGNIFICANCE
wilcox.test(ad_t4c4_t3c3$beta, hc_t4c4_t3c3$beta)
wilcox.test(ad_t4c4_t3c3$alpha, hc_t4c4_t3c3$alpha)
wilcox.test(ad_t4c4_t3c3$gamma, hc_t4c4_t3c3$gamma)

# plotting the wilcoxon tests for every band
# delta
pairD_delta <- data.frame(AD = ad_t4c4_t3c3$delta,
                          HC = hc_t4c4_t3c3$delta)
boxplot(pairD_delta,
        ylab = "Phase-Locking Value (PLV)",
        main = "Delta",
        border = "blue")
# theta
pairD_theta <- data.frame(AD = ad_t4c4_t3c3$theta,
                          HC = hc_t4c4_t3c3$theta)
boxplot(pairD_theta,
        ylab = "Phase-Locking Value (PLV)",
        main = "Theta *",
        border = "red")
# alpha
pairD_alpha <- data.frame(AD = ad_t4c4_t3c3$alpha,
                          HC = hc_t4c4_t3c3$alpha)
boxplot(pairD_alpha,
        ylab = "Phase-Locking Value (PLV)",
        main = "Alpha",
        border = "dark green")
# beta
pairD_beta <- data.frame(AD = ad_t4c4_t3c3$beta,
                          HC = hc_t4c4_t3c3$beta)
boxplot(pairD_beta,
        ylab = "Phase-Locking Value (PLV)",
        main = "Beta",
        border = "magenta")
# gamma
pairD_gamma <- data.frame(AD = ad_t4c4_t3c3$gamma,
                          HC = hc_t4c4_t3c3$gamma)
boxplot(pairD_gamma,
        ylab = "Phase-Locking Value (PLV)",
        main = "Gamma",
        border = "orange")

###############################################################################
## STEP 10: comparing pairs between them
## AD
# delta band
wilcox.test(ad_c3p3_c4p4$delta, ad_f3c3_f4c4$delta,paired = TRUE)
wilcox.test(ad_c3p3_c4p4$delta, ad_p4o2_p3o1$delta, paired = TRUE) # A < C
wilcox.test(ad_c3p3_c4p4$delta, ad_t4c4_t3c3$delta, paired = TRUE)
wilcox.test(ad_f3c3_f4c4$delta, ad_p4o2_p3o1$delta, paired = TRUE)
wilcox.test(ad_f3c3_f4c4$delta, ad_t4c4_t3c3$delta, paired = TRUE) # B > D
wilcox.test(ad_p4o2_p3o1$delta, ad_t4c4_t3c3$delta, paired = TRUE) # C > D

# theta band
wilcox.test(ad_c3p3_c4p4$theta, ad_f3c3_f4c4$theta, paired = TRUE)
wilcox.test(ad_c3p3_c4p4$theta, ad_p4o2_p3o1$theta, paired = TRUE)
wilcox.test(ad_c3p3_c4p4$theta, ad_t4c4_t3c3$theta, paired = TRUE) # A > D
wilcox.test(ad_f3c3_f4c4$theta, ad_p4o2_p3o1$theta, paired = TRUE)
wilcox.test(ad_f3c3_f4c4$theta, ad_t4c4_t3c3$theta, paired = TRUE) # B > D
wilcox.test(ad_p4o2_p3o1$theta, ad_t4c4_t3c3$theta, paired = TRUE)

## HC

# delta band
wilcox.test(hc_c3p3_c4p4$delta, hc_f3c3_f4c4$delta, paired = TRUE)
wilcox.test(hc_c3p3_c4p4$delta, hc_p4o2_p3o1$delta, paired = TRUE)
wilcox.test(hc_c3p3_c4p4$delta, hc_t4c4_t3c3$delta, paired = TRUE)
wilcox.test(hc_f3c3_f4c4$delta, hc_p4o2_p3o1$delta, paired = TRUE)
wilcox.test(hc_f3c3_f4c4$delta, hc_t4c4_t3c3$delta, paired = TRUE)
wilcox.test(hc_p4o2_p3o1$delta, hc_t4c4_t3c3$delta, paired = TRUE) # C > D

# theta band
wilcox.test(hc_c3p3_c4p4$theta, hc_f3c3_f4c4$theta, paired = TRUE)
wilcox.test(hc_c3p3_c4p4$theta, hc_p4o2_p3o1$theta, paired = TRUE)
wilcox.test(hc_c3p3_c4p4$theta, hc_t4c4_t3c3$theta, paired = TRUE) # A > D
wilcox.test(hc_f3c3_f4c4$theta, hc_p4o2_p3o1$theta, paired = TRUE)
wilcox.test(hc_f3c3_f4c4$theta, hc_t4c4_t3c3$theta, paired = TRUE) # B > D
wilcox.test(hc_p4o2_p3o1$theta, hc_t4c4_t3c3$theta, paired = TRUE) # C > D
