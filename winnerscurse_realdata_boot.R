## WINNER'S CURSE REAL DATA - BOOT: 

## This script uses 6 sets of GWAS summary statistics, two related to each
## trait. The traits of interest are body mass index (BMI), type 2 diabetes
## (T2D) and height. Each data set contains six columns, named 'chr', 'pos',
## 'rsid', 'beta', 'se' and 'beta_rep'. Here 'beta_rep' contains the association
## estimate obtained for each SNP in the corresponding replication GWAS. This
## script evaluates a number of winner's curse correction methods by comparing
## their performance at two significance thresholds, 5e-8 and 5e-4, for each
## data set, using the estimated MSE among significant SNPs.


## Load required packages:
library(winnerscurse)
library(ggplot2)
library(dplyr)
library(scam)
library(tidyr)
library("RColorBrewer")
col <- brewer.pal(8,"Dark2")
col1 <- brewer.pal(11,"RdYlBu")



## alternative bootstrap function 
BR_ss_2 <- function(summary_data,seed_opt = FALSE,seed=1998){
  stopifnot(all(c("rsid", "beta","se") %in% names(summary_data)))
  stopifnot(!all(is.na(summary_data$rsid)) && !all(is.na(summary_data$beta)) && !all(is.na(summary_data$se)))
  stopifnot(nrow(summary_data) > 5)
  stopifnot(is.numeric(summary_data$beta) && is.numeric(summary_data$se))
  stopifnot(!any(duplicated(summary_data$rsid)))
  summary_data <- dplyr::arrange(summary_data, dplyr::desc((summary_data$beta/summary_data$se)))
  N <- nrow(summary_data)
  if(seed_opt==TRUE){set.seed(seed)}
  beta_boot <- matrix(stats::rnorm(1*N, mean = rep(summary_data$beta,1), sd = rep(summary_data$se,1)), nrow=N, ncol=1, byrow=FALSE)
  beta_mat <- matrix(rep(summary_data$beta,1), nrow=N, ncol=1, byrow=FALSE)
  se_mat <- matrix(rep(summary_data$se,1), nrow=N, ncol=1, byrow=FALSE)
  beta_oob <- beta_mat
  ordering <- apply(beta_boot/se_mat, 2,order,decreasing=TRUE)
  bias_correct <- matrix(nrow=N, ncol=1)
  bias_correct[,1] <- (beta_boot[ordering[,1],1] - beta_oob[ordering[,1],1])/summary_data$se[ordering[,1]]
  z <- summary_data$beta/summary_data$se
  bias_correct <- stats::predict(stats::smooth.spline(z,bias_correct)$fit, z)$y
  beta_BR_ss <- summary_data$beta - summary_data$se*bias_correct[rank(-1*summary_data$beta/summary_data$se)]
  beta_BR_ss[sign(beta_BR_ss) != sign(summary_data$beta)] <- 0
  summary_data <- cbind(summary_data, beta_BR_ss)
  ## following line removed in this version!
  ## for (i in 1:N){
  ##  if(abs(summary_data$beta[i]) < abs(summary_data$beta_BR_ss[i])){summary_data$beta_BR_ss[i] <- summary_data$beta[i]}
  ## }
  summary_data <- dplyr::arrange(summary_data,dplyr::desc(abs(summary_data$beta/summary_data$se)))
  return(summary_data)
}


## PART 1) COMPUTING ESTIMATED MSE FOR EACH DATASET, METHOD AND THRESHOLD

## BMI 1

summary_data_bmi_1 <-  read.table('data/summary_data_bmi_1.txt',header=TRUE)

# Apply methods
summary_stats <- summary_data_bmi_1[,3:5]
out_BR1 <- BR_ss(summary_data = summary_stats)
out_BR2 <- BR_ss_2(summary_data = summary_stats)

# Use threshold 5e-8
out_BR1_ss_sig <- out_BR1[2*(pnorm(abs(out_BR1$beta/out_BR1$se), lower.tail=FALSE)) < 5e-8,]
out_BR2_ss_sig <- out_BR2[2*(pnorm(abs(out_BR2$beta/out_BR2$se), lower.tail=FALSE)) < 5e-8,]

summary_data_bmi_1 <- dplyr::arrange(summary_data_bmi_1,dplyr::desc(abs(summary_data_bmi_1$beta/summary_data_bmi_1$se)))
summary_data_sig_bmi_1 <- summary_data_bmi_1[2*(pnorm(abs(summary_data_bmi_1$beta/summary_data_bmi_1$se), lower.tail=FALSE)) < 5e-8,]
true_beta <- summary_data_sig_bmi_1$beta_rep

# Compute estimated MSE for each method for 5e-8 threshold
mse_BMI1 <- data.frame(naive = mean((true_beta - out_BR1_ss_sig$beta)^2) - mean(out_BR1_ss_sig$se^2), boot1 = mean((true_beta - out_BR1_ss_sig$beta_BR_ss)^2) - mean(out_BR1_ss_sig$se^2), boot2 = mean((true_beta - out_BR2_ss_sig$beta_BR_ss)^2) - mean(out_BR2_ss_sig$se^2), prop_nochange = sum(out_BR1_ss_sig$beta==out_BR1_ss_sig$beta_BR_ss)/length(out_BR1_ss_sig$rsid))

# Use threshold 5e-4 
out_BR1_ss_sig <- out_BR1[2*(pnorm(abs(out_BR1$beta/out_BR1$se), lower.tail=FALSE)) < 5e-4,]
out_BR2_ss_sig <- out_BR2[2*(pnorm(abs(out_BR2$beta/out_BR2$se), lower.tail=FALSE)) < 5e-4,]
summary_data_sig2 <- summary_data_bmi_1[2*(pnorm(abs(summary_data_bmi_1$beta/summary_data_bmi_1$se), lower.tail=FALSE)) < 5e-4,]
true_beta2 <- summary_data_sig2$beta_rep

# Compute estimated MSE for each method for 5e-4 threshold
mse_BMI1_5e_4 <- data.frame(naive = mean((true_beta2 - out_BR1_ss_sig$beta)^2) - mean(out_BR1_ss_sig$se^2), boot1 = mean((true_beta2 - out_BR1_ss_sig$beta_BR_ss)^2) - mean(out_BR1_ss_sig$se^2), boot2 = mean((true_beta2 - out_BR2_ss_sig$beta_BR_ss)^2) - mean(out_BR2_ss_sig$se^2), prop_nochange = sum(out_BR1_ss_sig$beta==out_BR1_ss_sig$beta_BR_ss)/length(out_BR1_ss_sig$rsid))

###############################################################################

## BMI 2

summary_data_bmi_2 <-  read.table('data/summary_data_bmi_2.txt',header=TRUE)

# Apply methods
summary_stats <- summary_data_bmi_2[,3:5]
out_BR1 <- BR_ss(summary_data = summary_stats)
out_BR2 <- BR_ss_2(summary_data = summary_stats)

# Use threshold 5e-8
out_BR1_ss_sig <- out_BR1[2*(pnorm(abs(out_BR1$beta/out_BR1$se), lower.tail=FALSE)) < 5e-8,]
out_BR2_ss_sig <- out_BR2[2*(pnorm(abs(out_BR2$beta/out_BR2$se), lower.tail=FALSE)) < 5e-8,]

summary_data_bmi_2 <- dplyr::arrange(summary_data_bmi_2,dplyr::desc(abs(summary_data_bmi_2$beta/summary_data_bmi_2$se)))
summary_data_sig_bmi_2 <- summary_data_bmi_1[2*(pnorm(abs(summary_data_bmi_2$beta/summary_data_bmi_2$se), lower.tail=FALSE)) < 5e-8,]
true_beta <- summary_data_sig_bmi_2$beta_rep

# Compute estimated MSE for each method for 5e-8 threshold
mse_BMI2 <- data.frame(naive = mean((true_beta - out_BR1_ss_sig$beta)^2) - mean(out_BR1_ss_sig$se^2), boot1 = mean((true_beta - out_BR1_ss_sig$beta_BR_ss)^2) - mean(out_BR1_ss_sig$se^2), boot2 = mean((true_beta - out_BR2_ss_sig$beta_BR_ss)^2) - mean(out_BR2_ss_sig$se^2), prop_nochange = sum(out_BR1_ss_sig$beta==out_BR1_ss_sig$beta_BR_ss)/length(out_BR1_ss_sig$rsid))

# Use threshold 5e-4 
out_BR1_ss_sig <- out_BR1[2*(pnorm(abs(out_BR1$beta/out_BR1$se), lower.tail=FALSE)) < 5e-4,]
out_BR2_ss_sig <- out_BR2[2*(pnorm(abs(out_BR2$beta/out_BR2$se), lower.tail=FALSE)) < 5e-4,]
summary_data_sig2 <- summary_data_bmi_2[2*(pnorm(abs(summary_data_bmi_2$beta/summary_data_bmi_2$se), lower.tail=FALSE)) < 5e-4,]
true_beta2 <- summary_data_sig2$beta_rep

# Compute estimated MSE for each method for 5e-4 threshold
mse_BMI2_5e_4 <- data.frame(naive = mean((true_beta2 - out_BR1_ss_sig$beta)^2) - mean(out_BR1_ss_sig$se^2), boot1 = mean((true_beta2 - out_BR1_ss_sig$beta_BR_ss)^2) - mean(out_BR1_ss_sig$se^2), boot2 = mean((true_beta2 - out_BR2_ss_sig$beta_BR_ss)^2) - mean(out_BR2_ss_sig$se^2), prop_nochange = sum(out_BR1_ss_sig$beta==out_BR1_ss_sig$beta_BR_ss)/length(out_BR1_ss_sig$rsid))

###############################################################################

## T2D 1

summary_data_T2D_1 <-  read.table('data/summary_data_T2D_1.txt',header=TRUE)
summary_data_T2D_1$beta <- log(summary_data_T2D_1$beta)
summary_data_T2D_1$beta_rep <- log(summary_data_T2D_1$beta_rep)

# Apply methods
summary_stats <- summary_data_T2D_1[,3:5]
out_BR1 <- BR_ss(summary_data = summary_stats)
out_BR2 <- BR_ss_2(summary_data = summary_stats)

# Use threshold 5e-8
out_BR1_ss_sig <- out_BR1[2*(pnorm(abs(out_BR1$beta/out_BR1$se), lower.tail=FALSE)) < 5e-8,]
out_BR2_ss_sig <- out_BR2[2*(pnorm(abs(out_BR2$beta/out_BR2$se), lower.tail=FALSE)) < 5e-8,]

summary_data_T2D_1 <- dplyr::arrange(summary_data_T2D_1,dplyr::desc(abs(summary_data_T2D_1$beta/summary_data_T2D_1$se)))
summary_data_sig_T2D_1 <- summary_data_T2D_1[2*(pnorm(abs(summary_data_T2D_1$beta/summary_data_T2D_1$se), lower.tail=FALSE)) < 5e-8,]
true_beta <- summary_data_sig_T2D_1$beta_rep

# Compute estimated MSE for each method for 5e-8 threshold
mse_T2D1 <- data.frame(naive = mean((true_beta - out_BR1_ss_sig$beta)^2) - mean(out_BR1_ss_sig$se^2), boot1 = mean((true_beta - out_BR1_ss_sig$beta_BR_ss)^2) - mean(out_BR1_ss_sig$se^2), boot2 = mean((true_beta - out_BR2_ss_sig$beta_BR_ss)^2) - mean(out_BR2_ss_sig$se^2), prop_nochange = sum(out_BR1_ss_sig$beta==out_BR1_ss_sig$beta_BR_ss)/length(out_BR1_ss_sig$rsid))

# Use threshold 5e-4 
out_BR1_ss_sig <- out_BR1[2*(pnorm(abs(out_BR1$beta/out_BR1$se), lower.tail=FALSE)) < 5e-4,]
out_BR2_ss_sig <- out_BR2[2*(pnorm(abs(out_BR2$beta/out_BR2$se), lower.tail=FALSE)) < 5e-4,]
summary_data_sig2 <- summary_data_T2D_1[2*(pnorm(abs(summary_data_T2D_1$beta/summary_data_T2D_1$se), lower.tail=FALSE)) < 5e-4,]
true_beta2 <- summary_data_sig2$beta_rep

# Compute estimated MSE for each method for 5e-4 threshold
mse_T2D1_5e_4 <- data.frame(naive = mean((true_beta2 - out_BR1_ss_sig$beta)^2) - mean(out_BR1_ss_sig$se^2), boot1 = mean((true_beta2 - out_BR1_ss_sig$beta_BR_ss)^2) - mean(out_BR1_ss_sig$se^2), boot2 = mean((true_beta2 - out_BR2_ss_sig$beta_BR_ss)^2) - mean(out_BR2_ss_sig$se^2), prop_nochange = sum(out_BR1_ss_sig$beta==out_BR1_ss_sig$beta_BR_ss)/length(out_BR1_ss_sig$rsid))

###############################################################################

## T2D 2

summary_data_T2D_2 <-  read.table('data/summary_data_T2D_2.txt',header=TRUE)
summary_data_T2D_2$beta <- log(summary_data_T2D_2$beta)
summary_data_T2D_2$beta_rep <- log(summary_data_T2D_2$beta_rep)

# Apply methods
summary_stats <- summary_data_T2D_2[,3:5]
out_BR1 <- BR_ss(summary_data = summary_stats)
out_BR2 <- BR_ss_2(summary_data = summary_stats)

# Use threshold 5e-8
out_BR1_ss_sig <- out_BR1[2*(pnorm(abs(out_BR1$beta/out_BR1$se), lower.tail=FALSE)) < 5e-8,]
out_BR2_ss_sig <- out_BR2[2*(pnorm(abs(out_BR2$beta/out_BR2$se), lower.tail=FALSE)) < 5e-8,]

summary_data_T2D_2 <- dplyr::arrange(summary_data_T2D_2,dplyr::desc(abs(summary_data_T2D_2$beta/summary_data_T2D_1$se)))
summary_data_sig_T2D_2 <- summary_data_T2D_2[2*(pnorm(abs(summary_data_T2D_2$beta/summary_data_T2D_2$se), lower.tail=FALSE)) < 5e-8,]
true_beta <- summary_data_sig_T2D_2$beta_rep

# Compute estimated MSE for each method for 5e-8 threshold
mse_T2D2 <- data.frame(naive = mean((true_beta - out_BR1_ss_sig$beta)^2) - mean(out_BR1_ss_sig$se^2), boot1 = mean((true_beta - out_BR1_ss_sig$beta_BR_ss)^2) - mean(out_BR1_ss_sig$se^2), boot2 = mean((true_beta - out_BR2_ss_sig$beta_BR_ss)^2) - mean(out_BR2_ss_sig$se^2), prop_nochange = sum(out_BR1_ss_sig$beta==out_BR1_ss_sig$beta_BR_ss)/length(out_BR1_ss_sig$rsid))

# Use threshold 5e-4 
out_BR1_ss_sig <- out_BR1[2*(pnorm(abs(out_BR1$beta/out_BR1$se), lower.tail=FALSE)) < 5e-4,]
out_BR2_ss_sig <- out_BR2[2*(pnorm(abs(out_BR2$beta/out_BR2$se), lower.tail=FALSE)) < 5e-4,]
summary_data_sig2 <- summary_data_T2D_2[2*(pnorm(abs(summary_data_T2D_2$beta/summary_data_T2D_2$se), lower.tail=FALSE)) < 5e-4,]
true_beta2 <- summary_data_sig2$beta_rep

# Compute estimated MSE for each method for 5e-4 threshold
mse_T2D2_5e_4 <- data.frame(naive = mean((true_beta2 - out_BR1_ss_sig$beta)^2) - mean(out_BR1_ss_sig$se^2), boot1 = mean((true_beta2 - out_BR1_ss_sig$beta_BR_ss)^2) - mean(out_BR1_ss_sig$se^2), boot2 = mean((true_beta2 - out_BR2_ss_sig$beta_BR_ss)^2) - mean(out_BR2_ss_sig$se^2), prop_nochange = sum(out_BR1_ss_sig$beta==out_BR1_ss_sig$beta_BR_ss)/length(out_BR1_ss_sig$rsid))

###############################################################################

## Height 1

summary_data_height_1 <-  read.table('data/summary_data_height_1.txt',header=TRUE)

# Apply methods
summary_stats <- summary_data_height_1[,3:5]
out_BR1 <- BR_ss(summary_data = summary_stats)
out_BR2 <- BR_ss_2(summary_data = summary_stats)

# Use threshold 5e-8
out_BR1_ss_sig <- out_BR1[2*(pnorm(abs(out_BR1$beta/out_BR1$se), lower.tail=FALSE)) < 5e-8,]
out_BR2_ss_sig <- out_BR2[2*(pnorm(abs(out_BR2$beta/out_BR2$se), lower.tail=FALSE)) < 5e-8,]

summary_data_height_1 <- dplyr::arrange(summary_data_height_1,dplyr::desc(abs(summary_data_height_1$beta/summary_data_height_1$se)))
summary_data_sig_height_1 <- summary_data_height_1[2*(pnorm(abs(summary_data_height_1$beta/summary_data_height_1$se), lower.tail=FALSE)) < 5e-8,]
true_beta <- summary_data_sig_height_1$beta_rep

# Compute estimated MSE for each method for 5e-8 threshold
mse_height1 <- data.frame(naive = mean((true_beta - out_BR1_ss_sig$beta)^2) - mean(out_BR1_ss_sig$se^2), boot1 = mean((true_beta - out_BR1_ss_sig$beta_BR_ss)^2) - mean(out_BR1_ss_sig$se^2), boot2 = mean((true_beta - out_BR2_ss_sig$beta_BR_ss)^2) - mean(out_BR2_ss_sig$se^2), prop_nochange = sum(out_BR1_ss_sig$beta==out_BR1_ss_sig$beta_BR_ss)/length(out_BR1_ss_sig$rsid))

# Use threshold 5e-4 
out_BR1_ss_sig <- out_BR1[2*(pnorm(abs(out_BR1$beta/out_BR1$se), lower.tail=FALSE)) < 5e-4,]
out_BR2_ss_sig <- out_BR2[2*(pnorm(abs(out_BR2$beta/out_BR2$se), lower.tail=FALSE)) < 5e-4,]
summary_data_sig2 <- summary_data_height_1[2*(pnorm(abs(summary_data_height_1$beta/summary_data_height_1$se), lower.tail=FALSE)) < 5e-4,]
true_beta2 <- summary_data_sig2$beta_rep

# Compute estimated MSE for each method for 5e-4 threshold
mse_height1_5e_4 <- data.frame(naive = mean((true_beta2 - out_BR1_ss_sig$beta)^2) - mean(out_BR1_ss_sig$se^2), boot1 = mean((true_beta2 - out_BR1_ss_sig$beta_BR_ss)^2) - mean(out_BR1_ss_sig$se^2), boot2 = mean((true_beta2 - out_BR2_ss_sig$beta_BR_ss)^2) - mean(out_BR2_ss_sig$se^2), prop_nochange = sum(out_BR1_ss_sig$beta==out_BR1_ss_sig$beta_BR_ss)/length(out_BR1_ss_sig$rsid))

###############################################################################

## Height 2

summary_data_height_2 <-  read.table('data/summary_data_height_2.txt',header=TRUE)

# Apply methods
summary_stats <- summary_data_height_2[,3:5]
out_BR1 <- BR_ss(summary_data = summary_stats)
out_BR2 <- BR_ss_2(summary_data = summary_stats)

# Use threshold 5e-8
out_BR1_ss_sig <- out_BR1[2*(pnorm(abs(out_BR1$beta/out_BR1$se), lower.tail=FALSE)) < 5e-8,]
out_BR2_ss_sig <- out_BR2[2*(pnorm(abs(out_BR2$beta/out_BR2$se), lower.tail=FALSE)) < 5e-8,]

summary_data_height_2 <- dplyr::arrange(summary_data_height_2,dplyr::desc(abs(summary_data_height_2$beta/summary_data_height_2$se)))
summary_data_sig_height_2 <- summary_data_height_2[2*(pnorm(abs(summary_data_height_2$beta/summary_data_height_2$se), lower.tail=FALSE)) < 5e-8,]
true_beta <- summary_data_sig_height_2$beta_rep

# Compute estimated MSE for each method for 5e-8 threshold
mse_height2 <- data.frame(naive = mean((true_beta - out_BR1_ss_sig$beta)^2) - mean(out_BR1_ss_sig$se^2), boot1 = mean((true_beta - out_BR1_ss_sig$beta_BR_ss)^2) - mean(out_BR1_ss_sig$se^2), boot2 = mean((true_beta - out_BR2_ss_sig$beta_BR_ss)^2) - mean(out_BR2_ss_sig$se^2), prop_nochange = sum(out_BR1_ss_sig$beta==out_BR1_ss_sig$beta_BR_ss)/length(out_BR1_ss_sig$rsid))

# Use threshold 5e-4 
out_BR1_ss_sig <- out_BR1[2*(pnorm(abs(out_BR1$beta/out_BR1$se), lower.tail=FALSE)) < 5e-4,]
out_BR2_ss_sig <- out_BR2[2*(pnorm(abs(out_BR2$beta/out_BR2$se), lower.tail=FALSE)) < 5e-4,]
summary_data_sig2 <- summary_data_height_2[2*(pnorm(abs(summary_data_height_2$beta/summary_data_height_2$se), lower.tail=FALSE)) < 5e-4,]
true_beta2 <- summary_data_sig2$beta_rep

# Compute estimated MSE for each method for 5e-4 threshold
mse_height2_5e_4 <- data.frame(naive = mean((true_beta2 - out_BR1_ss_sig$beta)^2) - mean(out_BR1_ss_sig$se^2), boot1 = mean((true_beta2 - out_BR1_ss_sig$beta_BR_ss)^2) - mean(out_BR1_ss_sig$se^2), boot2 = mean((true_beta2 - out_BR2_ss_sig$beta_BR_ss)^2) - mean(out_BR2_ss_sig$se^2), prop_nochange = sum(out_BR1_ss_sig$beta==out_BR1_ss_sig$beta_BR_ss)/length(out_BR1_ss_sig$rsid))

###############################################################################

## Combine all above results and save them in external results files 

mse_5e_8 <- rbind(mse_BMI1,mse_BMI2, mse_T2D1, mse_T2D2, mse_height1, mse_height2)
mse_5e_8 <- round(mse_5e_8[,1:11],5)
GWAS <- c("BMI 1", "BMI 2", "T2D 1", "T2D 2", "Height 1", "Height 2")
mse_5e_8 <- cbind(GWAS,mse_5e_8)
write.csv(mse_5e_8,"results/mse_5e_8_boot.txt", row.names = FALSE)

mse_5e_4 <- rbind(mse_BMI1_5e_4,mse_BMI2_5e_4, mse_T2D1_5e_4, mse_T2D2_5e_4, mse_height1_5e_4, mse_height2_5e_4)
mse_5e_4 <- round(mse_5e_4[,1:11],5)
GWAS <- c("BMI 1", "BMI 2", "T2D 1", "T2D 2", "Height 1", "Height 2")
mse_5e_4 <- cbind(GWAS,mse_5e_4)
write.csv(mse_5e_4,"results/mse_5e_4_boot.txt", row.names = FALSE)