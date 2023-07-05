## WINNER'S CURSE REAL DATA: 

## This script uses 6 sets of GWAS summary statistics, two related to each
## trait. The traits of interest are body mass index (T2D), type 2 diabetes
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


## PART 1) COMPUTING ESTIMATED MSE FOR EACH DATASET, METHOD AND THRESHOLD

## BMI 1

summary_data_bmi_1 <-  read.table('data/summary_data_bmi_1.txt',header=TRUE)

# Apply methods
summary_stats <- summary_data_bmi_1[,3:5]
out_CL <- conditional_likelihood(summary_data = summary_stats)
out_EB <- empirical_bayes(summary_data = summary_stats)
out_EB_df <- empirical_bayes(summary_data=summary_stats, method="fix_df")
out_EB_scam <- empirical_bayes(summary_data=summary_stats,method="scam")
out_EB_gam_po <- empirical_bayes(summary_data=summary_stats,method="gam_po")
out_EB_gam_nb <- empirical_bayes(summary_data=summary_stats,method="gam_nb")
out_FIQT <- FDR_IQT(summary_data = summary_stats)
out_BR <- BR_ss(summary_data = summary_stats)

# Use threshold 5e-8
out_EB_sig <- out_EB[2*(pnorm(abs(out_EB$beta/out_EB$se), lower.tail=FALSE)) < 5e-8,]
out_EB_df_sig <- out_EB_df[2*(pnorm(abs(out_EB_df$beta/out_EB_df$se), lower.tail=FALSE)) < 5e-8,]
out_EB_scam_sig <- out_EB_scam[2*(pnorm(abs(out_EB_scam$beta/out_EB_scam$se), lower.tail=FALSE)) < 5e-8,]
out_EB_gampo_sig <- out_EB_gam_po[2*(pnorm(abs(out_EB_gam_po$beta/out_EB_gam_po$se), lower.tail=FALSE)) < 5e-8,]
out_EB_gamnb_sig <- out_EB_gam_nb[2*(pnorm(abs(out_EB_gam_nb$beta/out_EB_gam_nb$se), lower.tail=FALSE)) < 5e-8,]
out_FIQT_sig <- out_FIQT[2*(pnorm(abs(out_FIQT$beta/out_FIQT$se), lower.tail=FALSE)) < 5e-8,]
out_BR_ss_sig <- out_BR[2*(pnorm(abs(out_BR$beta/out_BR$se), lower.tail=FALSE)) < 5e-8,]

summary_data_bmi_1 <- dplyr::arrange(summary_data_bmi_1,dplyr::desc(abs(summary_data_bmi_1$beta/summary_data_bmi_1$se)))
summary_data_sig_bmi_1 <- summary_data_bmi_1[2*(pnorm(abs(summary_data_bmi_1$beta/summary_data_bmi_1$se), lower.tail=FALSE)) < 5e-8,]
true_beta <- summary_data_sig_bmi_1$beta_rep

# Compute estimated MSE for each method for 5e-8 threshold
mse_BMI1 <- data.frame(naive = mean((true_beta - out_EB_sig$beta)^2) - mean(out_EB_sig$se^2), CL1 = mean((true_beta - out_CL$beta.cl1)^2) - mean(out_EB_sig$se^2), CL2 = mean((true_beta - out_CL$beta.cl2)^2) - mean(out_EB_sig$se^2), CL3 = mean((true_beta - out_CL$beta.cl3)^2) - mean(out_EB_sig$se^2), EB = mean((true_beta - out_EB_sig$beta_EB)^2) - mean(out_EB_sig$se^2),  EB_df = mean((true_beta - out_EB_df_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_scam = mean((true_beta - out_EB_scam_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_po = mean((true_beta - out_EB_gampo_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_nb = mean((true_beta - out_EB_gamnb_sig$beta_EB)^2) - mean(out_EB_sig$se^2), boot = mean((true_beta - out_BR_ss_sig$beta_BR_ss)^2) - mean(out_EB_sig$se^2), FIQT = mean((true_beta - out_FIQT_sig$beta_FIQT)^2) - mean(out_EB_sig$se^2))

# Use threshold 5e-4 (need to re-apply conditional likelihood with change of threshold)
out_EB_sig <- out_EB[2*(pnorm(abs(out_EB$beta/out_EB$se), lower.tail=FALSE)) < 5e-4,]
out_EB_df_sig <- out_EB_df[2*(pnorm(abs(out_EB_df$beta/out_EB_df$se), lower.tail=FALSE)) < 5e-4,]
out_EB_scam_sig <- out_EB_scam[2*(pnorm(abs(out_EB_scam$beta/out_EB_scam$se), lower.tail=FALSE)) < 5e-4,]
out_EB_gampo_sig <- out_EB_gam_po[2*(pnorm(abs(out_EB_gam_po$beta/out_EB_gam_po$se), lower.tail=FALSE)) < 5e-4,]
out_EB_gamnb_sig <- out_EB_gam_nb[2*(pnorm(abs(out_EB_gam_nb$beta/out_EB_gam_nb$se), lower.tail=FALSE)) < 5e-4,]
out_FIQT_sig <- out_FIQT[2*(pnorm(abs(out_FIQT$beta/out_FIQT$se), lower.tail=FALSE)) < 5e-4,]
out_BR_ss_sig <- out_BR[2*(pnorm(abs(out_BR$beta/out_BR$se), lower.tail=FALSE)) < 5e-4,]

out_CL <- conditional_likelihood(summary_data = summary_stats, alpha=5e-4)

summary_data_sig2 <- summary_data_bmi_1[2*(pnorm(abs(summary_data_bmi_1$beta/summary_data_bmi_1$se), lower.tail=FALSE)) < 5e-4,]
true_beta2 <- summary_data_sig2$beta_rep

# Compute estimated MSE for each method for 5e-4 threshold
mse_BMI1_5e_4 <- data.frame(naive = mean((true_beta2 - out_EB_sig$beta)^2) - mean(out_EB_sig$se^2), CL1 = mean((true_beta2 - out_CL$beta.cl1)^2) - mean(out_EB_sig$se^2), CL2 = mean((true_beta2 - out_CL$beta.cl2)^2) - mean(out_EB_sig$se^2), CL3 = mean((true_beta2 - out_CL$beta.cl3)^2) - mean(out_EB_sig$se^2), EB = mean((true_beta2 - out_EB_sig$beta_EB)^2) - mean(out_EB_sig$se^2),  EB_df = mean((true_beta2 - out_EB_df_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_scam = mean((true_beta2 - out_EB_scam_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_po = mean((true_beta2 - out_EB_gampo_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_nb = mean((true_beta2 - out_EB_gamnb_sig$beta_EB)^2) - mean(out_EB_sig$se^2), boot = mean((true_beta2 - out_BR_ss_sig$beta_BR_ss)^2) - mean(out_EB_sig$se^2), FIQT = mean((true_beta2 - out_FIQT_sig$beta_FIQT)^2) - mean(out_EB_sig$se^2))

###############################################################################

## BMI 2

summary_data_bmi_2 <-  read.table('data/summary_data_bmi_2.txt',header=TRUE)

# Apply methods
summary_stats <- summary_data_bmi_2[,3:5]
out_CL <- conditional_likelihood(summary_data = summary_stats)
out_EB <- empirical_bayes(summary_data = summary_stats)
out_EB_df <- empirical_bayes(summary_data=summary_stats, method="fix_df")
out_EB_scam <- empirical_bayes(summary_data=summary_stats,method="scam")
out_EB_gam_po <- empirical_bayes(summary_data=summary_stats,method="gam_po")
out_EB_gam_nb <- empirical_bayes(summary_data=summary_stats,method="gam_nb")
out_FIQT <- FDR_IQT(summary_data = summary_stats)
out_BR <- BR_ss(summary_data = summary_stats)

# Use threshold 5e-8
out_EB_sig <- out_EB[2*(pnorm(abs(out_EB$beta/out_EB$se), lower.tail=FALSE)) < 5e-8,]
out_EB_df_sig <- out_EB_df[2*(pnorm(abs(out_EB_df$beta/out_EB_df$se), lower.tail=FALSE)) < 5e-8,]
out_EB_scam_sig <- out_EB_scam[2*(pnorm(abs(out_EB_scam$beta/out_EB_scam$se), lower.tail=FALSE)) < 5e-8,]
out_EB_gampo_sig <- out_EB_gam_po[2*(pnorm(abs(out_EB_gam_po$beta/out_EB_gam_po$se), lower.tail=FALSE)) < 5e-8,]
out_EB_gamnb_sig <- out_EB_gam_nb[2*(pnorm(abs(out_EB_gam_nb$beta/out_EB_gam_nb$se), lower.tail=FALSE)) < 5e-8,]
out_FIQT_sig <- out_FIQT[2*(pnorm(abs(out_FIQT$beta/out_FIQT$se), lower.tail=FALSE)) < 5e-8,]
out_BR_ss_sig <- out_BR[2*(pnorm(abs(out_BR$beta/out_BR$se), lower.tail=FALSE)) < 5e-8,]

summary_data_bmi_2 <- dplyr::arrange(summary_data_bmi_2,dplyr::desc(abs(summary_data_bmi_2$beta/summary_data_bmi_2$se)))
summary_data_sig_bmi_2 <- summary_data_bmi_2[2*(pnorm(abs(summary_data_bmi_2$beta/summary_data_bmi_2$se), lower.tail=FALSE)) < 5e-8,]
true_beta <- summary_data_sig_bmi_2$beta_rep

# Compute estimated MSE for each method for 5e-8 threshold
mse_BMI2 <- data.frame(naive = mean((true_beta - out_EB_sig$beta)^2) - mean(out_EB_sig$se^2), CL1 = mean((true_beta - out_CL$beta.cl1)^2) - mean(out_EB_sig$se^2), CL2 = mean((true_beta - out_CL$beta.cl2)^2) - mean(out_EB_sig$se^2), CL3 = mean((true_beta - out_CL$beta.cl3)^2) - mean(out_EB_sig$se^2), EB = mean((true_beta - out_EB_sig$beta_EB)^2) - mean(out_EB_sig$se^2),  EB_df = mean((true_beta - out_EB_df_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_scam = mean((true_beta - out_EB_scam_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_po = mean((true_beta - out_EB_gampo_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_nb = mean((true_beta - out_EB_gamnb_sig$beta_EB)^2) - mean(out_EB_sig$se^2), boot = mean((true_beta - out_BR_ss_sig$beta_BR_ss)^2) - mean(out_EB_sig$se^2), FIQT = mean((true_beta - out_FIQT_sig$beta_FIQT)^2) - mean(out_EB_sig$se^2))

# Use threshold 5e-4 (need to re-apply conditional likelihood with change of threshold)
out_EB_sig <- out_EB[2*(pnorm(abs(out_EB$beta/out_EB$se), lower.tail=FALSE)) < 5e-4,]
out_EB_df_sig <- out_EB_df[2*(pnorm(abs(out_EB_df$beta/out_EB_df$se), lower.tail=FALSE)) < 5e-4,]
out_EB_scam_sig <- out_EB_scam[2*(pnorm(abs(out_EB_scam$beta/out_EB_scam$se), lower.tail=FALSE)) < 5e-4,]
out_EB_gampo_sig <- out_EB_gam_po[2*(pnorm(abs(out_EB_gam_po$beta/out_EB_gam_po$se), lower.tail=FALSE)) < 5e-4,]
out_EB_gamnb_sig <- out_EB_gam_nb[2*(pnorm(abs(out_EB_gam_nb$beta/out_EB_gam_nb$se), lower.tail=FALSE)) < 5e-4,]
out_FIQT_sig <- out_FIQT[2*(pnorm(abs(out_FIQT$beta/out_FIQT$se), lower.tail=FALSE)) < 5e-4,]
out_BR_ss_sig <- out_BR[2*(pnorm(abs(out_BR$beta/out_BR$se), lower.tail=FALSE)) < 5e-4,]

out_CL <- conditional_likelihood(summary_data = summary_stats, alpha=5e-4)

summary_data_sig2 <- summary_data_bmi_2[2*(pnorm(abs(summary_data_bmi_2$beta/summary_data_bmi_2$se), lower.tail=FALSE)) < 5e-4,]
true_beta2 <- summary_data_sig2$beta_rep

# Compute estimated MSE for each method for 5e-4 threshold
mse_BMI2_5e_4 <- data.frame(naive = mean((true_beta2 - out_EB_sig$beta)^2) - mean(out_EB_sig$se^2), CL1 = mean((true_beta2 - out_CL$beta.cl1)^2) - mean(out_EB_sig$se^2), CL2 = mean((true_beta2 - out_CL$beta.cl2)^2) - mean(out_EB_sig$se^2), CL3 = mean((true_beta2 - out_CL$beta.cl3)^2) - mean(out_EB_sig$se^2), EB = mean((true_beta2 - out_EB_sig$beta_EB)^2) - mean(out_EB_sig$se^2),  EB_df = mean((true_beta2 - out_EB_df_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_scam = mean((true_beta2 - out_EB_scam_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_po = mean((true_beta2 - out_EB_gampo_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_nb = mean((true_beta2 - out_EB_gamnb_sig$beta_EB)^2) - mean(out_EB_sig$se^2), boot = mean((true_beta2 - out_BR_ss_sig$beta_BR_ss)^2) - mean(out_EB_sig$se^2), FIQT = mean((true_beta2 - out_FIQT_sig$beta_FIQT)^2) - mean(out_EB_sig$se^2))

###############################################################################

## T2D 1

summary_data_T2D_1 <-  read.table('data/summary_data_T2D_1.txt',header=TRUE)
summary_data_T2D_1$beta <- log(summary_data_T2D_1$beta)
summary_data_T2D_1$beta_rep <- log(summary_data_T2D_1$beta_rep)

# Apply methods
summary_stats <- summary_data_T2D_1[,3:5]
out_CL <- conditional_likelihood(summary_data = summary_stats)
out_EB <- empirical_bayes(summary_data = summary_stats)
out_EB_df <- empirical_bayes(summary_data=summary_stats, method="fix_df")
out_EB_scam <- empirical_bayes(summary_data=summary_stats,method="scam")
out_EB_gam_po <- empirical_bayes(summary_data=summary_stats,method="gam_po")
out_EB_gam_nb <- empirical_bayes(summary_data=summary_stats,method="gam_nb")
out_FIQT <- FDR_IQT(summary_data = summary_stats)
out_BR <- BR_ss(summary_data = summary_stats)

# Use threshold 5e-8
out_EB_sig <- out_EB[2*(pnorm(abs(out_EB$beta/out_EB$se), lower.tail=FALSE)) < 5e-8,]
out_EB_df_sig <- out_EB_df[2*(pnorm(abs(out_EB_df$beta/out_EB_df$se), lower.tail=FALSE)) < 5e-8,]
out_EB_scam_sig <- out_EB_scam[2*(pnorm(abs(out_EB_scam$beta/out_EB_scam$se), lower.tail=FALSE)) < 5e-8,]
out_EB_gampo_sig <- out_EB_gam_po[2*(pnorm(abs(out_EB_gam_po$beta/out_EB_gam_po$se), lower.tail=FALSE)) < 5e-8,]
out_EB_gamnb_sig <- out_EB_gam_nb[2*(pnorm(abs(out_EB_gam_nb$beta/out_EB_gam_nb$se), lower.tail=FALSE)) < 5e-8,]
out_FIQT_sig <- out_FIQT[2*(pnorm(abs(out_FIQT$beta/out_FIQT$se), lower.tail=FALSE)) < 5e-8,]
out_BR_ss_sig <- out_BR[2*(pnorm(abs(out_BR$beta/out_BR$se), lower.tail=FALSE)) < 5e-8,]

summary_data_T2D_1 <- dplyr::arrange(summary_data_T2D_1,dplyr::desc(abs(summary_data_T2D_1$beta/summary_data_T2D_1$se)))
summary_data_sig_T2D_1 <- summary_data_T2D_1[2*(pnorm(abs(summary_data_T2D_1$beta/summary_data_T2D_1$se), lower.tail=FALSE)) < 5e-8,]
true_beta <- summary_data_sig_T2D_1$beta_rep

# Compute estimated MSE for each method at 5e-8 threshold
mse_T2D1 <- data.frame(naive = mean((true_beta - out_EB_sig$beta)^2) - mean(out_EB_sig$se^2), CL1 = mean((true_beta - out_CL$beta.cl1)^2) - mean(out_EB_sig$se^2), CL2 = mean((true_beta - out_CL$beta.cl2)^2) - mean(out_EB_sig$se^2), CL3 = mean((true_beta - out_CL$beta.cl3)^2) - mean(out_EB_sig$se^2), EB = mean((true_beta - out_EB_sig$beta_EB)^2) - mean(out_EB_sig$se^2),  EB_df = mean((true_beta - out_EB_df_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_scam = mean((true_beta - out_EB_scam_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_po = mean((true_beta - out_EB_gampo_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_nb = mean((true_beta - out_EB_gamnb_sig$beta_EB)^2) - mean(out_EB_sig$se^2), boot = mean((true_beta - out_BR_ss_sig$beta_BR_ss)^2) - mean(out_EB_sig$se^2), FIQT = mean((true_beta - out_FIQT_sig$beta_FIQT)^2) - mean(out_EB_sig$se^2))

# Use threshold 5e-4 (need to re-apply conditional likelihood with change of threshold)
out_EB_sig <- out_EB[2*(pnorm(abs(out_EB$beta/out_EB$se), lower.tail=FALSE)) < 5e-4,]
out_EB_df_sig <- out_EB_df[2*(pnorm(abs(out_EB_df$beta/out_EB_df$se), lower.tail=FALSE)) < 5e-4,]
out_EB_scam_sig <- out_EB_scam[2*(pnorm(abs(out_EB_scam$beta/out_EB_scam$se), lower.tail=FALSE)) < 5e-4,]
out_EB_gampo_sig <- out_EB_gam_po[2*(pnorm(abs(out_EB_gam_po$beta/out_EB_gam_po$se), lower.tail=FALSE)) < 5e-4,]
out_EB_gamnb_sig <- out_EB_gam_nb[2*(pnorm(abs(out_EB_gam_nb$beta/out_EB_gam_nb$se), lower.tail=FALSE)) < 5e-4,]
out_FIQT_sig <- out_FIQT[2*(pnorm(abs(out_FIQT$beta/out_FIQT$se), lower.tail=FALSE)) < 5e-4,]
out_BR_ss_sig <- out_BR[2*(pnorm(abs(out_BR$beta/out_BR$se), lower.tail=FALSE)) < 5e-4,]

out_CL <- conditional_likelihood(summary_data = summary_stats, alpha=5e-4)

summary_data_sig2 <- summary_data_T2D_1[2*(pnorm(abs(summary_data_T2D_1$beta/summary_data_T2D_1$se), lower.tail=FALSE)) < 5e-4,]
true_beta2 <- summary_data_sig2$beta_rep

# Compute estimated MSE for each method at 5e-4 threshold
mse_T2D1_5e_4 <- data.frame(naive = mean((true_beta2 - out_EB_sig$beta)^2) - mean(out_EB_sig$se^2), CL1 = mean((true_beta2 - out_CL$beta.cl1)^2) - mean(out_EB_sig$se^2), CL2 = mean((true_beta2 - out_CL$beta.cl2)^2) - mean(out_EB_sig$se^2), CL3 = mean((true_beta2 - out_CL$beta.cl3)^2) - mean(out_EB_sig$se^2), EB = mean((true_beta2 - out_EB_sig$beta_EB)^2) - mean(out_EB_sig$se^2),  EB_df = mean((true_beta2 - out_EB_df_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_scam = mean((true_beta2 - out_EB_scam_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_po = mean((true_beta2 - out_EB_gampo_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_nb = mean((true_beta2 - out_EB_gamnb_sig$beta_EB)^2) - mean(out_EB_sig$se^2), boot = mean((true_beta2 - out_BR_ss_sig$beta_BR_ss)^2) - mean(out_EB_sig$se^2), FIQT = mean((true_beta2 - out_FIQT_sig$beta_FIQT)^2) - mean(out_EB_sig$se^2))

###############################################################################

## T2D 2

summary_data_T2D_2 <-  read.table('data/summary_data_T2D_2.txt',header=TRUE)
summary_data_T2D_2$beta <- log(summary_data_T2D_2$beta)
summary_data_T2D_2$beta_rep <- log(summary_data_T2D_2$beta_rep)

# Apply methods
summary_stats <- summary_data_T2D_2[,3:5]
out_CL <- conditional_likelihood(summary_data = summary_stats)
out_EB <- empirical_bayes(summary_data = summary_stats)
out_EB_df <- empirical_bayes(summary_data=summary_stats, method="fix_df")
out_EB_scam <- empirical_bayes(summary_data=summary_stats,method="scam")
out_EB_gam_po <- empirical_bayes(summary_data=summary_stats,method="gam_po")
out_EB_gam_nb <- empirical_bayes(summary_data=summary_stats,method="gam_nb")
out_FIQT <- FDR_IQT(summary_data = summary_stats)
out_BR <- BR_ss(summary_data = summary_stats)

# Use threshold 5e-8
out_EB_sig <- out_EB[2*(pnorm(abs(out_EB$beta/out_EB$se), lower.tail=FALSE)) < 5e-8,]
out_EB_df_sig <- out_EB_df[2*(pnorm(abs(out_EB_df$beta/out_EB_df$se), lower.tail=FALSE)) < 5e-8,]
out_EB_scam_sig <- out_EB_scam[2*(pnorm(abs(out_EB_scam$beta/out_EB_scam$se), lower.tail=FALSE)) < 5e-8,]
out_EB_gampo_sig <- out_EB_gam_po[2*(pnorm(abs(out_EB_gam_po$beta/out_EB_gam_po$se), lower.tail=FALSE)) < 5e-8,]
out_EB_gamnb_sig <- out_EB_gam_nb[2*(pnorm(abs(out_EB_gam_nb$beta/out_EB_gam_nb$se), lower.tail=FALSE)) < 5e-8,]
out_FIQT_sig <- out_FIQT[2*(pnorm(abs(out_FIQT$beta/out_FIQT$se), lower.tail=FALSE)) < 5e-8,]
out_BR_ss_sig <- out_BR[2*(pnorm(abs(out_BR$beta/out_BR$se), lower.tail=FALSE)) < 5e-8,]

summary_data_T2D_2 <- dplyr::arrange(summary_data_T2D_2,dplyr::desc(abs(summary_data_T2D_2$beta/summary_data_T2D_2$se)))
summary_data_sig_T2D_2 <- summary_data_T2D_2[2*(pnorm(abs(summary_data_T2D_2$beta/summary_data_T2D_2$se), lower.tail=FALSE)) < 5e-8,]
true_beta <- summary_data_sig_T2D_2$beta_rep

# Compute estimated MSE for each method at 5e-8 threshold
mse_T2D2 <- data.frame(naive = mean((true_beta - out_EB_sig$beta)^2) - mean(out_EB_sig$se^2), CL1 = mean((true_beta - out_CL$beta.cl1)^2) - mean(out_EB_sig$se^2), CL2 = mean((true_beta - out_CL$beta.cl2)^2) - mean(out_EB_sig$se^2), CL3 = mean((true_beta - out_CL$beta.cl3)^2) - mean(out_EB_sig$se^2), EB = mean((true_beta - out_EB_sig$beta_EB)^2) - mean(out_EB_sig$se^2),  EB_df = mean((true_beta - out_EB_df_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_scam = mean((true_beta - out_EB_scam_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_po = mean((true_beta - out_EB_gampo_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_nb = mean((true_beta - out_EB_gamnb_sig$beta_EB)^2) - mean(out_EB_sig$se^2), boot = mean((true_beta - out_BR_ss_sig$beta_BR_ss)^2) - mean(out_EB_sig$se^2), FIQT = mean((true_beta - out_FIQT_sig$beta_FIQT)^2) - mean(out_EB_sig$se^2))

# Use threshold 5e-4 (need to re-apply conditional likelihood with change of threshold)
out_EB_sig <- out_EB[2*(pnorm(abs(out_EB$beta/out_EB$se), lower.tail=FALSE)) < 5e-4,]
out_EB_df_sig <- out_EB_df[2*(pnorm(abs(out_EB_df$beta/out_EB_df$se), lower.tail=FALSE)) < 5e-4,]
out_EB_scam_sig <- out_EB_scam[2*(pnorm(abs(out_EB_scam$beta/out_EB_scam$se), lower.tail=FALSE)) < 5e-4,]
out_EB_gampo_sig <- out_EB_gam_po[2*(pnorm(abs(out_EB_gam_po$beta/out_EB_gam_po$se), lower.tail=FALSE)) < 5e-4,]
out_EB_gamnb_sig <- out_EB_gam_nb[2*(pnorm(abs(out_EB_gam_nb$beta/out_EB_gam_nb$se), lower.tail=FALSE)) < 5e-4,]
out_FIQT_sig <- out_FIQT[2*(pnorm(abs(out_FIQT$beta/out_FIQT$se), lower.tail=FALSE)) < 5e-4,]
out_BR_ss_sig <- out_BR[2*(pnorm(abs(out_BR$beta/out_BR$se), lower.tail=FALSE)) < 5e-4,]

out_CL <- conditional_likelihood(summary_data = summary_stats, alpha=5e-4)

summary_data_sig2 <- summary_data_T2D_2[2*(pnorm(abs(summary_data_T2D_2$beta/summary_data_T2D_2$se), lower.tail=FALSE)) < 5e-4,]
true_beta2 <- summary_data_sig2$beta_rep

# Compute estimated MSE for each method at 5e-4 threshold
mse_T2D2_5e_4 <- data.frame(naive = mean((true_beta2 - out_EB_sig$beta)^2) - mean(out_EB_sig$se^2), CL1 = mean((true_beta2 - out_CL$beta.cl1)^2) - mean(out_EB_sig$se^2), CL2 = mean((true_beta2 - out_CL$beta.cl2)^2) - mean(out_EB_sig$se^2), CL3 = mean((true_beta2 - out_CL$beta.cl3)^2) - mean(out_EB_sig$se^2), EB = mean((true_beta2 - out_EB_sig$beta_EB)^2) - mean(out_EB_sig$se^2),  EB_df = mean((true_beta2 - out_EB_df_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_scam = mean((true_beta2 - out_EB_scam_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_po = mean((true_beta2 - out_EB_gampo_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_nb = mean((true_beta2 - out_EB_gamnb_sig$beta_EB)^2) - mean(out_EB_sig$se^2), boot = mean((true_beta2 - out_BR_ss_sig$beta_BR_ss)^2) - mean(out_EB_sig$se^2), FIQT = mean((true_beta2 - out_FIQT_sig$beta_FIQT)^2) - mean(out_EB_sig$se^2))

###############################################################################

## Height 1

summary_data_height_1 <-  read.table('data/summary_data_height_1.txt',header=TRUE)

# Apply methods
summary_stats <- summary_data_height_1[,3:5]
out_CL <- conditional_likelihood(summary_data = summary_stats)
out_EB <- empirical_bayes(summary_data = summary_stats)
out_EB_df <- empirical_bayes(summary_data=summary_stats, method="fix_df")
out_EB_scam <- empirical_bayes(summary_data=summary_stats,method="scam")
out_EB_gam_po <- empirical_bayes(summary_data=summary_stats,method="gam_po")
out_EB_gam_nb <- empirical_bayes(summary_data=summary_stats,method="gam_nb")
out_FIQT <- FDR_IQT(summary_data = summary_stats)
out_BR <- BR_ss(summary_data = summary_stats)

# Use threshold 5e-8
out_EB_sig <- out_EB[2*(pnorm(abs(out_EB$beta/out_EB$se), lower.tail=FALSE)) < 5e-8,]
out_EB_df_sig <- out_EB_df[2*(pnorm(abs(out_EB_df$beta/out_EB_df$se), lower.tail=FALSE)) < 5e-8,]
out_EB_scam_sig <- out_EB_scam[2*(pnorm(abs(out_EB_scam$beta/out_EB_scam$se), lower.tail=FALSE)) < 5e-8,]
out_EB_gampo_sig <- out_EB_gam_po[2*(pnorm(abs(out_EB_gam_po$beta/out_EB_gam_po$se), lower.tail=FALSE)) < 5e-8,]
out_EB_gamnb_sig <- out_EB_gam_nb[2*(pnorm(abs(out_EB_gam_nb$beta/out_EB_gam_nb$se), lower.tail=FALSE)) < 5e-8,]
out_FIQT_sig <- out_FIQT[2*(pnorm(abs(out_FIQT$beta/out_FIQT$se), lower.tail=FALSE)) < 5e-8,]
out_BR_ss_sig <- out_BR[2*(pnorm(abs(out_BR$beta/out_BR$se), lower.tail=FALSE)) < 5e-8,]

summary_data_height_1 <- dplyr::arrange(summary_data_height_1,dplyr::desc(abs(summary_data_height_1$beta/summary_data_height_1$se)))
summary_data_sig_height_1 <- summary_data_height_1[2*(pnorm(abs(summary_data_height_1$beta/summary_data_height_1$se), lower.tail=FALSE)) < 5e-8,]
true_beta <- summary_data_sig_height_1$beta_rep

# Compute estimated MSE for each method at 5e-8 threshold
mse_height1 <- data.frame(naive = mean((true_beta - out_EB_sig$beta)^2) - mean(out_EB_sig$se^2), CL1 = mean((true_beta - out_CL$beta.cl1)^2) - mean(out_EB_sig$se^2), CL2 = mean((true_beta - out_CL$beta.cl2)^2) - mean(out_EB_sig$se^2), CL3 = mean((true_beta - out_CL$beta.cl3)^2) - mean(out_EB_sig$se^2), EB = mean((true_beta - out_EB_sig$beta_EB)^2) - mean(out_EB_sig$se^2),  EB_df = mean((true_beta - out_EB_df_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_scam = mean((true_beta - out_EB_scam_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_po = mean((true_beta - out_EB_gampo_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_nb = mean((true_beta - out_EB_gamnb_sig$beta_EB)^2) - mean(out_EB_sig$se^2), boot = mean((true_beta - out_BR_ss_sig$beta_BR_ss)^2) - mean(out_EB_sig$se^2), FIQT = mean((true_beta - out_FIQT_sig$beta_FIQT)^2) - mean(out_EB_sig$se^2))

# Use threshold 5e-4 (need to re-apply conditional likelihood with change of threshold)
out_EB_sig <- out_EB[2*(pnorm(abs(out_EB$beta/out_EB$se), lower.tail=FALSE)) < 5e-4,]
out_EB_df_sig <- out_EB_df[2*(pnorm(abs(out_EB_df$beta/out_EB_df$se), lower.tail=FALSE)) < 5e-4,]
out_EB_scam_sig <- out_EB_scam[2*(pnorm(abs(out_EB_scam$beta/out_EB_scam$se), lower.tail=FALSE)) < 5e-4,]
out_EB_gampo_sig <- out_EB_gam_po[2*(pnorm(abs(out_EB_gam_po$beta/out_EB_gam_po$se), lower.tail=FALSE)) < 5e-4,]
out_EB_gamnb_sig <- out_EB_gam_nb[2*(pnorm(abs(out_EB_gam_nb$beta/out_EB_gam_nb$se), lower.tail=FALSE)) < 5e-4,]
out_FIQT_sig <- out_FIQT[2*(pnorm(abs(out_FIQT$beta/out_FIQT$se), lower.tail=FALSE)) < 5e-4,]
out_BR_ss_sig <- out_BR[2*(pnorm(abs(out_BR$beta/out_BR$se), lower.tail=FALSE)) < 5e-4,]

out_CL <- conditional_likelihood(summary_data = summary_stats, alpha=5e-4)

summary_data_sig2 <- summary_data_height_1[2*(pnorm(abs(summary_data_height_1$beta/summary_data_height_1$se), lower.tail=FALSE)) < 5e-4,]
true_beta2 <- summary_data_sig2$beta_rep

# Compute estimated MSE for each method at 5e-4 threshold
mse_height1_5e_4 <- data.frame(naive = mean((true_beta2 - out_EB_sig$beta)^2) - mean(out_EB_sig$se^2), CL1 = mean((true_beta2 - out_CL$beta.cl1)^2) - mean(out_EB_sig$se^2), CL2 = mean((true_beta2 - out_CL$beta.cl2)^2) - mean(out_EB_sig$se^2), CL3 = mean((true_beta2 - out_CL$beta.cl3)^2) - mean(out_EB_sig$se^2), EB = mean((true_beta2 - out_EB_sig$beta_EB)^2) - mean(out_EB_sig$se^2),  EB_df = mean((true_beta2 - out_EB_df_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_scam = mean((true_beta2 - out_EB_scam_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_po = mean((true_beta2 - out_EB_gampo_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_nb = mean((true_beta2 - out_EB_gamnb_sig$beta_EB)^2) - mean(out_EB_sig$se^2), boot = mean((true_beta2 - out_BR_ss_sig$beta_BR_ss)^2) - mean(out_EB_sig$se^2), FIQT = mean((true_beta2 - out_FIQT_sig$beta_FIQT)^2) - mean(out_EB_sig$se^2))

###############################################################################

## Height 2

summary_data_height_2 <-  read.table('data/summary_data_height_2.txt',header=TRUE)

# Apply methods
summary_stats <- summary_data_height_2[,3:5]
out_CL <- conditional_likelihood(summary_data = summary_stats)
out_EB <- empirical_bayes(summary_data = summary_stats)
out_EB_df <- empirical_bayes(summary_data=summary_stats, method="fix_df")
out_EB_scam <- empirical_bayes(summary_data=summary_stats,method="scam")
out_EB_gam_po <- empirical_bayes(summary_data=summary_stats,method="gam_po")
out_EB_gam_nb <- empirical_bayes(summary_data=summary_stats,method="gam_nb")
out_FIQT <- FDR_IQT(summary_data = summary_stats)
out_BR <- BR_ss(summary_data = summary_stats)

# Use threshold 5e-8
out_EB_sig <- out_EB[2*(pnorm(abs(out_EB$beta/out_EB$se), lower.tail=FALSE)) < 5e-8,]
out_EB_df_sig <- out_EB_df[2*(pnorm(abs(out_EB_df$beta/out_EB_df$se), lower.tail=FALSE)) < 5e-8,]
out_EB_scam_sig <- out_EB_scam[2*(pnorm(abs(out_EB_scam$beta/out_EB_scam$se), lower.tail=FALSE)) < 5e-8,]
out_EB_gampo_sig <- out_EB_gam_po[2*(pnorm(abs(out_EB_gam_po$beta/out_EB_gam_po$se), lower.tail=FALSE)) < 5e-8,]
out_EB_gamnb_sig <- out_EB_gam_nb[2*(pnorm(abs(out_EB_gam_nb$beta/out_EB_gam_nb$se), lower.tail=FALSE)) < 5e-8,]
out_FIQT_sig <- out_FIQT[2*(pnorm(abs(out_FIQT$beta/out_FIQT$se), lower.tail=FALSE)) < 5e-8,]
out_BR_ss_sig <- out_BR[2*(pnorm(abs(out_BR$beta/out_BR$se), lower.tail=FALSE)) < 5e-8,]

summary_data_height_2 <- dplyr::arrange(summary_data_height_2,dplyr::desc(abs(summary_data_height_2$beta/summary_data_height_2$se)))
summary_data_sig_height_2 <- summary_data_height_2[2*(pnorm(abs(summary_data_height_2$beta/summary_data_height_2$se), lower.tail=FALSE)) < 5e-8,]
true_beta <- summary_data_sig_height_2$beta_rep

# Compute estimated MSE for each method at 5e-8 threshold
mse_height2 <- data.frame(naive = mean((true_beta - out_EB_sig$beta)^2) - mean(out_EB_sig$se^2), CL1 = mean((true_beta - out_CL$beta.cl1)^2) - mean(out_EB_sig$se^2), CL2 = mean((true_beta - out_CL$beta.cl2)^2) - mean(out_EB_sig$se^2), CL3 = mean((true_beta - out_CL$beta.cl3)^2) - mean(out_EB_sig$se^2), EB = mean((true_beta - out_EB_sig$beta_EB)^2) - mean(out_EB_sig$se^2),  EB_df = mean((true_beta - out_EB_df_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_scam = mean((true_beta - out_EB_scam_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_po = mean((true_beta - out_EB_gampo_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_nb = mean((true_beta - out_EB_gamnb_sig$beta_EB)^2) - mean(out_EB_sig$se^2), boot = mean((true_beta - out_BR_ss_sig$beta_BR_ss)^2) - mean(out_EB_sig$se^2), FIQT = mean((true_beta - out_FIQT_sig$beta_FIQT)^2) - mean(out_EB_sig$se^2))

# Use threshold 5e-4 (need to re-apply conditional likelihood with change of threshold)
out_EB_sig <- out_EB[2*(pnorm(abs(out_EB$beta/out_EB$se), lower.tail=FALSE)) < 5e-4,]
out_EB_df_sig <- out_EB_df[2*(pnorm(abs(out_EB_df$beta/out_EB_df$se), lower.tail=FALSE)) < 5e-4,]
out_EB_scam_sig <- out_EB_scam[2*(pnorm(abs(out_EB_scam$beta/out_EB_scam$se), lower.tail=FALSE)) < 5e-4,]
out_EB_gampo_sig <- out_EB_gam_po[2*(pnorm(abs(out_EB_gam_po$beta/out_EB_gam_po$se), lower.tail=FALSE)) < 5e-4,]
out_EB_gamnb_sig <- out_EB_gam_nb[2*(pnorm(abs(out_EB_gam_nb$beta/out_EB_gam_nb$se), lower.tail=FALSE)) < 5e-4,]
out_FIQT_sig <- out_FIQT[2*(pnorm(abs(out_FIQT$beta/out_FIQT$se), lower.tail=FALSE)) < 5e-4,]
out_BR_ss_sig <- out_BR[2*(pnorm(abs(out_BR$beta/out_BR$se), lower.tail=FALSE)) < 5e-4,]

out_CL <- conditional_likelihood(summary_data = summary_stats, alpha=5e-4)

summary_data_sig2 <- summary_data_height_2[2*(pnorm(abs(summary_data_height_2$beta/summary_data_height_2$se), lower.tail=FALSE)) < 5e-4,]
true_beta2 <- summary_data_sig2$beta_rep

# Compute estimated MSE for each method at 5e-4 threshold
mse_height2_5e_4 <- data.frame(naive = mean((true_beta2 - out_EB_sig$beta)^2) - mean(out_EB_sig$se^2), CL1 = mean((true_beta2 - out_CL$beta.cl1)^2) - mean(out_EB_sig$se^2), CL2 = mean((true_beta2 - out_CL$beta.cl2)^2) - mean(out_EB_sig$se^2), CL3 = mean((true_beta2 - out_CL$beta.cl3)^2) - mean(out_EB_sig$se^2), EB = mean((true_beta2 - out_EB_sig$beta_EB)^2) - mean(out_EB_sig$se^2),  EB_df = mean((true_beta2 - out_EB_df_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_scam = mean((true_beta2 - out_EB_scam_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_po = mean((true_beta2 - out_EB_gampo_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_nb = mean((true_beta2 - out_EB_gamnb_sig$beta_EB)^2) - mean(out_EB_sig$se^2), boot = mean((true_beta2 - out_BR_ss_sig$beta_BR_ss)^2) - mean(out_EB_sig$se^2), FIQT = mean((true_beta2 - out_FIQT_sig$beta_FIQT)^2) - mean(out_EB_sig$se^2))

###############################################################################

## Combine all above results and save them in external results files 

mse_5e_8 <- rbind(mse_BMI1,mse_BMI2, mse_T2D1, mse_T2D2, mse_height1, mse_height2)
mse_5e_8 <- round(mse_5e_8[,1:11],5)
GWAS <- c("BMI 1", "BMI 2", "T2D 1", "T2D 2", "Height 1", "Height 2")
mse_5e_8 <- cbind(GWAS,mse_5e_8)
write.csv(mse_5e_8,"results/mse_5e_8.txt", row.names = FALSE)

mse_5e_4 <- rbind(mse_BMI1_5e_4,mse_BMI2_5e_4, mse_T2D1_5e_4, mse_T2D2_5e_4, mse_height1_5e_4, mse_height2_5e_4)
mse_5e_4 <- round(mse_5e_4[,1:11],5)
GWAS <- c("BMI 1", "BMI 2", "T2D 1", "T2D 2", "Height 1", "Height 2")
mse_5e_4 <- cbind(GWAS,mse_5e_4)
write.csv(mse_5e_4,"results/mse_5e_4.txt", row.names = FALSE)

###############################################################################
###############################################################################

## PART 2) PLOT RESULTS APPROPRIATELY


new_mse_5e_8 <- pivot_longer(mse_5e_8, -c(GWAS), values_to = "MSE", names_to = "Method")
new_mse_5e_4 <- pivot_longer(mse_5e_4, -c(GWAS), values_to = "MSE", names_to = "Method")

new_mse_5e_8$Method <- factor(new_mse_5e_8$Method, levels=c("naive", "CL1", "CL2", "CL3", "EB", "EB_df", "EB_gam_nb", "EB_gam_po", "EB_scam", "boot", "FIQT"))
new_mse_5e_4$Method <- factor(new_mse_5e_4$Method, levels=c("naive", "CL1", "CL2", "CL3", "EB", "EB_df", "EB_gam_nb", "EB_gam_po", "EB_scam", "boot", "FIQT"))
new_mse_5e_8$GWAS <- factor(new_mse_5e_8$GWAS, levels=c("BMI 1", "BMI 2", "T2D 1", "T2D 2", "Height 1", "Height 2"))
new_mse_5e_4$GWAS <- factor(new_mse_5e_4$GWAS, levels=c("BMI 1", "BMI 2", "T2D 1", "T2D 2", "Height 1", "Height 2"))


plot_1 <- ggplot(new_mse_5e_8,aes(x=MSE,y=Method,fill=Method)) + geom_bar(stat="identity")+ geom_col(width=1) + facet_wrap(.~GWAS,nrow=6, strip.position = "right") + scale_shape_manual(values=c(6,20,20,20,15,18,18,18,18,10,15)) + scale_fill_manual(values=c(col[1],col[2],col[8],col[4],col[5],col[3],col[6],col[7],col1[1],col1[11],col1[9])) + ylab("Method") + 
  xlab(expression(paste(italic("Estimated"), " MSE at ", 5%*%10^-8))) + theme_bw() + theme(text = element_text(size=12),legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text.y = element_text(angle = 0)) + geom_vline(data = new_mse_5e_8[new_mse_5e_8$Method=="naive",], aes(xintercept = MSE), linetype=2, size=0.5) + geom_vline(xintercept=0)

plot_2 <- ggplot(new_mse_5e_4,aes(x=MSE,y=Method,fill=Method)) + geom_bar(stat="identity")+ geom_col(width=1) + facet_wrap(GWAS~.,nrow=6, strip.position = "right") + scale_shape_manual(values=c(6,20,20,20,15,18,18,18,18,10,15)) + scale_fill_manual(values=c(col[1],col[2],col[8],col[4],col[5],col[3],col[6],col[7],col1[1],col1[11],col1[9])) + ylab("Method") + 
  xlab(expression(paste(italic("Estimated"), " MSE at ", 5%*%10^-4))) + theme_bw() + theme(text = element_text(size=12),legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text.y = element_text(angle=0)) + geom_vline(data = new_mse_5e_4[new_mse_5e_4$Method=="naive",], aes(xintercept = MSE), linetype=2, size=0.5) + geom_vline(xintercept=0)


plot_1
plot_2

###############################################################################
###############################################################################

## PART 3) EXPLORATORY Z vs BIAS PLOTS 

## BMI GWAS 1
out <- summary_data_bmi_1
title <- "BMI GWAS 1"
title_col <- col[1]
out$bias <- out$beta - out$beta_rep
out$z <- out$beta/out$se
subout <- out[(abs(out$beta) > (abs(out$beta_rep) + 1.96*out$se)) & (abs(out$beta/out$se) > qnorm(1-(5e-4)/2)),]
BMI1 <- ggplot(out,aes(x=z,y=bias)) + geom_point(size=1) + geom_point(data=subout, 
                                                                aes(x=z,y=bias), 
                                                                color='grey40',
                                                                size=1) + xlab("z") +
    ylab("bias") + ggtitle(title) +  theme_classic() + geom_vline(xintercept=qnorm(1-(5e-8)/2), colour="red", linetype="dashed",size=1) + geom_vline(xintercept=-qnorm(1-(5e-8)/2), colour="red", linetype="dashed",size=1) +
    geom_vline(xintercept=qnorm(1-(5e-4)/2), colour="darkred", linetype="dashed",size=1) + geom_vline(xintercept=-qnorm(1-(5e-4)/2), colour="darkred", linetype="dashed",size=1) + geom_hline(yintercept=0) +
    theme(text = element_text(size=12),
      plot.title = element_text(hjust = 0.5, face="bold", colour=title_col),
    )


## BMI GWAS 2
out <- summary_data_bmi_2
title <- "BMI GWAS 2"
title_col <- col[2]
out$bias <- out$beta - out$beta_rep
out$z <- out$beta/out$se
subout <- out[(abs(out$beta) > (abs(out$beta_rep) + 1.96*out$se)) & (abs(out$beta/out$se) > qnorm(1-(5e-4)/2)),]
BMI2 <- ggplot(out,aes(x=z,y=bias)) + geom_point(size=1) + geom_point(data=subout, 
                                                                      aes(x=z,y=bias), 
                                                                      color='grey40',
                                                                      size=1) + xlab("z") +
  ylab("bias") + ggtitle(title) +  theme_classic() + geom_vline(xintercept=qnorm(1-(5e-8)/2), colour="red", linetype="dashed",size=1) + geom_vline(xintercept=-qnorm(1-(5e-8)/2), colour="red", linetype="dashed",size=1) +
  geom_vline(xintercept=qnorm(1-(5e-4)/2), colour="darkred", linetype="dashed",size=1) + geom_vline(xintercept=-qnorm(1-(5e-4)/2), colour="darkred", linetype="dashed",size=1) + geom_hline(yintercept=0) +
  theme(text = element_text(size=12),
        plot.title = element_text(hjust = 0.5, face="bold", colour=title_col),
  )

## T2D GWAS 1
out <- summary_data_T2D_1
title <- "T2D GWAS 1"
title_col <- col[3]
out$bias <- out$beta - out$beta_rep
out$z <- out$beta/out$se
subout <- out[(abs(out$beta) > (abs(out$beta_rep) + 1.96*out$se)) & (abs(out$beta/out$se) > qnorm(1-(5e-4)/2)),]
T2D1 <- ggplot(out,aes(x=z,y=bias)) + geom_point(size=1) + geom_point(data=subout, 
                                                                      aes(x=z,y=bias), 
                                                                      color='grey40',
                                                                      size=1) + xlab("z") +
  ylab("bias") + ggtitle(title) +  theme_classic() + geom_vline(xintercept=qnorm(1-(5e-8)/2), colour="red", linetype="dashed",size=1) + geom_vline(xintercept=-qnorm(1-(5e-8)/2), colour="red", linetype="dashed",size=1) +
  geom_vline(xintercept=qnorm(1-(5e-4)/2), colour="darkred", linetype="dashed",size=1) + geom_vline(xintercept=-qnorm(1-(5e-4)/2), colour="darkred", linetype="dashed",size=1) + geom_hline(yintercept=0) +
  theme(text = element_text(size=12),
        plot.title = element_text(hjust = 0.5, face="bold", colour=title_col),
  )

## T2D GWAS 2
out <- summary_data_T2D_2
title <- "T2D GWAS 2"
title_col <- col[4]
out$bias <- out$beta - out$beta_rep
out$z <- out$beta/out$se
subout <- out[(abs(out$beta) > (abs(out$beta_rep) + 1.96*out$se)) & (abs(out$beta/out$se) > qnorm(1-(5e-4)/2)),]
T2D2 <- ggplot(out,aes(x=z,y=bias)) + geom_point(size=1) + geom_point(data=subout, 
                                                                      aes(x=z,y=bias), 
                                                                      color='grey40',
                                                                      size=1) + xlab("z") +
  ylab("bias") + ggtitle(title) +  theme_classic() + geom_vline(xintercept=qnorm(1-(5e-8)/2), colour="red", linetype="dashed",size=1) + geom_vline(xintercept=-qnorm(1-(5e-8)/2), colour="red", linetype="dashed",size=1) +
  geom_vline(xintercept=qnorm(1-(5e-4)/2), colour="darkred", linetype="dashed",size=1) + geom_vline(xintercept=-qnorm(1-(5e-4)/2), colour="darkred", linetype="dashed",size=1) + geom_hline(yintercept=0) +
  theme(text = element_text(size=12),
        plot.title = element_text(hjust = 0.5, face="bold", colour=title_col),
  )


## Height GWAS 1
out <- summary_data_height_1
title <- "Height GWAS 1"
title_col <- col[5]
out$bias <- out$beta - out$beta_rep
out$z <- out$beta/out$se
subout <- out[(abs(out$beta) > (abs(out$beta_rep) + 1.96*out$se)) & (abs(out$beta/out$se) > qnorm(1-(5e-4)/2)),]
Height1 <- ggplot(out,aes(x=z,y=bias)) + geom_point(size=1) + geom_point(data=subout, 
                                                                      aes(x=z,y=bias), 
                                                                      color='grey40',
                                                                      size=1) + xlab("z") +
  ylab("bias") + ggtitle(title) +  theme_classic() + geom_vline(xintercept=qnorm(1-(5e-8)/2), colour="red", linetype="dashed",size=1) + geom_vline(xintercept=-qnorm(1-(5e-8)/2), colour="red", linetype="dashed",size=1) +
  geom_vline(xintercept=qnorm(1-(5e-4)/2), colour="darkred", linetype="dashed",size=1) + geom_vline(xintercept=-qnorm(1-(5e-4)/2), colour="darkred", linetype="dashed",size=1) + geom_hline(yintercept=0) +
  theme(text = element_text(size=12),
        plot.title = element_text(hjust = 0.5, face="bold", colour=title_col),
  )


## Height GWAS 2
out <- summary_data_height_2
title <- "Height GWAS 2"
title_col <- col[6]
out$bias <- out$beta - out$beta_rep
out$z <- out$beta/out$se
subout <- out[(abs(out$beta) > (abs(out$beta_rep) + 1.96*out$se)) & (abs(out$beta/out$se) > qnorm(1-(5e-4)/2)),]
Height2 <- ggplot(out,aes(x=z,y=bias)) + geom_point(size=1) + geom_point(data=subout, 
                                                                      aes(x=z,y=bias), 
                                                                      color='grey40',
                                                                      size=1) + xlab("z") +
  ylab("bias") + ggtitle(title) +  theme_classic() + geom_vline(xintercept=qnorm(1-(5e-8)/2), colour="red", linetype="dashed",size=1) + geom_vline(xintercept=-qnorm(1-(5e-8)/2), colour="red", linetype="dashed",size=1) +
  geom_vline(xintercept=qnorm(1-(5e-4)/2), colour="darkred", linetype="dashed",size=1) + geom_vline(xintercept=-qnorm(1-(5e-4)/2), colour="darkred", linetype="dashed",size=1) + geom_hline(yintercept=0) +
  theme(text = element_text(size=12),
        plot.title = element_text(hjust = 0.5, face="bold", colour=title_col),
  )


## Combine above six plots into a figure together
ggarrange(BMI1, BMI2, T2D1, T2D2, Height1, Height2, labels = c("A", "B", "C", "D", "E", "F"), ncol = 2, nrow = 3)

################################################################################
################################################################################

## PART 4) COMPUTING AVERAGE BIAS FOR EACH DATASET, METHOD AND THRESHOLD


## Average bias is simply defined as the mean of adjusted discovery association
## estimate minus the replication association estimate among significant SNPs

## In order for bias not to cancel each other out - we need to compute
## separately for both up and down significant SNPs, i.e. significant SNPs with
## negative and positive association estimates

## BMI 1

summary_data_bmi_1 <-  read.table('data/summary_data_bmi_1.txt',header=TRUE)

# Apply methods
summary_stats <- summary_data_bmi_1[,3:5]
out_CL <- conditional_likelihood(summary_data = summary_stats)
out_EB <- empirical_bayes(summary_data = summary_stats)
out_EB_df <- empirical_bayes(summary_data=summary_stats, method="fix_df")
out_EB_scam <- empirical_bayes(summary_data=summary_stats,method="scam")
out_EB_gam_po <- empirical_bayes(summary_data=summary_stats,method="gam_po")
out_EB_gam_nb <- empirical_bayes(summary_data=summary_stats,method="gam_nb")
out_FIQT <- FDR_IQT(summary_data = summary_stats)
out_BR <- BR_ss(summary_data = summary_stats)

# Use threshold 5e-8 (positive)
out_CL_sig <- out_CL[out_CL$beta/out_CL$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_EB_sig <- out_EB[out_EB$beta/out_EB$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_EB_df_sig <- out_EB_df[out_EB_df$beta/out_EB_df$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_EB_scam_sig <- out_EB_scam[out_EB_scam$beta/out_EB_scam$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_EB_gampo_sig <- out_EB_gam_po[out_EB_gam_po$beta/out_EB_gam_po$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_EB_gamnb_sig <- out_EB_gam_nb[out_EB_gam_nb$beta/out_EB_gam_nb$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_FIQT_sig <- out_FIQT[out_FIQT$beta/out_FIQT$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_BR_ss_sig <- out_BR[out_BR$beta/out_BR$se > qnorm(5e-8/2, lower.tail=FALSE),]

summary_data_bmi_1 <- dplyr::arrange(summary_data_bmi_1,dplyr::desc(abs(summary_data_bmi_1$beta/summary_data_bmi_1$se)))
summary_data_sig_bmi_1 <- summary_data_bmi_1[summary_data_bmi_1$beta/summary_data_bmi_1$se > qnorm(5e-8/2, lower.tail=FALSE),]
true_beta <- summary_data_sig_bmi_1$beta_rep

# Compute average bias for each method for 5e-8 threshold (positive)
bias_BMI1_up <- data.frame(naive = mean((out_EB_sig$beta - true_beta)), CL1 = mean((out_CL_sig$beta.cl1 - true_beta)), CL2 = mean((out_CL_sig$beta.cl2 - true_beta)), CL3 = mean((out_CL_sig$beta.cl3 - true_beta)), EB = mean((out_EB_sig$beta_EB - true_beta)),  EB_df = mean((out_EB_df_sig$beta_EB - true_beta)), EB_scam = mean((out_EB_scam_sig$beta_EB - true_beta)), EB_gam_po = mean((out_EB_gampo_sig$beta_EB - true_beta)), EB_gam_nb = mean((out_EB_gamnb_sig$beta_EB - true_beta)), boot = mean((out_BR_ss_sig$beta_BR_ss - true_beta)), FIQT = mean((out_FIQT_sig$beta_FIQT - true_beta)))


# Use threshold 5e-8 (negative)
out_CL_sig <- out_CL[out_CL$beta/out_CL$se < qnorm(5e-8/2,),]
out_EB_sig <- out_EB[out_EB$beta/out_EB$se < qnorm(5e-8/2),]
out_EB_df_sig <- out_EB_df[out_EB_df$beta/out_EB_df$se < qnorm(5e-8/2),]
out_EB_scam_sig <- out_EB_scam[out_EB_scam$beta/out_EB_scam$se < qnorm(5e-8/2),]
out_EB_gampo_sig <- out_EB_gam_po[out_EB_gam_po$beta/out_EB_gam_po$se < qnorm(5e-8/2),]
out_EB_gamnb_sig <- out_EB_gam_nb[out_EB_gam_nb$beta/out_EB_gam_nb$se < qnorm(5e-8/2),]
out_FIQT_sig <- out_FIQT[out_FIQT$beta/out_FIQT$se < qnorm(5e-8/2),]
out_BR_ss_sig <- out_BR[out_BR$beta/out_BR$se < qnorm(5e-8/2),]

summary_data_bmi_1 <- dplyr::arrange(summary_data_bmi_1,dplyr::desc(abs(summary_data_bmi_1$beta/summary_data_bmi_1$se)))
summary_data_sig_bmi_1 <- summary_data_bmi_1[summary_data_bmi_1$beta/summary_data_bmi_1$se < qnorm(5e-8/2),]
true_beta <- summary_data_sig_bmi_1$beta_rep

# Compute average bias for each method for 5e-8 threshold (negative)
bias_BMI1_down <- data.frame(naive = mean((out_EB_sig$beta - true_beta)), CL1 = mean((out_CL_sig$beta.cl1 - true_beta)), CL2 = mean((out_CL_sig$beta.cl2 - true_beta)), CL3 = mean((out_CL_sig$beta.cl3 - true_beta)), EB = mean((out_EB_sig$beta_EB - true_beta)),  EB_df = mean((out_EB_df_sig$beta_EB - true_beta)), EB_scam = mean((out_EB_scam_sig$beta_EB - true_beta)), EB_gam_po = mean((out_EB_gampo_sig$beta_EB - true_beta)), EB_gam_nb = mean((out_EB_gamnb_sig$beta_EB - true_beta)), boot = mean((out_BR_ss_sig$beta_BR_ss - true_beta)), FIQT = mean((out_FIQT_sig$beta_FIQT - true_beta)))


# Use threshold 5e-4 (need to re-apply conditional likelihood with change of threshold)

out_CL <- conditional_likelihood(summary_data = summary_stats, alpha=5e-4)

# Use threshold 5e-4 (positive)
out_CL_sig <- out_CL[out_CL$beta/out_CL$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_EB_sig <- out_EB[out_EB$beta/out_EB$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_EB_df_sig <- out_EB_df[out_EB_df$beta/out_EB_df$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_EB_scam_sig <- out_EB_scam[out_EB_scam$beta/out_EB_scam$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_EB_gampo_sig <- out_EB_gam_po[out_EB_gam_po$beta/out_EB_gam_po$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_EB_gamnb_sig <- out_EB_gam_nb[out_EB_gam_nb$beta/out_EB_gam_nb$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_FIQT_sig <- out_FIQT[out_FIQT$beta/out_FIQT$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_BR_ss_sig <- out_BR[out_BR$beta/out_BR$se > qnorm(5e-4/2, lower.tail=FALSE),]

summary_data_bmi_1 <- dplyr::arrange(summary_data_bmi_1,dplyr::desc(abs(summary_data_bmi_1$beta/summary_data_bmi_1$se)))
summary_data_sig_bmi_1 <- summary_data_bmi_1[summary_data_bmi_1$beta/summary_data_bmi_1$se > qnorm(5e-4/2, lower.tail=FALSE),]
true_beta <- summary_data_sig_bmi_1$beta_rep

# Compute average bias for each method for 5e-4 threshold (positive)
bias_BMI1_5e_4_up <- data.frame(naive = mean((out_EB_sig$beta - true_beta)), CL1 = mean((out_CL_sig$beta.cl1 - true_beta)), CL2 = mean((out_CL_sig$beta.cl2 - true_beta)), CL3 = mean((out_CL_sig$beta.cl3 - true_beta)), EB = mean((out_EB_sig$beta_EB - true_beta)),  EB_df = mean((out_EB_df_sig$beta_EB - true_beta)), EB_scam = mean((out_EB_scam_sig$beta_EB - true_beta)), EB_gam_po = mean((out_EB_gampo_sig$beta_EB - true_beta)), EB_gam_nb = mean((out_EB_gamnb_sig$beta_EB - true_beta)), boot = mean((out_BR_ss_sig$beta_BR_ss - true_beta)), FIQT = mean((out_FIQT_sig$beta_FIQT - true_beta)))


# Use threshold 5e-4 (negative)
out_CL_sig <- out_CL[out_CL$beta/out_CL$se < qnorm(5e-4/2,),]
out_EB_sig <- out_EB[out_EB$beta/out_EB$se < qnorm(5e-4/2),]
out_EB_df_sig <- out_EB_df[out_EB_df$beta/out_EB_df$se < qnorm(5e-4/2),]
out_EB_scam_sig <- out_EB_scam[out_EB_scam$beta/out_EB_scam$se < qnorm(5e-4/2),]
out_EB_gampo_sig <- out_EB_gam_po[out_EB_gam_po$beta/out_EB_gam_po$se < qnorm(5e-4/2),]
out_EB_gamnb_sig <- out_EB_gam_nb[out_EB_gam_nb$beta/out_EB_gam_nb$se < qnorm(5e-4/2),]
out_FIQT_sig <- out_FIQT[out_FIQT$beta/out_FIQT$se < qnorm(5e-4/2),]
out_BR_ss_sig <- out_BR[out_BR$beta/out_BR$se < qnorm(5e-4/2),]

summary_data_bmi_1 <- dplyr::arrange(summary_data_bmi_1,dplyr::desc(abs(summary_data_bmi_1$beta/summary_data_bmi_1$se)))
summary_data_sig_bmi_1 <- summary_data_bmi_1[summary_data_bmi_1$beta/summary_data_bmi_1$se < qnorm(5e-4/2),]
true_beta <- summary_data_sig_bmi_1$beta_rep

# Compute average bias for each method for 5e-4 threshold (negative)
bias_BMI1_5e_4_down <- data.frame(naive = mean((out_EB_sig$beta - true_beta)), CL1 = mean((out_CL_sig$beta.cl1 - true_beta)), CL2 = mean((out_CL_sig$beta.cl2 - true_beta)), CL3 = mean((out_CL_sig$beta.cl3 - true_beta)), EB = mean((out_EB_sig$beta_EB - true_beta)),  EB_df = mean((out_EB_df_sig$beta_EB - true_beta)), EB_scam = mean((out_EB_scam_sig$beta_EB - true_beta)), EB_gam_po = mean((out_EB_gampo_sig$beta_EB - true_beta)), EB_gam_nb = mean((out_EB_gamnb_sig$beta_EB - true_beta)), boot = mean((out_BR_ss_sig$beta_BR_ss - true_beta)), FIQT = mean((out_FIQT_sig$beta_FIQT - true_beta)))

###############################################################################

## BMI 2

summary_data_bmi_2 <-  read.table('data/summary_data_bmi_2.txt',header=TRUE)

# Apply methods
summary_stats <- summary_data_bmi_2[,3:5]
out_CL <- conditional_likelihood(summary_data = summary_stats)
out_EB <- empirical_bayes(summary_data = summary_stats)
out_EB_df <- empirical_bayes(summary_data=summary_stats, method="fix_df")
out_EB_scam <- empirical_bayes(summary_data=summary_stats,method="scam")
out_EB_gam_po <- empirical_bayes(summary_data=summary_stats,method="gam_po")
out_EB_gam_nb <- empirical_bayes(summary_data=summary_stats,method="gam_nb")
out_FIQT <- FDR_IQT(summary_data = summary_stats)
out_BR <- BR_ss(summary_data = summary_stats)

# Use threshold 5e-8 (positive)
out_CL_sig <- out_CL[out_CL$beta/out_CL$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_EB_sig <- out_EB[out_EB$beta/out_EB$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_EB_df_sig <- out_EB_df[out_EB_df$beta/out_EB_df$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_EB_scam_sig <- out_EB_scam[out_EB_scam$beta/out_EB_scam$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_EB_gampo_sig <- out_EB_gam_po[out_EB_gam_po$beta/out_EB_gam_po$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_EB_gamnb_sig <- out_EB_gam_nb[out_EB_gam_nb$beta/out_EB_gam_nb$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_FIQT_sig <- out_FIQT[out_FIQT$beta/out_FIQT$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_BR_ss_sig <- out_BR[out_BR$beta/out_BR$se > qnorm(5e-8/2, lower.tail=FALSE),]

summary_data_bmi_2 <- dplyr::arrange(summary_data_bmi_2,dplyr::desc(abs(summary_data_bmi_2$beta/summary_data_bmi_2$se)))
summary_data_sig_bmi_2 <- summary_data_bmi_2[summary_data_bmi_2$beta/summary_data_bmi_2$se > qnorm(5e-8/2, lower.tail=FALSE),]
true_beta <- summary_data_sig_bmi_2$beta_rep

# Compute average bias for each method for 5e-8 threshold (positive)
bias_BMI2_up <- data.frame(naive = mean((out_EB_sig$beta - true_beta)), CL1 = mean((out_CL_sig$beta.cl1 - true_beta)), CL2 = mean((out_CL_sig$beta.cl2 - true_beta)), CL3 = mean((out_CL_sig$beta.cl3 - true_beta)), EB = mean((out_EB_sig$beta_EB - true_beta)),  EB_df = mean((out_EB_df_sig$beta_EB - true_beta)), EB_scam = mean((out_EB_scam_sig$beta_EB - true_beta)), EB_gam_po = mean((out_EB_gampo_sig$beta_EB - true_beta)), EB_gam_nb = mean((out_EB_gamnb_sig$beta_EB - true_beta)), boot = mean((out_BR_ss_sig$beta_BR_ss - true_beta)), FIQT = mean((out_FIQT_sig$beta_FIQT - true_beta)))


# Use threshold 5e-8 (negative)
out_CL_sig <- out_CL[out_CL$beta/out_CL$se < qnorm(5e-8/2,),]
out_EB_sig <- out_EB[out_EB$beta/out_EB$se < qnorm(5e-8/2),]
out_EB_df_sig <- out_EB_df[out_EB_df$beta/out_EB_df$se < qnorm(5e-8/2),]
out_EB_scam_sig <- out_EB_scam[out_EB_scam$beta/out_EB_scam$se < qnorm(5e-8/2),]
out_EB_gampo_sig <- out_EB_gam_po[out_EB_gam_po$beta/out_EB_gam_po$se < qnorm(5e-8/2),]
out_EB_gamnb_sig <- out_EB_gam_nb[out_EB_gam_nb$beta/out_EB_gam_nb$se < qnorm(5e-8/2),]
out_FIQT_sig <- out_FIQT[out_FIQT$beta/out_FIQT$se < qnorm(5e-8/2),]
out_BR_ss_sig <- out_BR[out_BR$beta/out_BR$se < qnorm(5e-8/2),]

summary_data_bmi_2 <- dplyr::arrange(summary_data_bmi_2,dplyr::desc(abs(summary_data_bmi_2$beta/summary_data_bmi_2$se)))
summary_data_sig_bmi_2 <- summary_data_bmi_2[summary_data_bmi_2$beta/summary_data_bmi_2$se < qnorm(5e-8/2),]
true_beta <- summary_data_sig_bmi_2$beta_rep

# Compute average bias for each method for 5e-8 threshold (negative)
bias_BMI2_down <- data.frame(naive = mean((out_EB_sig$beta - true_beta)), CL1 = mean((out_CL_sig$beta.cl1 - true_beta)), CL2 = mean((out_CL_sig$beta.cl2 - true_beta)), CL3 = mean((out_CL_sig$beta.cl3 - true_beta)), EB = mean((out_EB_sig$beta_EB - true_beta)),  EB_df = mean((out_EB_df_sig$beta_EB - true_beta)), EB_scam = mean((out_EB_scam_sig$beta_EB - true_beta)), EB_gam_po = mean((out_EB_gampo_sig$beta_EB - true_beta)), EB_gam_nb = mean((out_EB_gamnb_sig$beta_EB - true_beta)), boot = mean((out_BR_ss_sig$beta_BR_ss - true_beta)), FIQT = mean((out_FIQT_sig$beta_FIQT - true_beta)))


# Use threshold 5e-4 (need to re-apply conditional likelihood with change of threshold)

out_CL <- conditional_likelihood(summary_data = summary_stats, alpha=5e-4)

# Use threshold 5e-4 (positive)
out_CL_sig <- out_CL[out_CL$beta/out_CL$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_EB_sig <- out_EB[out_EB$beta/out_EB$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_EB_df_sig <- out_EB_df[out_EB_df$beta/out_EB_df$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_EB_scam_sig <- out_EB_scam[out_EB_scam$beta/out_EB_scam$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_EB_gampo_sig <- out_EB_gam_po[out_EB_gam_po$beta/out_EB_gam_po$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_EB_gamnb_sig <- out_EB_gam_nb[out_EB_gam_nb$beta/out_EB_gam_nb$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_FIQT_sig <- out_FIQT[out_FIQT$beta/out_FIQT$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_BR_ss_sig <- out_BR[out_BR$beta/out_BR$se > qnorm(5e-4/2, lower.tail=FALSE),]

summary_data_bmi_2 <- dplyr::arrange(summary_data_bmi_2,dplyr::desc(abs(summary_data_bmi_2$beta/summary_data_bmi_2$se)))
summary_data_sig_bmi_2 <- summary_data_bmi_2[summary_data_bmi_2$beta/summary_data_bmi_2$se > qnorm(5e-4/2, lower.tail=FALSE),]
true_beta <- summary_data_sig_bmi_2$beta_rep

# Compute average bias for each method for 5e-4 threshold (positive)
bias_BMI2_5e_4_up <- data.frame(naive = mean((out_EB_sig$beta - true_beta)), CL1 = mean((out_CL_sig$beta.cl1 - true_beta)), CL2 = mean((out_CL_sig$beta.cl2 - true_beta)), CL3 = mean((out_CL_sig$beta.cl3 - true_beta)), EB = mean((out_EB_sig$beta_EB - true_beta)),  EB_df = mean((out_EB_df_sig$beta_EB - true_beta)), EB_scam = mean((out_EB_scam_sig$beta_EB - true_beta)), EB_gam_po = mean((out_EB_gampo_sig$beta_EB - true_beta)), EB_gam_nb = mean((out_EB_gamnb_sig$beta_EB - true_beta)), boot = mean((out_BR_ss_sig$beta_BR_ss - true_beta)), FIQT = mean((out_FIQT_sig$beta_FIQT - true_beta)))

# Use threshold 5e-4 (negative)
out_CL_sig <- out_CL[out_CL$beta/out_CL$se < qnorm(5e-4/2,),]
out_EB_sig <- out_EB[out_EB$beta/out_EB$se < qnorm(5e-4/2),]
out_EB_df_sig <- out_EB_df[out_EB_df$beta/out_EB_df$se < qnorm(5e-4/2),]
out_EB_scam_sig <- out_EB_scam[out_EB_scam$beta/out_EB_scam$se < qnorm(5e-4/2),]
out_EB_gampo_sig <- out_EB_gam_po[out_EB_gam_po$beta/out_EB_gam_po$se < qnorm(5e-4/2),]
out_EB_gamnb_sig <- out_EB_gam_nb[out_EB_gam_nb$beta/out_EB_gam_nb$se < qnorm(5e-4/2),]
out_FIQT_sig <- out_FIQT[out_FIQT$beta/out_FIQT$se < qnorm(5e-4/2),]
out_BR_ss_sig <- out_BR[out_BR$beta/out_BR$se < qnorm(5e-4/2),]

summary_data_bmi_2 <- dplyr::arrange(summary_data_bmi_2,dplyr::desc(abs(summary_data_bmi_2$beta/summary_data_bmi_2$se)))
summary_data_sig_bmi_2 <- summary_data_bmi_2[summary_data_bmi_2$beta/summary_data_bmi_2$se < qnorm(5e-4/2),]
true_beta <- summary_data_sig_bmi_2$beta_rep

# Compute average bias for each method for 5e-4 threshold (negative)
bias_BMI2_5e_4_down <- data.frame(naive = mean((out_EB_sig$beta - true_beta)), CL1 = mean((out_CL_sig$beta.cl1 - true_beta)), CL2 = mean((out_CL_sig$beta.cl2 - true_beta)), CL3 = mean((out_CL_sig$beta.cl3 - true_beta)), EB = mean((out_EB_sig$beta_EB - true_beta)),  EB_df = mean((out_EB_df_sig$beta_EB - true_beta)), EB_scam = mean((out_EB_scam_sig$beta_EB - true_beta)), EB_gam_po = mean((out_EB_gampo_sig$beta_EB - true_beta)), EB_gam_nb = mean((out_EB_gamnb_sig$beta_EB - true_beta)), boot = mean((out_BR_ss_sig$beta_BR_ss - true_beta)), FIQT = mean((out_FIQT_sig$beta_FIQT - true_beta)))

###############################################################################

## T2D 1

summary_data_T2D_1 <-  read.table('data/summary_data_T2D_1.txt',header=TRUE)
summary_data_T2D_1$beta <- log(summary_data_T2D_1$beta)
summary_data_T2D_1$beta_rep <- log(summary_data_T2D_1$beta_rep)

# Apply methods
summary_stats <- summary_data_T2D_1[,3:5]
out_CL <- conditional_likelihood(summary_data = summary_stats)
out_EB <- empirical_bayes(summary_data = summary_stats)
out_EB_df <- empirical_bayes(summary_data=summary_stats, method="fix_df")
out_EB_scam <- empirical_bayes(summary_data=summary_stats,method="scam")
out_EB_gam_po <- empirical_bayes(summary_data=summary_stats,method="gam_po")
out_EB_gam_nb <- empirical_bayes(summary_data=summary_stats,method="gam_nb")
out_FIQT <- FDR_IQT(summary_data = summary_stats)
out_BR <- BR_ss(summary_data = summary_stats)

# Use threshold 5e-8 (positive)
out_CL_sig <- out_CL[out_CL$beta/out_CL$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_EB_sig <- out_EB[out_EB$beta/out_EB$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_EB_df_sig <- out_EB_df[out_EB_df$beta/out_EB_df$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_EB_scam_sig <- out_EB_scam[out_EB_scam$beta/out_EB_scam$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_EB_gampo_sig <- out_EB_gam_po[out_EB_gam_po$beta/out_EB_gam_po$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_EB_gamnb_sig <- out_EB_gam_nb[out_EB_gam_nb$beta/out_EB_gam_nb$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_FIQT_sig <- out_FIQT[out_FIQT$beta/out_FIQT$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_BR_ss_sig <- out_BR[out_BR$beta/out_BR$se > qnorm(5e-8/2, lower.tail=FALSE),]

summary_data_T2D_1 <- dplyr::arrange(summary_data_T2D_1,dplyr::desc(abs(summary_data_T2D_1$beta/summary_data_T2D_1$se)))
summary_data_sig_T2D_1 <- summary_data_T2D_1[summary_data_T2D_1$beta/summary_data_T2D_1$se > qnorm(5e-8/2, lower.tail=FALSE),]
true_beta <- summary_data_sig_T2D_1$beta_rep

# Compute average bias for each method for 5e-8 threshold (positive)
bias_T2D1_up <- data.frame(naive = mean((out_EB_sig$beta - true_beta)), CL1 = mean((out_CL_sig$beta.cl1 - true_beta)), CL2 = mean((out_CL_sig$beta.cl2 - true_beta)), CL3 = mean((out_CL_sig$beta.cl3 - true_beta)), EB = mean((out_EB_sig$beta_EB - true_beta)),  EB_df = mean((out_EB_df_sig$beta_EB - true_beta)), EB_scam = mean((out_EB_scam_sig$beta_EB - true_beta)), EB_gam_po = mean((out_EB_gampo_sig$beta_EB - true_beta)), EB_gam_nb = mean((out_EB_gamnb_sig$beta_EB - true_beta)), boot = mean((out_BR_ss_sig$beta_BR_ss - true_beta)), FIQT = mean((out_FIQT_sig$beta_FIQT - true_beta)))


# Use threshold 5e-8 (negative)
out_CL_sig <- out_CL[out_CL$beta/out_CL$se < qnorm(5e-8/2,),]
out_EB_sig <- out_EB[out_EB$beta/out_EB$se < qnorm(5e-8/2),]
out_EB_df_sig <- out_EB_df[out_EB_df$beta/out_EB_df$se < qnorm(5e-8/2),]
out_EB_scam_sig <- out_EB_scam[out_EB_scam$beta/out_EB_scam$se < qnorm(5e-8/2),]
out_EB_gampo_sig <- out_EB_gam_po[out_EB_gam_po$beta/out_EB_gam_po$se < qnorm(5e-8/2),]
out_EB_gamnb_sig <- out_EB_gam_nb[out_EB_gam_nb$beta/out_EB_gam_nb$se < qnorm(5e-8/2),]
out_FIQT_sig <- out_FIQT[out_FIQT$beta/out_FIQT$se < qnorm(5e-8/2),]
out_BR_ss_sig <- out_BR[out_BR$beta/out_BR$se < qnorm(5e-8/2),]

summary_data_T2D_1 <- dplyr::arrange(summary_data_T2D_1,dplyr::desc(abs(summary_data_T2D_1$beta/summary_data_T2D_1$se)))
summary_data_sig_T2D_1 <- summary_data_T2D_1[summary_data_T2D_1$beta/summary_data_T2D_1$se < qnorm(5e-8/2),]
true_beta <- summary_data_sig_T2D_1$beta_rep

# Compute average bias for each method for 5e-8 threshold (negative)
bias_T2D1_down <- data.frame(naive = mean((out_EB_sig$beta - true_beta)), CL1 = mean((out_CL_sig$beta.cl1 - true_beta)), CL2 = mean((out_CL_sig$beta.cl2 - true_beta)), CL3 = mean((out_CL_sig$beta.cl3 - true_beta)), EB = mean((out_EB_sig$beta_EB - true_beta)),  EB_df = mean((out_EB_df_sig$beta_EB - true_beta)), EB_scam = mean((out_EB_scam_sig$beta_EB - true_beta)), EB_gam_po = mean((out_EB_gampo_sig$beta_EB - true_beta)), EB_gam_nb = mean((out_EB_gamnb_sig$beta_EB - true_beta)), boot = mean((out_BR_ss_sig$beta_BR_ss - true_beta)), FIQT = mean((out_FIQT_sig$beta_FIQT - true_beta)))


# Use threshold 5e-4 (need to re-apply conditional likelihood with change of threshold)

out_CL <- conditional_likelihood(summary_data = summary_stats, alpha=5e-4)

# Use threshold 5e-4 (positive)
out_CL_sig <- out_CL[out_CL$beta/out_CL$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_EB_sig <- out_EB[out_EB$beta/out_EB$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_EB_df_sig <- out_EB_df[out_EB_df$beta/out_EB_df$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_EB_scam_sig <- out_EB_scam[out_EB_scam$beta/out_EB_scam$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_EB_gampo_sig <- out_EB_gam_po[out_EB_gam_po$beta/out_EB_gam_po$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_EB_gamnb_sig <- out_EB_gam_nb[out_EB_gam_nb$beta/out_EB_gam_nb$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_FIQT_sig <- out_FIQT[out_FIQT$beta/out_FIQT$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_BR_ss_sig <- out_BR[out_BR$beta/out_BR$se > qnorm(5e-4/2, lower.tail=FALSE),]

summary_data_T2D_1 <- dplyr::arrange(summary_data_T2D_1,dplyr::desc(abs(summary_data_T2D_1$beta/summary_data_T2D_1$se)))
summary_data_sig_T2D_1 <- summary_data_T2D_1[summary_data_T2D_1$beta/summary_data_T2D_1$se > qnorm(5e-4/2, lower.tail=FALSE),]
true_beta <- summary_data_sig_T2D_1$beta_rep

# Compute average bias for each method for 5e-4 threshold (positive)
bias_T2D1_5e_4_up <- data.frame(naive = mean((out_EB_sig$beta - true_beta)), CL1 = mean((out_CL_sig$beta.cl1 - true_beta)), CL2 = mean((out_CL_sig$beta.cl2 - true_beta)), CL3 = mean((out_CL_sig$beta.cl3 - true_beta)), EB = mean((out_EB_sig$beta_EB - true_beta)),  EB_df = mean((out_EB_df_sig$beta_EB - true_beta)), EB_scam = mean((out_EB_scam_sig$beta_EB - true_beta)), EB_gam_po = mean((out_EB_gampo_sig$beta_EB - true_beta)), EB_gam_nb = mean((out_EB_gamnb_sig$beta_EB - true_beta)), boot = mean((out_BR_ss_sig$beta_BR_ss - true_beta)), FIQT = mean((out_FIQT_sig$beta_FIQT - true_beta)))


# Use threshold 5e-4 (negative)
out_CL_sig <- out_CL[out_CL$beta/out_CL$se < qnorm(5e-4/2,),]
out_EB_sig <- out_EB[out_EB$beta/out_EB$se < qnorm(5e-4/2),]
out_EB_df_sig <- out_EB_df[out_EB_df$beta/out_EB_df$se < qnorm(5e-4/2),]
out_EB_scam_sig <- out_EB_scam[out_EB_scam$beta/out_EB_scam$se < qnorm(5e-4/2),]
out_EB_gampo_sig <- out_EB_gam_po[out_EB_gam_po$beta/out_EB_gam_po$se < qnorm(5e-4/2),]
out_EB_gamnb_sig <- out_EB_gam_nb[out_EB_gam_nb$beta/out_EB_gam_nb$se < qnorm(5e-4/2),]
out_FIQT_sig <- out_FIQT[out_FIQT$beta/out_FIQT$se < qnorm(5e-4/2),]
out_BR_ss_sig <- out_BR[out_BR$beta/out_BR$se < qnorm(5e-4/2),]

summary_data_T2D_1 <- dplyr::arrange(summary_data_T2D_1,dplyr::desc(abs(summary_data_T2D_1$beta/summary_data_T2D_1$se)))
summary_data_sig_T2D_1 <- summary_data_T2D_1[summary_data_T2D_1$beta/summary_data_T2D_1$se < qnorm(5e-4/2),]
true_beta <- summary_data_sig_T2D_1$beta_rep

# Compute average bias for each method for 5e-4 threshold (negative)
bias_T2D1_5e_4_down <- data.frame(naive = mean((out_EB_sig$beta - true_beta)), CL1 = mean((out_CL_sig$beta.cl1 - true_beta)), CL2 = mean((out_CL_sig$beta.cl2 - true_beta)), CL3 = mean((out_CL_sig$beta.cl3 - true_beta)), EB = mean((out_EB_sig$beta_EB - true_beta)),  EB_df = mean((out_EB_df_sig$beta_EB - true_beta)), EB_scam = mean((out_EB_scam_sig$beta_EB - true_beta)), EB_gam_po = mean((out_EB_gampo_sig$beta_EB - true_beta)), EB_gam_nb = mean((out_EB_gamnb_sig$beta_EB - true_beta)), boot = mean((out_BR_ss_sig$beta_BR_ss - true_beta)), FIQT = mean((out_FIQT_sig$beta_FIQT - true_beta)))

###############################################################################

## T2D 2

summary_data_T2D_2 <-  read.table('data/summary_data_T2D_2.txt',header=TRUE)
summary_data_T2D_2$beta <- log(summary_data_T2D_2$beta)
summary_data_T2D_2$beta_rep <- log(summary_data_T2D_2$beta_rep)

# Apply methods
summary_stats <- summary_data_T2D_2[,3:5]
out_CL <- conditional_likelihood(summary_data = summary_stats)
out_EB <- empirical_bayes(summary_data = summary_stats)
out_EB_df <- empirical_bayes(summary_data=summary_stats, method="fix_df")
out_EB_scam <- empirical_bayes(summary_data=summary_stats,method="scam")
out_EB_gam_po <- empirical_bayes(summary_data=summary_stats,method="gam_po")
out_EB_gam_nb <- empirical_bayes(summary_data=summary_stats,method="gam_nb")
out_FIQT <- FDR_IQT(summary_data = summary_stats)
out_BR <- BR_ss(summary_data = summary_stats)

# Use threshold 5e-8 (positive)
out_CL_sig <- out_CL[out_CL$beta/out_CL$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_EB_sig <- out_EB[out_EB$beta/out_EB$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_EB_df_sig <- out_EB_df[out_EB_df$beta/out_EB_df$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_EB_scam_sig <- out_EB_scam[out_EB_scam$beta/out_EB_scam$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_EB_gampo_sig <- out_EB_gam_po[out_EB_gam_po$beta/out_EB_gam_po$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_EB_gamnb_sig <- out_EB_gam_nb[out_EB_gam_nb$beta/out_EB_gam_nb$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_FIQT_sig <- out_FIQT[out_FIQT$beta/out_FIQT$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_BR_ss_sig <- out_BR[out_BR$beta/out_BR$se > qnorm(5e-8/2, lower.tail=FALSE),]

summary_data_T2D_2 <- dplyr::arrange(summary_data_T2D_2,dplyr::desc(abs(summary_data_T2D_2$beta/summary_data_T2D_2$se)))
summary_data_sig_T2D_2 <- summary_data_T2D_2[summary_data_T2D_2$beta/summary_data_T2D_2$se > qnorm(5e-8/2, lower.tail=FALSE),]
true_beta <- summary_data_sig_T2D_2$beta_rep

# Compute average bias for each method for 5e-8 threshold (positive)
bias_T2D2_up <- data.frame(naive = mean((out_EB_sig$beta - true_beta)), CL1 = mean((out_CL_sig$beta.cl1 - true_beta)), CL2 = mean((out_CL_sig$beta.cl2 - true_beta)), CL3 = mean((out_CL_sig$beta.cl3 - true_beta)), EB = mean((out_EB_sig$beta_EB - true_beta)),  EB_df = mean((out_EB_df_sig$beta_EB - true_beta)), EB_scam = mean((out_EB_scam_sig$beta_EB - true_beta)), EB_gam_po = mean((out_EB_gampo_sig$beta_EB - true_beta)), EB_gam_nb = mean((out_EB_gamnb_sig$beta_EB - true_beta)), boot = mean((out_BR_ss_sig$beta_BR_ss - true_beta)), FIQT = mean((out_FIQT_sig$beta_FIQT - true_beta)))


# Use threshold 5e-8 (negative)
out_CL_sig <- out_CL[out_CL$beta/out_CL$se < qnorm(5e-8/2,),]
out_EB_sig <- out_EB[out_EB$beta/out_EB$se < qnorm(5e-8/2),]
out_EB_df_sig <- out_EB_df[out_EB_df$beta/out_EB_df$se < qnorm(5e-8/2),]
out_EB_scam_sig <- out_EB_scam[out_EB_scam$beta/out_EB_scam$se < qnorm(5e-8/2),]
out_EB_gampo_sig <- out_EB_gam_po[out_EB_gam_po$beta/out_EB_gam_po$se < qnorm(5e-8/2),]
out_EB_gamnb_sig <- out_EB_gam_nb[out_EB_gam_nb$beta/out_EB_gam_nb$se < qnorm(5e-8/2),]
out_FIQT_sig <- out_FIQT[out_FIQT$beta/out_FIQT$se < qnorm(5e-8/2),]
out_BR_ss_sig <- out_BR[out_BR$beta/out_BR$se < qnorm(5e-8/2),]

summary_data_T2D_2 <- dplyr::arrange(summary_data_T2D_2,dplyr::desc(abs(summary_data_T2D_2$beta/summary_data_T2D_2$se)))
summary_data_sig_T2D_2 <- summary_data_T2D_2[summary_data_T2D_2$beta/summary_data_T2D_2$se < qnorm(5e-8/2),]
true_beta <- summary_data_sig_T2D_2$beta_rep

# Compute average bias for each method for 5e-8 threshold (negative)
bias_T2D2_down <- data.frame(naive = mean((out_EB_sig$beta - true_beta)), CL1 = mean((out_CL_sig$beta.cl1 - true_beta)), CL2 = mean((out_CL_sig$beta.cl2 - true_beta)), CL3 = mean((out_CL_sig$beta.cl3 - true_beta)), EB = mean((out_EB_sig$beta_EB - true_beta)),  EB_df = mean((out_EB_df_sig$beta_EB - true_beta)), EB_scam = mean((out_EB_scam_sig$beta_EB - true_beta)), EB_gam_po = mean((out_EB_gampo_sig$beta_EB - true_beta)), EB_gam_nb = mean((out_EB_gamnb_sig$beta_EB - true_beta)), boot = mean((out_BR_ss_sig$beta_BR_ss - true_beta)), FIQT = mean((out_FIQT_sig$beta_FIQT - true_beta)))


# Use threshold 5e-4 (need to re-apply conditional likelihood with change of threshold)

out_CL <- conditional_likelihood(summary_data = summary_stats, alpha=5e-4)

# Use threshold 5e-4 (positive)
out_CL_sig <- out_CL[out_CL$beta/out_CL$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_EB_sig <- out_EB[out_EB$beta/out_EB$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_EB_df_sig <- out_EB_df[out_EB_df$beta/out_EB_df$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_EB_scam_sig <- out_EB_scam[out_EB_scam$beta/out_EB_scam$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_EB_gampo_sig <- out_EB_gam_po[out_EB_gam_po$beta/out_EB_gam_po$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_EB_gamnb_sig <- out_EB_gam_nb[out_EB_gam_nb$beta/out_EB_gam_nb$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_FIQT_sig <- out_FIQT[out_FIQT$beta/out_FIQT$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_BR_ss_sig <- out_BR[out_BR$beta/out_BR$se > qnorm(5e-4/2, lower.tail=FALSE),]

summary_data_T2D_2 <- dplyr::arrange(summary_data_T2D_2,dplyr::desc(abs(summary_data_T2D_2$beta/summary_data_T2D_2$se)))
summary_data_sig_T2D_2 <- summary_data_T2D_2[summary_data_T2D_2$beta/summary_data_T2D_2$se > qnorm(5e-4/2, lower.tail=FALSE),]
true_beta <- summary_data_sig_T2D_2$beta_rep

# Compute average bias for each method for 5e-4 threshold (positive)
bias_T2D2_5e_4_up <- data.frame(naive = mean((out_EB_sig$beta - true_beta)), CL1 = mean((out_CL_sig$beta.cl1 - true_beta)), CL2 = mean((out_CL_sig$beta.cl2 - true_beta)), CL3 = mean((out_CL_sig$beta.cl3 - true_beta)), EB = mean((out_EB_sig$beta_EB - true_beta)),  EB_df = mean((out_EB_df_sig$beta_EB - true_beta)), EB_scam = mean((out_EB_scam_sig$beta_EB - true_beta)), EB_gam_po = mean((out_EB_gampo_sig$beta_EB - true_beta)), EB_gam_nb = mean((out_EB_gamnb_sig$beta_EB - true_beta)), boot = mean((out_BR_ss_sig$beta_BR_ss - true_beta)), FIQT = mean((out_FIQT_sig$beta_FIQT - true_beta)))


# Use threshold 5e-4 (negative)
out_CL_sig <- out_CL[out_CL$beta/out_CL$se < qnorm(5e-4/2,),]
out_EB_sig <- out_EB[out_EB$beta/out_EB$se < qnorm(5e-4/2),]
out_EB_df_sig <- out_EB_df[out_EB_df$beta/out_EB_df$se < qnorm(5e-4/2),]
out_EB_scam_sig <- out_EB_scam[out_EB_scam$beta/out_EB_scam$se < qnorm(5e-4/2),]
out_EB_gampo_sig <- out_EB_gam_po[out_EB_gam_po$beta/out_EB_gam_po$se < qnorm(5e-4/2),]
out_EB_gamnb_sig <- out_EB_gam_nb[out_EB_gam_nb$beta/out_EB_gam_nb$se < qnorm(5e-4/2),]
out_FIQT_sig <- out_FIQT[out_FIQT$beta/out_FIQT$se < qnorm(5e-4/2),]
out_BR_ss_sig <- out_BR[out_BR$beta/out_BR$se < qnorm(5e-4/2),]

summary_data_T2D_2 <- dplyr::arrange(summary_data_T2D_2,dplyr::desc(abs(summary_data_T2D_2$beta/summary_data_T2D_2$se)))
summary_data_sig_T2D_2 <- summary_data_T2D_2[summary_data_T2D_2$beta/summary_data_T2D_2$se < qnorm(5e-4/2),]
true_beta <- summary_data_sig_T2D_2$beta_rep

# Compute average bias for each method for 5e-4 threshold (negative)
bias_T2D2_5e_4_down <- data.frame(naive = mean((out_EB_sig$beta - true_beta)), CL1 = mean((out_CL_sig$beta.cl1 - true_beta)), CL2 = mean((out_CL_sig$beta.cl2 - true_beta)), CL3 = mean((out_CL_sig$beta.cl3 - true_beta)), EB = mean((out_EB_sig$beta_EB - true_beta)),  EB_df = mean((out_EB_df_sig$beta_EB - true_beta)), EB_scam = mean((out_EB_scam_sig$beta_EB - true_beta)), EB_gam_po = mean((out_EB_gampo_sig$beta_EB - true_beta)), EB_gam_nb = mean((out_EB_gamnb_sig$beta_EB - true_beta)), boot = mean((out_BR_ss_sig$beta_BR_ss - true_beta)), FIQT = mean((out_FIQT_sig$beta_FIQT - true_beta)))

###############################################################################

## Height 1

summary_data_height_1 <-  read.table('data/summary_data_height_1.txt',header=TRUE)

# Apply methods
summary_stats <- summary_data_height_1[,3:5]
out_CL <- conditional_likelihood(summary_data = summary_stats)
out_EB <- empirical_bayes(summary_data = summary_stats)
out_EB_df <- empirical_bayes(summary_data=summary_stats, method="fix_df")
out_EB_scam <- empirical_bayes(summary_data=summary_stats,method="scam")
out_EB_gam_po <- empirical_bayes(summary_data=summary_stats,method="gam_po")
out_EB_gam_nb <- empirical_bayes(summary_data=summary_stats,method="gam_nb")
out_FIQT <- FDR_IQT(summary_data = summary_stats)
out_BR <- BR_ss(summary_data = summary_stats)

# Use threshold 5e-8 (positive)
out_CL_sig <- out_CL[out_CL$beta/out_CL$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_EB_sig <- out_EB[out_EB$beta/out_EB$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_EB_df_sig <- out_EB_df[out_EB_df$beta/out_EB_df$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_EB_scam_sig <- out_EB_scam[out_EB_scam$beta/out_EB_scam$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_EB_gampo_sig <- out_EB_gam_po[out_EB_gam_po$beta/out_EB_gam_po$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_EB_gamnb_sig <- out_EB_gam_nb[out_EB_gam_nb$beta/out_EB_gam_nb$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_FIQT_sig <- out_FIQT[out_FIQT$beta/out_FIQT$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_BR_ss_sig <- out_BR[out_BR$beta/out_BR$se > qnorm(5e-8/2, lower.tail=FALSE),]

summary_data_height_1 <- dplyr::arrange(summary_data_height_1,dplyr::desc(abs(summary_data_height_1$beta/summary_data_height_1$se)))
summary_data_sig_height_1 <- summary_data_height_1[summary_data_height_1$beta/summary_data_height_1$se > qnorm(5e-8/2, lower.tail=FALSE),]
true_beta <- summary_data_sig_height_1$beta_rep

# Compute average bias for each method for 5e-8 threshold (positive)
bias_height1_up <- data.frame(naive = mean((out_EB_sig$beta - true_beta)), CL1 = mean((out_CL_sig$beta.cl1 - true_beta)), CL2 = mean((out_CL_sig$beta.cl2 - true_beta)), CL3 = mean((out_CL_sig$beta.cl3 - true_beta)), EB = mean((out_EB_sig$beta_EB - true_beta)),  EB_df = mean((out_EB_df_sig$beta_EB - true_beta)), EB_scam = mean((out_EB_scam_sig$beta_EB - true_beta)), EB_gam_po = mean((out_EB_gampo_sig$beta_EB - true_beta)), EB_gam_nb = mean((out_EB_gamnb_sig$beta_EB - true_beta)), boot = mean((out_BR_ss_sig$beta_BR_ss - true_beta)), FIQT = mean((out_FIQT_sig$beta_FIQT - true_beta)))

# Use threshold 5e-8 (negative)
out_CL_sig <- out_CL[out_CL$beta/out_CL$se < qnorm(5e-8/2,),]
out_EB_sig <- out_EB[out_EB$beta/out_EB$se < qnorm(5e-8/2),]
out_EB_df_sig <- out_EB_df[out_EB_df$beta/out_EB_df$se < qnorm(5e-8/2),]
out_EB_scam_sig <- out_EB_scam[out_EB_scam$beta/out_EB_scam$se < qnorm(5e-8/2),]
out_EB_gampo_sig <- out_EB_gam_po[out_EB_gam_po$beta/out_EB_gam_po$se < qnorm(5e-8/2),]
out_EB_gamnb_sig <- out_EB_gam_nb[out_EB_gam_nb$beta/out_EB_gam_nb$se < qnorm(5e-8/2),]
out_FIQT_sig <- out_FIQT[out_FIQT$beta/out_FIQT$se < qnorm(5e-8/2),]
out_BR_ss_sig <- out_BR[out_BR$beta/out_BR$se < qnorm(5e-8/2),]

summary_data_height_1 <- dplyr::arrange(summary_data_height_1,dplyr::desc(abs(summary_data_height_1$beta/summary_data_height_1$se)))
summary_data_sig_height_1 <- summary_data_height_1[summary_data_height_1$beta/summary_data_height_1$se < qnorm(5e-8/2),]
true_beta <- summary_data_sig_height_1$beta_rep

# Compute average bias for each method for 5e-8 threshold (negative)
bias_height1_down <- data.frame(naive = mean((out_EB_sig$beta - true_beta)), CL1 = mean((out_CL_sig$beta.cl1 - true_beta)), CL2 = mean((out_CL_sig$beta.cl2 - true_beta)), CL3 = mean((out_CL_sig$beta.cl3 - true_beta)), EB = mean((out_EB_sig$beta_EB - true_beta)),  EB_df = mean((out_EB_df_sig$beta_EB - true_beta)), EB_scam = mean((out_EB_scam_sig$beta_EB - true_beta)), EB_gam_po = mean((out_EB_gampo_sig$beta_EB - true_beta)), EB_gam_nb = mean((out_EB_gamnb_sig$beta_EB - true_beta)), boot = mean((out_BR_ss_sig$beta_BR_ss - true_beta)), FIQT = mean((out_FIQT_sig$beta_FIQT - true_beta)))


# Use threshold 5e-4 (need to re-apply conditional likelihood with change of threshold)

out_CL <- conditional_likelihood(summary_data = summary_stats, alpha=5e-4)

# Use threshold 5e-4 (positive)
out_CL_sig <- out_CL[out_CL$beta/out_CL$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_EB_sig <- out_EB[out_EB$beta/out_EB$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_EB_df_sig <- out_EB_df[out_EB_df$beta/out_EB_df$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_EB_scam_sig <- out_EB_scam[out_EB_scam$beta/out_EB_scam$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_EB_gampo_sig <- out_EB_gam_po[out_EB_gam_po$beta/out_EB_gam_po$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_EB_gamnb_sig <- out_EB_gam_nb[out_EB_gam_nb$beta/out_EB_gam_nb$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_FIQT_sig <- out_FIQT[out_FIQT$beta/out_FIQT$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_BR_ss_sig <- out_BR[out_BR$beta/out_BR$se > qnorm(5e-4/2, lower.tail=FALSE),]

summary_data_height_1 <- dplyr::arrange(summary_data_height_1,dplyr::desc(abs(summary_data_height_1$beta/summary_data_height_1$se)))
summary_data_sig_height_1 <- summary_data_height_1[summary_data_height_1$beta/summary_data_height_1$se > qnorm(5e-4/2, lower.tail=FALSE),]
true_beta <- summary_data_sig_height_1$beta_rep

# Compute average bias for each method for 5e-4 threshold (positive)
bias_height1_5e_4_up <- data.frame(naive = mean((out_EB_sig$beta - true_beta)), CL1 = mean((out_CL_sig$beta.cl1 - true_beta)), CL2 = mean((out_CL_sig$beta.cl2 - true_beta)), CL3 = mean((out_CL_sig$beta.cl3 - true_beta)), EB = mean((out_EB_sig$beta_EB - true_beta)),  EB_df = mean((out_EB_df_sig$beta_EB - true_beta)), EB_scam = mean((out_EB_scam_sig$beta_EB - true_beta)), EB_gam_po = mean((out_EB_gampo_sig$beta_EB - true_beta)), EB_gam_nb = mean((out_EB_gamnb_sig$beta_EB - true_beta)), boot = mean((out_BR_ss_sig$beta_BR_ss - true_beta)), FIQT = mean((out_FIQT_sig$beta_FIQT - true_beta)))

# Use threshold 5e-4 (negative)
out_CL_sig <- out_CL[out_CL$beta/out_CL$se < qnorm(5e-4/2,),]
out_EB_sig <- out_EB[out_EB$beta/out_EB$se < qnorm(5e-4/2),]
out_EB_df_sig <- out_EB_df[out_EB_df$beta/out_EB_df$se < qnorm(5e-4/2),]
out_EB_scam_sig <- out_EB_scam[out_EB_scam$beta/out_EB_scam$se < qnorm(5e-4/2),]
out_EB_gampo_sig <- out_EB_gam_po[out_EB_gam_po$beta/out_EB_gam_po$se < qnorm(5e-4/2),]
out_EB_gamnb_sig <- out_EB_gam_nb[out_EB_gam_nb$beta/out_EB_gam_nb$se < qnorm(5e-4/2),]
out_FIQT_sig <- out_FIQT[out_FIQT$beta/out_FIQT$se < qnorm(5e-4/2),]
out_BR_ss_sig <- out_BR[out_BR$beta/out_BR$se < qnorm(5e-4/2),]

summary_data_height_1 <- dplyr::arrange(summary_data_height_1,dplyr::desc(abs(summary_data_height_1$beta/summary_data_height_1$se)))
summary_data_sig_height_1 <- summary_data_T2D_1[summary_data_height_1$beta/summary_data_height_1$se < qnorm(5e-4/2),]
true_beta <- summary_data_sig_height_1$beta_rep

# Compute average bias for each method for 5e-4 threshold (negative)
bias_height1_5e_4_down <- data.frame(naive = mean((out_EB_sig$beta - true_beta)), CL1 = mean((out_CL_sig$beta.cl1 - true_beta)), CL2 = mean((out_CL_sig$beta.cl2 - true_beta)), CL3 = mean((out_CL_sig$beta.cl3 - true_beta)), EB = mean((out_EB_sig$beta_EB - true_beta)),  EB_df = mean((out_EB_df_sig$beta_EB - true_beta)), EB_scam = mean((out_EB_scam_sig$beta_EB - true_beta)), EB_gam_po = mean((out_EB_gampo_sig$beta_EB - true_beta)), EB_gam_nb = mean((out_EB_gamnb_sig$beta_EB - true_beta)), boot = mean((out_BR_ss_sig$beta_BR_ss - true_beta)), FIQT = mean((out_FIQT_sig$beta_FIQT - true_beta)))

###############################################################################

## Height 2

summary_data_height_2 <-  read.table('data/summary_data_height_2.txt',header=TRUE)

# Apply methods
summary_stats <- summary_data_height_2[,3:5]
out_CL <- conditional_likelihood(summary_data = summary_stats)
out_EB <- empirical_bayes(summary_data = summary_stats)
out_EB_df <- empirical_bayes(summary_data=summary_stats, method="fix_df")
out_EB_scam <- empirical_bayes(summary_data=summary_stats,method="scam")
out_EB_gam_po <- empirical_bayes(summary_data=summary_stats,method="gam_po")
out_EB_gam_nb <- empirical_bayes(summary_data=summary_stats,method="gam_nb")
out_FIQT <- FDR_IQT(summary_data = summary_stats)
out_BR <- BR_ss(summary_data = summary_stats)

# Use threshold 5e-8 (positive)
out_CL_sig <- out_CL[out_CL$beta/out_CL$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_EB_sig <- out_EB[out_EB$beta/out_EB$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_EB_df_sig <- out_EB_df[out_EB_df$beta/out_EB_df$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_EB_scam_sig <- out_EB_scam[out_EB_scam$beta/out_EB_scam$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_EB_gampo_sig <- out_EB_gam_po[out_EB_gam_po$beta/out_EB_gam_po$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_EB_gamnb_sig <- out_EB_gam_nb[out_EB_gam_nb$beta/out_EB_gam_nb$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_FIQT_sig <- out_FIQT[out_FIQT$beta/out_FIQT$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_BR_ss_sig <- out_BR[out_BR$beta/out_BR$se > qnorm(5e-8/2, lower.tail=FALSE),]

summary_data_height_2 <- dplyr::arrange(summary_data_height_2,dplyr::desc(abs(summary_data_height_2$beta/summary_data_height_2$se)))
summary_data_sig_height_2 <- summary_data_height_2[summary_data_height_2$beta/summary_data_height_2$se > qnorm(5e-8/2, lower.tail=FALSE),]
true_beta <- summary_data_sig_height_2$beta_rep

# Compute average bias for each method for 5e-8 threshold (positive)
bias_height2_up <- data.frame(naive = mean((out_EB_sig$beta - true_beta)), CL1 = mean((out_CL_sig$beta.cl1 - true_beta)), CL2 = mean((out_CL_sig$beta.cl2 - true_beta)), CL3 = mean((out_CL_sig$beta.cl3 - true_beta)), EB = mean((out_EB_sig$beta_EB - true_beta)),  EB_df = mean((out_EB_df_sig$beta_EB - true_beta)), EB_scam = mean((out_EB_scam_sig$beta_EB - true_beta)), EB_gam_po = mean((out_EB_gampo_sig$beta_EB - true_beta)), EB_gam_nb = mean((out_EB_gamnb_sig$beta_EB - true_beta)), boot = mean((out_BR_ss_sig$beta_BR_ss - true_beta)), FIQT = mean((out_FIQT_sig$beta_FIQT - true_beta)))

# Use threshold 5e-8 (negative)
out_CL_sig <- out_CL[out_CL$beta/out_CL$se < qnorm(5e-8/2,),]
out_EB_sig <- out_EB[out_EB$beta/out_EB$se < qnorm(5e-8/2),]
out_EB_df_sig <- out_EB_df[out_EB_df$beta/out_EB_df$se < qnorm(5e-8/2),]
out_EB_scam_sig <- out_EB_scam[out_EB_scam$beta/out_EB_scam$se < qnorm(5e-8/2),]
out_EB_gampo_sig <- out_EB_gam_po[out_EB_gam_po$beta/out_EB_gam_po$se < qnorm(5e-8/2),]
out_EB_gamnb_sig <- out_EB_gam_nb[out_EB_gam_nb$beta/out_EB_gam_nb$se < qnorm(5e-8/2),]
out_FIQT_sig <- out_FIQT[out_FIQT$beta/out_FIQT$se < qnorm(5e-8/2),]
out_BR_ss_sig <- out_BR[out_BR$beta/out_BR$se < qnorm(5e-8/2),]

summary_data_height_2 <- dplyr::arrange(summary_data_height_2,dplyr::desc(abs(summary_data_height_2$beta/summary_data_height_2$se)))
summary_data_sig_height_2 <- summary_data_height_2[summary_data_height_2$beta/summary_data_height_2$se < qnorm(5e-8/2),]
true_beta <- summary_data_sig_height_2$beta_rep

# Compute average bias for each method for 5e-8 threshold (negative)
bias_height2_down <- data.frame(naive = mean((out_EB_sig$beta - true_beta)), CL1 = mean((out_CL_sig$beta.cl1 - true_beta)), CL2 = mean((out_CL_sig$beta.cl2 - true_beta)), CL3 = mean((out_CL_sig$beta.cl3 - true_beta)), EB = mean((out_EB_sig$beta_EB - true_beta)),  EB_df = mean((out_EB_df_sig$beta_EB - true_beta)), EB_scam = mean((out_EB_scam_sig$beta_EB - true_beta)), EB_gam_po = mean((out_EB_gampo_sig$beta_EB - true_beta)), EB_gam_nb = mean((out_EB_gamnb_sig$beta_EB - true_beta)), boot = mean((out_BR_ss_sig$beta_BR_ss - true_beta)), FIQT = mean((out_FIQT_sig$beta_FIQT - true_beta)))


# Use threshold 5e-4 (need to re-apply conditional likelihood with change of threshold)

out_CL <- conditional_likelihood(summary_data = summary_stats, alpha=5e-4)

# Use threshold 5e-4 (positive)
out_CL_sig <- out_CL[out_CL$beta/out_CL$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_EB_sig <- out_EB[out_EB$beta/out_EB$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_EB_df_sig <- out_EB_df[out_EB_df$beta/out_EB_df$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_EB_scam_sig <- out_EB_scam[out_EB_scam$beta/out_EB_scam$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_EB_gampo_sig <- out_EB_gam_po[out_EB_gam_po$beta/out_EB_gam_po$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_EB_gamnb_sig <- out_EB_gam_nb[out_EB_gam_nb$beta/out_EB_gam_nb$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_FIQT_sig <- out_FIQT[out_FIQT$beta/out_FIQT$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_BR_ss_sig <- out_BR[out_BR$beta/out_BR$se > qnorm(5e-4/2, lower.tail=FALSE),]

summary_data_height_2 <- dplyr::arrange(summary_data_height_2,dplyr::desc(abs(summary_data_height_2$beta/summary_data_height_2$se)))
summary_data_sig_height_2 <- summary_data_height_2[summary_data_height_2$beta/summary_data_height_2$se > qnorm(5e-4/2, lower.tail=FALSE),]
true_beta <- summary_data_sig_height_2$beta_rep

# Compute average bias for each method for 5e-4 threshold (positive)
bias_height2_5e_4_up <- data.frame(naive = mean((out_EB_sig$beta - true_beta)), CL1 = mean((out_CL_sig$beta.cl1 - true_beta)), CL2 = mean((out_CL_sig$beta.cl2 - true_beta)), CL3 = mean((out_CL_sig$beta.cl3 - true_beta)), EB = mean((out_EB_sig$beta_EB - true_beta)),  EB_df = mean((out_EB_df_sig$beta_EB - true_beta)), EB_scam = mean((out_EB_scam_sig$beta_EB - true_beta)), EB_gam_po = mean((out_EB_gampo_sig$beta_EB - true_beta)), EB_gam_nb = mean((out_EB_gamnb_sig$beta_EB - true_beta)), boot = mean((out_BR_ss_sig$beta_BR_ss - true_beta)), FIQT = mean((out_FIQT_sig$beta_FIQT - true_beta)))

# Use threshold 5e-4 (negative)
out_CL_sig <- out_CL[out_CL$beta/out_CL$se < qnorm(5e-4/2,),]
out_EB_sig <- out_EB[out_EB$beta/out_EB$se < qnorm(5e-4/2),]
out_EB_df_sig <- out_EB_df[out_EB_df$beta/out_EB_df$se < qnorm(5e-4/2),]
out_EB_scam_sig <- out_EB_scam[out_EB_scam$beta/out_EB_scam$se < qnorm(5e-4/2),]
out_EB_gampo_sig <- out_EB_gam_po[out_EB_gam_po$beta/out_EB_gam_po$se < qnorm(5e-4/2),]
out_EB_gamnb_sig <- out_EB_gam_nb[out_EB_gam_nb$beta/out_EB_gam_nb$se < qnorm(5e-4/2),]
out_FIQT_sig <- out_FIQT[out_FIQT$beta/out_FIQT$se < qnorm(5e-4/2),]
out_BR_ss_sig <- out_BR[out_BR$beta/out_BR$se < qnorm(5e-4/2),]

summary_data_height_2 <- dplyr::arrange(summary_data_height_2,dplyr::desc(abs(summary_data_height_2$beta/summary_data_height_2$se)))
summary_data_sig_height_2 <- summary_data_T2D_2[summary_data_height_2$beta/summary_data_height_2$se < qnorm(5e-4/2),]
true_beta <- summary_data_sig_height_2$beta_rep

# Compute average bias for each method for 5e-4 threshold (negative)
bias_height2_5e_4_down <- data.frame(naive = mean((out_EB_sig$beta - true_beta)), CL1 = mean((out_CL_sig$beta.cl1 - true_beta)), CL2 = mean((out_CL_sig$beta.cl2 - true_beta)), CL3 = mean((out_CL_sig$beta.cl3 - true_beta)), EB = mean((out_EB_sig$beta_EB - true_beta)),  EB_df = mean((out_EB_df_sig$beta_EB - true_beta)), EB_scam = mean((out_EB_scam_sig$beta_EB - true_beta)), EB_gam_po = mean((out_EB_gampo_sig$beta_EB - true_beta)), EB_gam_nb = mean((out_EB_gamnb_sig$beta_EB - true_beta)), boot = mean((out_BR_ss_sig$beta_BR_ss - true_beta)), FIQT = mean((out_FIQT_sig$beta_FIQT - true_beta)))

###############################################################################

## Combine all above results and save them in external results files 

bias_5e_8_up <- rbind(bias_BMI1_up,bias_BMI2_up, bias_T2D1_up, bias_T2D2_up, bias_height1_up, bias_height2_up)
bias_5e_8_up <- round(bias_5e_8_up[,1:11],5)
GWAS <- c("BMI 1", "BMI 2", "T2D 1", "T2D 2", "Height 1", "Height 2")
bias_5e_8_up <- cbind(GWAS,bias_5e_8_up)
write.csv(bias_5e_8_up,"bias_5e_8_positive.txt", row.names = FALSE)

bias_5e_8_down <- rbind(bias_BMI1_down,bias_BMI2_down, bias_T2D1_down, bias_T2D2_down, bias_height1_down, bias_height2_down)
bias_5e_8_down <- round(bias_5e_8_down[,1:11],5)
GWAS <- c("BMI 1", "BMI 2", "T2D 1", "T2D 2", "Height 1", "Height 2")
bias_5e_8_down <- cbind(GWAS,bias_5e_8_down)
write.csv(bias_5e_8_down,"bias_5e_8_negative.txt", row.names = FALSE)

bias_5e_4_up <- rbind(bias_BMI1_5e_4_up,bias_BMI2_5e_4_up, bias_T2D1_5e_4_up, bias_T2D2_5e_4_up, bias_height1_5e_4_up, bias_height2_5e_4_up)
bias_5e_4_up <- round(bias_5e_4_up[,1:11],5)
GWAS <- c("BMI 1", "BMI 2", "T2D 1", "T2D 2", "Height 1", "Height 2")
bias_5e_4_up <- cbind(GWAS,bias_5e_4_up)
write.csv(bias_5e_4_up,"bias_5e_4_positive.txt", row.names = FALSE)

bias_5e_4_down <- rbind(bias_BMI1_5e_4_down,bias_BMI2_5e_4_down, bias_T2D1_5e_4_down, bias_T2D2_5e_4_down, bias_height1_5e_4_down, bias_height2_5e_4_down)
bias_5e_4_down <- round(bias_5e_4_down[,1:11],5)
GWAS <- c("BMI 1", "BMI 2", "T2D 1", "T2D 2", "Height 1", "Height 2")
bias_5e_4_down <- cbind(GWAS,bias_5e_4_down)
write.csv(bias_5e_4_down,"bias_5e_4_negative.txt", row.names = FALSE)

###############################################################################
###############################################################################

## PART 5) PLOT BIAS RESULTS APPROPRIATELY


new_bias_5e_8_up <- pivot_longer(bias_5e_8_up, -c(GWAS), values_to = "Bias", names_to = "Method")
new_bias_5e_8_down <- pivot_longer(bias_5e_8_down, -c(GWAS), values_to = "Bias", names_to = "Method")
new_bias_5e_4_up <- pivot_longer(bias_5e_4_up, -c(GWAS), values_to = "Bias", names_to = "Method")
new_bias_5e_4_down <- pivot_longer(bias_5e_4_down, -c(GWAS), values_to = "Bias", names_to = "Method")

new_bias_5e_8_up$Method <- factor(new_bias_5e_8_up$Method, levels=c("naive", "CL1", "CL2", "CL3", "EB", "EB_df", "EB_gam_nb", "EB_gam_po", "EB_scam", "boot", "FIQT"))
new_bias_5e_8_down$Method <- factor(new_bias_5e_8_down$Method, levels=c("naive", "CL1", "CL2", "CL3", "EB", "EB_df", "EB_gam_nb", "EB_gam_po", "EB_scam", "boot", "FIQT"))
new_bias_5e_4_up$Method <- factor(new_bias_5e_4_up$Method, levels=c("naive", "CL1", "CL2", "CL3", "EB", "EB_df", "EB_gam_nb", "EB_gam_po", "EB_scam", "boot", "FIQT"))
new_bias_5e_4_down$Method <- factor(new_bias_5e_4_down$Method, levels=c("naive", "CL1", "CL2", "CL3", "EB", "EB_df", "EB_gam_nb", "EB_gam_po", "EB_scam", "boot", "FIQT"))
new_bias_5e_8_up$GWAS <- factor(new_bias_5e_8_up$GWAS, levels=c("BMI 1", "BMI 2", "T2D 1", "T2D 2", "Height 1", "Height 2"))
new_bias_5e_8_down$GWAS <- factor(new_bias_5e_8_down$GWAS, levels=c("BMI 1", "BMI 2", "T2D 1", "T2D 2", "Height 1", "Height 2"))
new_bias_5e_4_up$GWAS <- factor(new_bias_5e_4_up$GWAS, levels=c("BMI 1", "BMI 2", "T2D 1", "T2D 2", "Height 1", "Height 2"))
new_bias_5e_4_down$GWAS <- factor(new_bias_5e_4_down$GWAS, levels=c("BMI 1", "BMI 2", "T2D 1", "T2D 2", "Height 1", "Height 2"))



plot_1 <- ggplot(new_bias_5e_8_up,aes(x=Bias,y=Method,fill=Method)) + geom_bar(stat="identity")+ geom_col(width=1) + facet_wrap(.~GWAS,nrow=6, strip.position = "right") + scale_shape_manual(values=c(6,20,20,20,15,18,18,18,18,10,15)) + scale_fill_manual(values=c(col[1],col[2],col[8],col[4],col[5],col[3],col[6],col[7],col1[1],col1[11],col1[9])) + ylab("Method") + 
  xlab(expression(paste(italic("Average"), " Bias at ", 5%*%10^-8, " for positive SNPs "))) + theme_bw() + theme(text = element_text(size=12),legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text.y = element_text(angle = 0)) + geom_vline(data = new_bias_5e_8_up[new_bias_5e_8_up$Method=="naive",], aes(xintercept = Bias), linetype=2, size=0.5) + geom_vline(xintercept=0)

plot_2 <- ggplot(new_bias_5e_8_down,aes(x=Bias,y=Method,fill=Method)) + geom_bar(stat="identity")+ geom_col(width=1) + facet_wrap(.~GWAS,nrow=6, strip.position = "right") + scale_shape_manual(values=c(6,20,20,20,15,18,18,18,18,10,15)) + scale_fill_manual(values=c(col[1],col[2],col[8],col[4],col[5],col[3],col[6],col[7],col1[1],col1[11],col1[9])) + ylab("Method") + 
  xlab(expression(paste(italic("Average"), " Bias at ", 5%*%10^-8, " for negative SNPs "))) + theme_bw() + theme(text = element_text(size=12),legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text.y = element_text(angle = 0)) + geom_vline(data = new_bias_5e_8_down[new_bias_5e_8_down$Method=="naive",], aes(xintercept = Bias), linetype=2, size=0.5) + geom_vline(xintercept=0)

plot_3 <- ggplot(new_bias_5e_4_up,aes(x=Bias,y=Method,fill=Method)) + geom_bar(stat="identity")+ geom_col(width=1) + facet_wrap(GWAS~.,nrow=6, strip.position = "right") + scale_shape_manual(values=c(6,20,20,20,15,18,18,18,18,10,15)) + scale_fill_manual(values=c(col[1],col[2],col[8],col[4],col[5],col[3],col[6],col[7],col1[1],col1[11],col1[9])) + ylab("Method") + 
  xlab(expression(paste(italic("Average"), " Bias at ", 5%*%10^-4, " for positive SNPs "))) + theme_bw() + theme(text = element_text(size=12),legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text.y = element_text(angle=0)) + geom_vline(data = new_bias_5e_4_up[new_bias_5e_4_up$Method=="naive",], aes(xintercept = Bias), linetype=2, size=0.5) + geom_vline(xintercept=0)

plot_4 <- ggplot(new_bias_5e_4_down,aes(x=Bias,y=Method,fill=Method)) + geom_bar(stat="identity")+ geom_col(width=1) + facet_wrap(GWAS~.,nrow=6, strip.position = "right") + scale_shape_manual(values=c(6,20,20,20,15,18,18,18,18,10,15)) + scale_fill_manual(values=c(col[1],col[2],col[8],col[4],col[5],col[3],col[6],col[7],col1[1],col1[11],col1[9])) + ylab("Method") + 
  xlab(expression(paste(italic("Average"), " Bias at ", 5%*%10^-4, " for negative SNPs "))) + theme_bw() + theme(text = element_text(size=12),legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text.y = element_text(angle=0)) + geom_vline(data = new_bias_5e_4_down[new_bias_5e_4_down$Method=="naive",], aes(xintercept = Bias), linetype=2, size=0.5) + geom_vline(xintercept=0)


plot_1
plot_2
plot_3
plot_4

###############################################################################
###############################################################################

## COMPUTING ESTIMATED MSE FOR BMI, METHOD AT LOWER THRESHOLDS, 5e-10, 5e-12, 5e-14, 

###############################################################################

## BMI 1

# Apply methods
summary_data_bmi_1 <-  read.table('data/summary_data_bmi_1.txt',header=TRUE)
summary_stats <- summary_data_bmi_1[,3:5]
out_CL <- conditional_likelihood(summary_data = summary_stats, alpha=5e-10)
out_EB <- empirical_bayes(summary_data = summary_stats)
out_EB_df <- empirical_bayes(summary_data=summary_stats, method="fix_df")
out_EB_scam <- empirical_bayes(summary_data=summary_stats,method="scam")
out_EB_gam_po <- empirical_bayes(summary_data=summary_stats,method="gam_po")
out_EB_gam_nb <- empirical_bayes(summary_data=summary_stats,method="gam_nb")
out_FIQT <- FDR_IQT(summary_data = summary_stats)
out_BR <- BR_ss(summary_data = summary_stats)

# Use threshold 5e-10
out_EB_sig <- out_EB[2*(pnorm(abs(out_EB$beta/out_EB$se), lower.tail=FALSE)) < 5e-10,]
out_EB_df_sig <- out_EB_df[2*(pnorm(abs(out_EB_df$beta/out_EB_df$se), lower.tail=FALSE)) < 5e-10,]
out_EB_scam_sig <- out_EB_scam[2*(pnorm(abs(out_EB_scam$beta/out_EB_scam$se), lower.tail=FALSE)) < 5e-10,]
out_EB_gampo_sig <- out_EB_gam_po[2*(pnorm(abs(out_EB_gam_po$beta/out_EB_gam_po$se), lower.tail=FALSE)) < 5e-10,]
out_EB_gamnb_sig <- out_EB_gam_nb[2*(pnorm(abs(out_EB_gam_nb$beta/out_EB_gam_nb$se), lower.tail=FALSE)) < 5e-10,]
out_FIQT_sig <- out_FIQT[2*(pnorm(abs(out_FIQT$beta/out_FIQT$se), lower.tail=FALSE)) < 5e-10,]
out_BR_ss_sig <- out_BR[2*(pnorm(abs(out_BR$beta/out_BR$se), lower.tail=FALSE)) < 5e-10,]
summary_data_bmi_1 <- dplyr::arrange(summary_data_bmi_1,dplyr::desc(abs(summary_data_bmi_1$beta/summary_data_bmi_1$se)))
summary_data_sig_bmi_1 <- summary_data_bmi_1[2*(pnorm(abs(summary_data_bmi_1$beta/summary_data_bmi_1$se), lower.tail=FALSE)) < 5e-10,]
true_beta <- summary_data_sig_bmi_1$beta_rep
# Compute estimated MSE for each method for 5e-10 threshold
mse_bmi1_5e_10 <- data.frame(naive = mean((true_beta - out_EB_sig$beta)^2) - mean(out_EB_sig$se^2), CL1 = mean((true_beta - out_CL$beta.cl1)^2) - mean(out_EB_sig$se^2), CL2 = mean((true_beta - out_CL$beta.cl2)^2) - mean(out_EB_sig$se^2), CL3 = mean((true_beta - out_CL$beta.cl3)^2) - mean(out_EB_sig$se^2), EB = mean((true_beta - out_EB_sig$beta_EB)^2) - mean(out_EB_sig$se^2),  EB_df = mean((true_beta - out_EB_df_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_scam = mean((true_beta - out_EB_scam_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_po = mean((true_beta - out_EB_gampo_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_nb = mean((true_beta - out_EB_gamnb_sig$beta_EB)^2) - mean(out_EB_sig$se^2), boot = mean((true_beta - out_BR_ss_sig$beta_BR_ss)^2) - mean(out_EB_sig$se^2), FIQT = mean((true_beta - out_FIQT_sig$beta_FIQT)^2) - mean(out_EB_sig$se^2))

# Use threshold 5e-12 (need to re-apply conditional likelihood with change of threshold)
out_EB_sig <- out_EB[2*(pnorm(abs(out_EB$beta/out_EB$se), lower.tail=FALSE)) < 5e-12,]
out_EB_df_sig <- out_EB_df[2*(pnorm(abs(out_EB_df$beta/out_EB_df$se), lower.tail=FALSE)) < 5e-12,]
out_EB_scam_sig <- out_EB_scam[2*(pnorm(abs(out_EB_scam$beta/out_EB_scam$se), lower.tail=FALSE)) < 5e-12,]
out_EB_gampo_sig <- out_EB_gam_po[2*(pnorm(abs(out_EB_gam_po$beta/out_EB_gam_po$se), lower.tail=FALSE)) < 5e-12,]
out_EB_gamnb_sig <- out_EB_gam_nb[2*(pnorm(abs(out_EB_gam_nb$beta/out_EB_gam_nb$se), lower.tail=FALSE)) < 5e-12,]
out_FIQT_sig <- out_FIQT[2*(pnorm(abs(out_FIQT$beta/out_FIQT$se), lower.tail=FALSE)) < 5e-12,]
out_BR_ss_sig <- out_BR[2*(pnorm(abs(out_BR$beta/out_BR$se), lower.tail=FALSE)) < 5e-12,]
out_CL <- conditional_likelihood(summary_data = summary_stats, alpha=5e-12)
summary_data_sig2 <- summary_data_bmi_1[2*(pnorm(abs(summary_data_bmi_1$beta/summary_data_bmi_1$se), lower.tail=FALSE)) < 5e-12,]
true_beta2 <- summary_data_sig2$beta_rep
# Compute estimated MSE for each method for 5e-12 threshold
mse_bmi1_5e_12 <- data.frame(naive = mean((true_beta2 - out_EB_sig$beta)^2) - mean(out_EB_sig$se^2), CL1 = mean((true_beta2 - out_CL$beta.cl1)^2) - mean(out_EB_sig$se^2), CL2 = mean((true_beta2 - out_CL$beta.cl2)^2) - mean(out_EB_sig$se^2), CL3 = mean((true_beta2 - out_CL$beta.cl3)^2) - mean(out_EB_sig$se^2), EB = mean((true_beta2 - out_EB_sig$beta_EB)^2) - mean(out_EB_sig$se^2),  EB_df = mean((true_beta2 - out_EB_df_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_scam = mean((true_beta2 - out_EB_scam_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_po = mean((true_beta2 - out_EB_gampo_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_nb = mean((true_beta2 - out_EB_gamnb_sig$beta_EB)^2) - mean(out_EB_sig$se^2), boot = mean((true_beta2 - out_BR_ss_sig$beta_BR_ss)^2) - mean(out_EB_sig$se^2), FIQT = mean((true_beta2 - out_FIQT_sig$beta_FIQT)^2) - mean(out_EB_sig$se^2))

# Use threshold 5e-14 (need to re-apply conditional likelihood with change of threshold)
out_EB_sig <- out_EB[2*(pnorm(abs(out_EB$beta/out_EB$se), lower.tail=FALSE)) < 5e-14,]
out_EB_df_sig <- out_EB_df[2*(pnorm(abs(out_EB_df$beta/out_EB_df$se), lower.tail=FALSE)) < 5e-14,]
out_EB_scam_sig <- out_EB_scam[2*(pnorm(abs(out_EB_scam$beta/out_EB_scam$se), lower.tail=FALSE)) < 5e-14,]
out_EB_gampo_sig <- out_EB_gam_po[2*(pnorm(abs(out_EB_gam_po$beta/out_EB_gam_po$se), lower.tail=FALSE)) < 5e-14,]
out_EB_gamnb_sig <- out_EB_gam_nb[2*(pnorm(abs(out_EB_gam_nb$beta/out_EB_gam_nb$se), lower.tail=FALSE)) < 5e-14,]
out_FIQT_sig <- out_FIQT[2*(pnorm(abs(out_FIQT$beta/out_FIQT$se), lower.tail=FALSE)) < 5e-14,]
out_BR_ss_sig <- out_BR[2*(pnorm(abs(out_BR$beta/out_BR$se), lower.tail=FALSE)) < 5e-14,]
out_CL <- conditional_likelihood(summary_data = summary_stats, alpha=5e-14)
summary_data_sig2 <- summary_data_bmi_1[2*(pnorm(abs(summary_data_bmi_1$beta/summary_data_bmi_1$se), lower.tail=FALSE)) < 5e-14,]
true_beta2 <- summary_data_sig2$beta_rep
# Compute estimated MSE for each method for 5e-14 threshold
mse_bmi1_5e_14 <- data.frame(naive = mean((true_beta2 - out_EB_sig$beta)^2) - mean(out_EB_sig$se^2), CL1 = mean((true_beta2 - out_CL$beta.cl1)^2) - mean(out_EB_sig$se^2), CL2 = mean((true_beta2 - out_CL$beta.cl2)^2) - mean(out_EB_sig$se^2), CL3 = mean((true_beta2 - out_CL$beta.cl3)^2) - mean(out_EB_sig$se^2), EB = mean((true_beta2 - out_EB_sig$beta_EB)^2) - mean(out_EB_sig$se^2),  EB_df = mean((true_beta2 - out_EB_df_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_scam = mean((true_beta2 - out_EB_scam_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_po = mean((true_beta2 - out_EB_gampo_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_nb = mean((true_beta2 - out_EB_gamnb_sig$beta_EB)^2) - mean(out_EB_sig$se^2), boot = mean((true_beta2 - out_BR_ss_sig$beta_BR_ss)^2) - mean(out_EB_sig$se^2), FIQT = mean((true_beta2 - out_FIQT_sig$beta_FIQT)^2) - mean(out_EB_sig$se^2))

###############################################################################

## bmi 2

# Apply methods
summary_data_bmi_2 <-  read.table('data/summary_data_bmi_2.txt',header=TRUE)
summary_stats <- summary_data_bmi_2[,3:5]
out_CL <- conditional_likelihood(summary_data = summary_stats, alpha=5e-10)
out_EB <- empirical_bayes(summary_data = summary_stats)
out_EB_df <- empirical_bayes(summary_data=summary_stats, method="fix_df")
out_EB_scam <- empirical_bayes(summary_data=summary_stats,method="scam")
out_EB_gam_po <- empirical_bayes(summary_data=summary_stats,method="gam_po")
out_EB_gam_nb <- empirical_bayes(summary_data=summary_stats,method="gam_nb")
out_FIQT <- FDR_IQT(summary_data = summary_stats)
out_BR <- BR_ss(summary_data = summary_stats)

# Use threshold 5e-10
out_EB_sig <- out_EB[2*(pnorm(abs(out_EB$beta/out_EB$se), lower.tail=FALSE)) < 5e-10,]
out_EB_df_sig <- out_EB_df[2*(pnorm(abs(out_EB_df$beta/out_EB_df$se), lower.tail=FALSE)) < 5e-10,]
out_EB_scam_sig <- out_EB_scam[2*(pnorm(abs(out_EB_scam$beta/out_EB_scam$se), lower.tail=FALSE)) < 5e-10,]
out_EB_gampo_sig <- out_EB_gam_po[2*(pnorm(abs(out_EB_gam_po$beta/out_EB_gam_po$se), lower.tail=FALSE)) < 5e-10,]
out_EB_gamnb_sig <- out_EB_gam_nb[2*(pnorm(abs(out_EB_gam_nb$beta/out_EB_gam_nb$se), lower.tail=FALSE)) < 5e-10,]
out_FIQT_sig <- out_FIQT[2*(pnorm(abs(out_FIQT$beta/out_FIQT$se), lower.tail=FALSE)) < 5e-10,]
out_BR_ss_sig <- out_BR[2*(pnorm(abs(out_BR$beta/out_BR$se), lower.tail=FALSE)) < 5e-10,]
summary_data_bmi_2 <- dplyr::arrange(summary_data_bmi_2,dplyr::desc(abs(summary_data_bmi_2$beta/summary_data_bmi_2$se)))
summary_data_sig_bmi_2 <- summary_data_bmi_2[2*(pnorm(abs(summary_data_bmi_2$beta/summary_data_bmi_2$se), lower.tail=FALSE)) < 5e-10,]
true_beta <- summary_data_sig_bmi_2$beta_rep
# Compute estimated MSE for each method for 5e-10 threshold
mse_bmi2_5e_10 <- data.frame(naive = mean((true_beta - out_EB_sig$beta)^2) - mean(out_EB_sig$se^2), CL1 = mean((true_beta - out_CL$beta.cl1)^2) - mean(out_EB_sig$se^2), CL2 = mean((true_beta - out_CL$beta.cl2)^2) - mean(out_EB_sig$se^2), CL3 = mean((true_beta - out_CL$beta.cl3)^2) - mean(out_EB_sig$se^2), EB = mean((true_beta - out_EB_sig$beta_EB)^2) - mean(out_EB_sig$se^2),  EB_df = mean((true_beta - out_EB_df_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_scam = mean((true_beta - out_EB_scam_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_po = mean((true_beta - out_EB_gampo_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_nb = mean((true_beta - out_EB_gamnb_sig$beta_EB)^2) - mean(out_EB_sig$se^2), boot = mean((true_beta - out_BR_ss_sig$beta_BR_ss)^2) - mean(out_EB_sig$se^2), FIQT = mean((true_beta - out_FIQT_sig$beta_FIQT)^2) - mean(out_EB_sig$se^2))

# Use threshold 5e-12 (need to re-apply conditional likelihood with change of threshold)
out_EB_sig <- out_EB[2*(pnorm(abs(out_EB$beta/out_EB$se), lower.tail=FALSE)) < 5e-12,]
out_EB_df_sig <- out_EB_df[2*(pnorm(abs(out_EB_df$beta/out_EB_df$se), lower.tail=FALSE)) < 5e-12,]
out_EB_scam_sig <- out_EB_scam[2*(pnorm(abs(out_EB_scam$beta/out_EB_scam$se), lower.tail=FALSE)) < 5e-12,]
out_EB_gampo_sig <- out_EB_gam_po[2*(pnorm(abs(out_EB_gam_po$beta/out_EB_gam_po$se), lower.tail=FALSE)) < 5e-12,]
out_EB_gamnb_sig <- out_EB_gam_nb[2*(pnorm(abs(out_EB_gam_nb$beta/out_EB_gam_nb$se), lower.tail=FALSE)) < 5e-12,]
out_FIQT_sig <- out_FIQT[2*(pnorm(abs(out_FIQT$beta/out_FIQT$se), lower.tail=FALSE)) < 5e-12,]
out_BR_ss_sig <- out_BR[2*(pnorm(abs(out_BR$beta/out_BR$se), lower.tail=FALSE)) < 5e-12,]
out_CL <- conditional_likelihood(summary_data = summary_stats, alpha=5e-12)
summary_data_sig2 <- summary_data_bmi_2[2*(pnorm(abs(summary_data_bmi_2$beta/summary_data_bmi_2$se), lower.tail=FALSE)) < 5e-12,]
true_beta2 <- summary_data_sig2$beta_rep
# Compute estimated MSE for each method for 5e-12 threshold
mse_bmi2_5e_12 <- data.frame(naive = mean((true_beta2 - out_EB_sig$beta)^2) - mean(out_EB_sig$se^2), CL1 = mean((true_beta2 - out_CL$beta.cl1)^2) - mean(out_EB_sig$se^2), CL2 = mean((true_beta2 - out_CL$beta.cl2)^2) - mean(out_EB_sig$se^2), CL3 = mean((true_beta2 - out_CL$beta.cl3)^2) - mean(out_EB_sig$se^2), EB = mean((true_beta2 - out_EB_sig$beta_EB)^2) - mean(out_EB_sig$se^2),  EB_df = mean((true_beta2 - out_EB_df_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_scam = mean((true_beta2 - out_EB_scam_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_po = mean((true_beta2 - out_EB_gampo_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_nb = mean((true_beta2 - out_EB_gamnb_sig$beta_EB)^2) - mean(out_EB_sig$se^2), boot = mean((true_beta2 - out_BR_ss_sig$beta_BR_ss)^2) - mean(out_EB_sig$se^2), FIQT = mean((true_beta2 - out_FIQT_sig$beta_FIQT)^2) - mean(out_EB_sig$se^2))

# Use threshold 5e-14 (need to re-apply conditional likelihood with change of threshold)
out_EB_sig <- out_EB[2*(pnorm(abs(out_EB$beta/out_EB$se), lower.tail=FALSE)) < 5e-14,]
out_EB_df_sig <- out_EB_df[2*(pnorm(abs(out_EB_df$beta/out_EB_df$se), lower.tail=FALSE)) < 5e-14,]
out_EB_scam_sig <- out_EB_scam[2*(pnorm(abs(out_EB_scam$beta/out_EB_scam$se), lower.tail=FALSE)) < 5e-14,]
out_EB_gampo_sig <- out_EB_gam_po[2*(pnorm(abs(out_EB_gam_po$beta/out_EB_gam_po$se), lower.tail=FALSE)) < 5e-14,]
out_EB_gamnb_sig <- out_EB_gam_nb[2*(pnorm(abs(out_EB_gam_nb$beta/out_EB_gam_nb$se), lower.tail=FALSE)) < 5e-14,]
out_FIQT_sig <- out_FIQT[2*(pnorm(abs(out_FIQT$beta/out_FIQT$se), lower.tail=FALSE)) < 5e-14,]
out_BR_ss_sig <- out_BR[2*(pnorm(abs(out_BR$beta/out_BR$se), lower.tail=FALSE)) < 5e-14,]
out_CL <- conditional_likelihood(summary_data = summary_stats, alpha=5e-14)
summary_data_sig2 <- summary_data_bmi_2[2*(pnorm(abs(summary_data_bmi_2$beta/summary_data_bmi_2$se), lower.tail=FALSE)) < 5e-14,]
true_beta2 <- summary_data_sig2$beta_rep
# Compute estimated MSE for each method for 5e-14 threshold
mse_bmi2_5e_14 <- data.frame(naive = mean((true_beta2 - out_EB_sig$beta)^2) - mean(out_EB_sig$se^2), CL1 = mean((true_beta2 - out_CL$beta.cl1)^2) - mean(out_EB_sig$se^2), CL2 = mean((true_beta2 - out_CL$beta.cl2)^2) - mean(out_EB_sig$se^2), CL3 = mean((true_beta2 - out_CL$beta.cl3)^2) - mean(out_EB_sig$se^2), EB = mean((true_beta2 - out_EB_sig$beta_EB)^2) - mean(out_EB_sig$se^2),  EB_df = mean((true_beta2 - out_EB_df_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_scam = mean((true_beta2 - out_EB_scam_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_po = mean((true_beta2 - out_EB_gampo_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_nb = mean((true_beta2 - out_EB_gamnb_sig$beta_EB)^2) - mean(out_EB_sig$se^2), boot = mean((true_beta2 - out_BR_ss_sig$beta_BR_ss)^2) - mean(out_EB_sig$se^2), FIQT = mean((true_beta2 - out_FIQT_sig$beta_FIQT)^2) - mean(out_EB_sig$se^2))

###############################################################################

## Combine all above results and save them in external results files 
mse_5e_10 <- rbind(mse_bmi1_5e_10, mse_bmi2_5e_10)
mse_5e_10 <- round(mse_5e_10[,1:11],5)
GWAS <- c("BMI 1", "BMI 2")
mse_5e_10 <- cbind(GWAS,mse_5e_10)
write.csv(mse_5e_10,"results/mse_5e_10_BMI.txt", row.names = FALSE)

mse_5e_12 <- rbind(mse_bmi1_5e_12, mse_bmi2_5e_12)
mse_5e_12 <- round(mse_5e_12[,1:11],5)
GWAS <- c("BMI 1", "BMI 2")
mse_5e_12 <- cbind(GWAS,mse_5e_12)
write.csv(mse_5e_12,"results/mse_5e_12_BMI.txt", row.names = FALSE)

mse_5e_14 <- rbind(mse_bmi1_5e_14, mse_bmi2_5e_14)
mse_5e_14 <- round(mse_5e_14[,1:11],5)
GWAS <- c("BMI 1", "BMI 2")
mse_5e_14 <- cbind(GWAS,mse_5e_14)
write.csv(mse_5e_14,"results/mse_5e_14_BMI.txt", row.names = FALSE)

###############################################################################

## PART 2) PLOT RESULTS APPROPRIATELY

new_mse_5e_10 <- pivot_longer(mse_5e_10, -c(GWAS), values_to = "MSE", names_to = "Method")
new_mse_5e_12 <- pivot_longer(mse_5e_12, -c(GWAS), values_to = "MSE", names_to = "Method")
new_mse_5e_14 <- pivot_longer(mse_5e_14, -c(GWAS), values_to = "MSE", names_to = "Method")

new_mse_5e_10$Method <- factor(new_mse_5e_10$Method, levels=c("naive", "CL1", "CL2", "CL3", "EB", "EB_df", "EB_gam_nb", "EB_gam_po", "EB_scam", "boot", "FIQT"))
new_mse_5e_10$GWAS <- factor(new_mse_5e_10$GWAS, levels=c("BMI 1", "BMI 2"))

new_mse_5e_12$Method <- factor(new_mse_5e_12$Method, levels=c("naive", "CL1", "CL2", "CL3", "EB", "EB_df", "EB_gam_nb", "EB_gam_po", "EB_scam", "boot", "FIQT"))
new_mse_5e_12$GWAS <- factor(new_mse_5e_12$GWAS, levels=c("BMI 1", "BMI 2"))

new_mse_5e_14$Method <- factor(new_mse_5e_14$Method, levels=c("naive", "CL1", "CL2", "CL3", "EB", "EB_df", "EB_gam_nb", "EB_gam_po", "EB_scam", "boot", "FIQT"))
new_mse_5e_14$GWAS <- factor(new_mse_5e_14$GWAS, levels=c("BMI 1", "BMI 2"))


sum(2*(pnorm(abs(summary_data_bmi_1$beta/summary_data_bmi_1$se), lower.tail=FALSE)) < 5e-10)
sum(2*(pnorm(abs(summary_data_bmi_1$beta/summary_data_bmi_1$se), lower.tail=FALSE)) < 5e-12)
sum(2*(pnorm(abs(summary_data_bmi_1$beta/summary_data_bmi_1$se), lower.tail=FALSE)) < 5e-14)

sum(2*(pnorm(abs(summary_data_bmi_2$beta/summary_data_bmi_2$se), lower.tail=FALSE)) < 5e-10)
sum(2*(pnorm(abs(summary_data_bmi_2$beta/summary_data_bmi_2$se), lower.tail=FALSE)) < 5e-12)
sum(2*(pnorm(abs(summary_data_bmi_2$beta/summary_data_bmi_2$se), lower.tail=FALSE)) < 5e-14)


ann_text1 <- data.frame(MSE = 0.00145,Method = "boot",lab = "label", GWAS = factor("BMI 1",levels = c("BMI 1", "BMI 2")))
ann_text2 <- data.frame(MSE = 0.00145,Method = "boot",lab = "label", GWAS = factor("BMI 2",levels = c("BMI 1", "BMI 2")))
plot_10 <- ggplot(new_mse_5e_10,aes(x=MSE,y=Method,fill=Method)) + geom_bar(stat="identity")+ geom_col(width=1) + facet_wrap(.~GWAS,nrow=6, strip.position = "right") + scale_shape_manual(values=c(6,20,20,20,15,18,18,18,18,10,15)) + scale_fill_manual(values=c(col[1],col[2],col[8],col[4],col[5],col[3],col[6],col[7],col1[1],col1[11],col1[9])) + ylab("Method") + 
  xlab(expression(paste(italic("Estimated"), " MSE at ", 5%*%10^-10))) + theme_bw() + theme(text = element_text(size=12),legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text.y = element_text(angle = 0)) + geom_vline(data = new_mse_5e_10[new_mse_5e_10$Method=="naive",], aes(xintercept = MSE), linetype=2, size=0.5) + geom_vline(xintercept=0) + 
  geom_label(data = ann_text1,label = "3333 sig SNPs",size=3.5,fill="#F4F4F4") +
  geom_label(data = ann_text2,label = "4107 sig SNPs",size=3.5,fill="#F4F4F4") 


ann_text1 <- data.frame(MSE = 0.0019,Method = "boot",lab = "label", GWAS = factor("BMI 1",levels = c("BMI 1", "BMI 2")))
ann_text2 <- data.frame(MSE = 0.0019,Method = "boot",lab = "label", GWAS = factor("BMI 2",levels = c("BMI 1", "BMI 2")))
plot_12 <- ggplot(new_mse_5e_12,aes(x=MSE,y=Method,fill=Method)) + geom_bar(stat="identity")+ geom_col(width=1) + facet_wrap(GWAS~.,nrow=6, strip.position = "right") + scale_shape_manual(values=c(6,20,20,20,15,18,18,18,18,10,15)) + scale_fill_manual(values=c(col[1],col[2],col[8],col[4],col[5],col[3],col[6],col[7],col1[1],col1[11],col1[9])) + ylab("Method") + 
  xlab(expression(paste(italic("Estimated"), " MSE at ", 5%*%10^-12))) + theme_bw() + theme(text = element_text(size=12),legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text.y = element_text(angle=0)) + geom_vline(data = new_mse_5e_12[new_mse_5e_12$Method=="naive",], aes(xintercept = MSE), linetype=2, size=0.5) + geom_vline(xintercept=0) + 
  geom_label(data = ann_text1,label = "2113 sig SNPs",size=3.5,fill="#F4F4F4") +
  geom_label(data = ann_text2,label = "2745 sig SNPs",size=3.5,fill="#F4F4F4") 

ann_text1 <- data.frame(MSE = 0.002,Method = "boot",lab = "label", GWAS = factor("BMI 1",levels = c("BMI 1", "BMI 2")))
ann_text2 <- data.frame(MSE = 0.002,Method = "boot",lab = "label", GWAS = factor("BMI 2",levels = c("BMI 1", "BMI 2")))
plot_14 <- ggplot(new_mse_5e_14,aes(x=MSE,y=Method,fill=Method)) + geom_bar(stat="identity")+ geom_col(width=1) + facet_wrap(GWAS~.,nrow=6, strip.position = "right") + scale_shape_manual(values=c(6,20,20,20,15,18,18,18,18,10,15)) + scale_fill_manual(values=c(col[1],col[2],col[8],col[4],col[5],col[3],col[6],col[7],col1[1],col1[11],col1[9])) + ylab("Method") + 
  xlab(expression(paste(italic("Estimated"), " MSE at ", 5%*%10^-14))) + theme_bw() + theme(text = element_text(size=12),legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text.y = element_text(angle=0)) + geom_vline(data = new_mse_5e_14[new_mse_5e_14$Method=="naive",], aes(xintercept = MSE), linetype=2, size=0.5) + geom_vline(xintercept=0) + 
  geom_label(data = ann_text1,label = "1519 sig SNPs",size=3.5,fill="#F4F4F4") +
  geom_label(data = ann_text2,label = "1742 sig SNPs",size=3.5,fill="#F4F4F4") 


library(patchwork)
library(ggpubr)
figure <- plot_10 + plot_12 + plot_14
figure + plot_layout(guides = "collect", ncol=1) + plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = "bold"))

##################################################################################
##################################################################################

## Height: thresholds consider: 5e-32,5e-34,5e-36

###############################################################################

## Height 1

# Apply methods
summary_data_height_1 <-  read.table('data/summary_data_height_1.txt',header=TRUE)
summary_stats <- summary_data_height_1[,3:5]
out_CL <- conditional_likelihood(summary_data = summary_stats, alpha=5e-32)
out_EB <- empirical_bayes(summary_data = summary_stats)
out_EB_df <- empirical_bayes(summary_data=summary_stats, method="fix_df")
out_EB_scam <- empirical_bayes(summary_data=summary_stats,method="scam")
out_EB_gam_po <- empirical_bayes(summary_data=summary_stats,method="gam_po")
out_EB_gam_nb <- empirical_bayes(summary_data=summary_stats,method="gam_nb")
out_FIQT <- FDR_IQT(summary_data = summary_stats)
out_BR <- BR_ss(summary_data = summary_stats)

# Use threshold 5e-32
out_EB_sig <- out_EB[2*(pnorm(abs(out_EB$beta/out_EB$se), lower.tail=FALSE)) < 5e-32,]
out_EB_df_sig <- out_EB_df[2*(pnorm(abs(out_EB_df$beta/out_EB_df$se), lower.tail=FALSE)) < 5e-32,]
out_EB_scam_sig <- out_EB_scam[2*(pnorm(abs(out_EB_scam$beta/out_EB_scam$se), lower.tail=FALSE)) < 5e-32,]
out_EB_gampo_sig <- out_EB_gam_po[2*(pnorm(abs(out_EB_gam_po$beta/out_EB_gam_po$se), lower.tail=FALSE)) < 5e-32,]
out_EB_gamnb_sig <- out_EB_gam_nb[2*(pnorm(abs(out_EB_gam_nb$beta/out_EB_gam_nb$se), lower.tail=FALSE)) < 5e-32,]
out_FIQT_sig <- out_FIQT[2*(pnorm(abs(out_FIQT$beta/out_FIQT$se), lower.tail=FALSE)) < 5e-32,]
out_BR_ss_sig <- out_BR[2*(pnorm(abs(out_BR$beta/out_BR$se), lower.tail=FALSE)) < 5e-32,]
summary_data_height_1 <- dplyr::arrange(summary_data_height_1,dplyr::desc(abs(summary_data_height_1$beta/summary_data_height_1$se)))
summary_data_sig_height_1 <- summary_data_height_1[2*(pnorm(abs(summary_data_height_1$beta/summary_data_height_1$se), lower.tail=FALSE)) < 5e-32,]
true_beta <- summary_data_sig_height_1$beta_rep
# Compute estimated MSE for each method for 5e-32 threshold
mse_height1_5e_32 <- data.frame(naive = mean((true_beta - out_EB_sig$beta)^2) - mean(out_EB_sig$se^2), CL1 = mean((true_beta - out_CL$beta.cl1)^2) - mean(out_EB_sig$se^2), CL2 = mean((true_beta - out_CL$beta.cl2)^2) - mean(out_EB_sig$se^2), CL3 = mean((true_beta - out_CL$beta.cl3)^2) - mean(out_EB_sig$se^2), EB = mean((true_beta - out_EB_sig$beta_EB)^2) - mean(out_EB_sig$se^2),  EB_df = mean((true_beta - out_EB_df_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_scam = mean((true_beta - out_EB_scam_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_po = mean((true_beta - out_EB_gampo_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_nb = mean((true_beta - out_EB_gamnb_sig$beta_EB)^2) - mean(out_EB_sig$se^2), boot = mean((true_beta - out_BR_ss_sig$beta_BR_ss)^2) - mean(out_EB_sig$se^2), FIQT = mean((true_beta - out_FIQT_sig$beta_FIQT)^2) - mean(out_EB_sig$se^2))

# Use threshold 5e-12 (need to re-apply conditional likelihood with change of threshold)
out_EB_sig <- out_EB[2*(pnorm(abs(out_EB$beta/out_EB$se), lower.tail=FALSE)) < 5e-12,]
out_EB_df_sig <- out_EB_df[2*(pnorm(abs(out_EB_df$beta/out_EB_df$se), lower.tail=FALSE)) < 5e-12,]
out_EB_scam_sig <- out_EB_scam[2*(pnorm(abs(out_EB_scam$beta/out_EB_scam$se), lower.tail=FALSE)) < 5e-12,]
out_EB_gampo_sig <- out_EB_gam_po[2*(pnorm(abs(out_EB_gam_po$beta/out_EB_gam_po$se), lower.tail=FALSE)) < 5e-12,]
out_EB_gamnb_sig <- out_EB_gam_nb[2*(pnorm(abs(out_EB_gam_nb$beta/out_EB_gam_nb$se), lower.tail=FALSE)) < 5e-12,]
out_FIQT_sig <- out_FIQT[2*(pnorm(abs(out_FIQT$beta/out_FIQT$se), lower.tail=FALSE)) < 5e-12,]
out_BR_ss_sig <- out_BR[2*(pnorm(abs(out_BR$beta/out_BR$se), lower.tail=FALSE)) < 5e-12,]
out_CL <- conditional_likelihood(summary_data = summary_stats, alpha=5e-12)
summary_data_sig2 <- summary_data_height_1[2*(pnorm(abs(summary_data_height_1$beta/summary_data_height_1$se), lower.tail=FALSE)) < 5e-12,]
true_beta2 <- summary_data_sig2$beta_rep
# Compute estimated MSE for each method for 5e-12 threshold
mse_height1_5e_34 <- data.frame(naive = mean((true_beta2 - out_EB_sig$beta)^2) - mean(out_EB_sig$se^2), CL1 = mean((true_beta2 - out_CL$beta.cl1)^2) - mean(out_EB_sig$se^2), CL2 = mean((true_beta2 - out_CL$beta.cl2)^2) - mean(out_EB_sig$se^2), CL3 = mean((true_beta2 - out_CL$beta.cl3)^2) - mean(out_EB_sig$se^2), EB = mean((true_beta2 - out_EB_sig$beta_EB)^2) - mean(out_EB_sig$se^2),  EB_df = mean((true_beta2 - out_EB_df_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_scam = mean((true_beta2 - out_EB_scam_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_po = mean((true_beta2 - out_EB_gampo_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_nb = mean((true_beta2 - out_EB_gamnb_sig$beta_EB)^2) - mean(out_EB_sig$se^2), boot = mean((true_beta2 - out_BR_ss_sig$beta_BR_ss)^2) - mean(out_EB_sig$se^2), FIQT = mean((true_beta2 - out_FIQT_sig$beta_FIQT)^2) - mean(out_EB_sig$se^2))

# Use threshold 5e-36 (need to re-apply conditional likelihood with change of threshold)
out_EB_sig <- out_EB[2*(pnorm(abs(out_EB$beta/out_EB$se), lower.tail=FALSE)) < 5e-36,]
out_EB_df_sig <- out_EB_df[2*(pnorm(abs(out_EB_df$beta/out_EB_df$se), lower.tail=FALSE)) < 5e-36,]
out_EB_scam_sig <- out_EB_scam[2*(pnorm(abs(out_EB_scam$beta/out_EB_scam$se), lower.tail=FALSE)) < 5e-36,]
out_EB_gampo_sig <- out_EB_gam_po[2*(pnorm(abs(out_EB_gam_po$beta/out_EB_gam_po$se), lower.tail=FALSE)) < 5e-36,]
out_EB_gamnb_sig <- out_EB_gam_nb[2*(pnorm(abs(out_EB_gam_nb$beta/out_EB_gam_nb$se), lower.tail=FALSE)) < 5e-36,]
out_FIQT_sig <- out_FIQT[2*(pnorm(abs(out_FIQT$beta/out_FIQT$se), lower.tail=FALSE)) < 5e-36,]
out_BR_ss_sig <- out_BR[2*(pnorm(abs(out_BR$beta/out_BR$se), lower.tail=FALSE)) < 5e-36,]
out_CL <- conditional_likelihood(summary_data = summary_stats, alpha=5e-36)
summary_data_sig2 <- summary_data_height_1[2*(pnorm(abs(summary_data_height_1$beta/summary_data_height_1$se), lower.tail=FALSE)) < 5e-36,]
true_beta2 <- summary_data_sig2$beta_rep
# Compute estimated MSE for each method for 5e-36 threshold
mse_height1_5e_36 <- data.frame(naive = mean((true_beta2 - out_EB_sig$beta)^2) - mean(out_EB_sig$se^2), CL1 = mean((true_beta2 - out_CL$beta.cl1)^2) - mean(out_EB_sig$se^2), CL2 = mean((true_beta2 - out_CL$beta.cl2)^2) - mean(out_EB_sig$se^2), CL3 = mean((true_beta2 - out_CL$beta.cl3)^2) - mean(out_EB_sig$se^2), EB = mean((true_beta2 - out_EB_sig$beta_EB)^2) - mean(out_EB_sig$se^2),  EB_df = mean((true_beta2 - out_EB_df_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_scam = mean((true_beta2 - out_EB_scam_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_po = mean((true_beta2 - out_EB_gampo_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_nb = mean((true_beta2 - out_EB_gamnb_sig$beta_EB)^2) - mean(out_EB_sig$se^2), boot = mean((true_beta2 - out_BR_ss_sig$beta_BR_ss)^2) - mean(out_EB_sig$se^2), FIQT = mean((true_beta2 - out_FIQT_sig$beta_FIQT)^2) - mean(out_EB_sig$se^2))

###############################################################################

## Height 2

# Apply methods
summary_data_height_2 <-  read.table('data/summary_data_height_2.txt',header=TRUE)
summary_stats <- summary_data_height_2[,3:5]
out_CL <- conditional_likelihood(summary_data = summary_stats, alpha=5e-32)
out_EB <- empirical_bayes(summary_data = summary_stats)
out_EB_df <- empirical_bayes(summary_data=summary_stats, method="fix_df")
out_EB_scam <- empirical_bayes(summary_data=summary_stats,method="scam")
out_EB_gam_po <- empirical_bayes(summary_data=summary_stats,method="gam_po")
out_EB_gam_nb <- empirical_bayes(summary_data=summary_stats,method="gam_nb")
out_FIQT <- FDR_IQT(summary_data = summary_stats)
out_BR <- BR_ss(summary_data = summary_stats)

# Use threshold 5e-32
out_EB_sig <- out_EB[2*(pnorm(abs(out_EB$beta/out_EB$se), lower.tail=FALSE)) < 5e-32,]
out_EB_df_sig <- out_EB_df[2*(pnorm(abs(out_EB_df$beta/out_EB_df$se), lower.tail=FALSE)) < 5e-32,]
out_EB_scam_sig <- out_EB_scam[2*(pnorm(abs(out_EB_scam$beta/out_EB_scam$se), lower.tail=FALSE)) < 5e-32,]
out_EB_gampo_sig <- out_EB_gam_po[2*(pnorm(abs(out_EB_gam_po$beta/out_EB_gam_po$se), lower.tail=FALSE)) < 5e-32,]
out_EB_gamnb_sig <- out_EB_gam_nb[2*(pnorm(abs(out_EB_gam_nb$beta/out_EB_gam_nb$se), lower.tail=FALSE)) < 5e-32,]
out_FIQT_sig <- out_FIQT[2*(pnorm(abs(out_FIQT$beta/out_FIQT$se), lower.tail=FALSE)) < 5e-32,]
out_BR_ss_sig <- out_BR[2*(pnorm(abs(out_BR$beta/out_BR$se), lower.tail=FALSE)) < 5e-32,]
summary_data_height_2 <- dplyr::arrange(summary_data_height_2,dplyr::desc(abs(summary_data_height_2$beta/summary_data_height_2$se)))
summary_data_sig_height_2 <- summary_data_height_2[2*(pnorm(abs(summary_data_height_2$beta/summary_data_height_2$se), lower.tail=FALSE)) < 5e-32,]
true_beta <- summary_data_sig_height_2$beta_rep
# Compute estimated MSE for each method for 5e-32 threshold
mse_height2_5e_32 <- data.frame(naive = mean((true_beta - out_EB_sig$beta)^2) - mean(out_EB_sig$se^2), CL1 = mean((true_beta - out_CL$beta.cl1)^2) - mean(out_EB_sig$se^2), CL2 = mean((true_beta - out_CL$beta.cl2)^2) - mean(out_EB_sig$se^2), CL3 = mean((true_beta - out_CL$beta.cl3)^2) - mean(out_EB_sig$se^2), EB = mean((true_beta - out_EB_sig$beta_EB)^2) - mean(out_EB_sig$se^2),  EB_df = mean((true_beta - out_EB_df_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_scam = mean((true_beta - out_EB_scam_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_po = mean((true_beta - out_EB_gampo_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_nb = mean((true_beta - out_EB_gamnb_sig$beta_EB)^2) - mean(out_EB_sig$se^2), boot = mean((true_beta - out_BR_ss_sig$beta_BR_ss)^2) - mean(out_EB_sig$se^2), FIQT = mean((true_beta - out_FIQT_sig$beta_FIQT)^2) - mean(out_EB_sig$se^2))

# Use threshold 5e-12 (need to re-apply conditional likelihood with change of threshold)
out_EB_sig <- out_EB[2*(pnorm(abs(out_EB$beta/out_EB$se), lower.tail=FALSE)) < 5e-12,]
out_EB_df_sig <- out_EB_df[2*(pnorm(abs(out_EB_df$beta/out_EB_df$se), lower.tail=FALSE)) < 5e-12,]
out_EB_scam_sig <- out_EB_scam[2*(pnorm(abs(out_EB_scam$beta/out_EB_scam$se), lower.tail=FALSE)) < 5e-12,]
out_EB_gampo_sig <- out_EB_gam_po[2*(pnorm(abs(out_EB_gam_po$beta/out_EB_gam_po$se), lower.tail=FALSE)) < 5e-12,]
out_EB_gamnb_sig <- out_EB_gam_nb[2*(pnorm(abs(out_EB_gam_nb$beta/out_EB_gam_nb$se), lower.tail=FALSE)) < 5e-12,]
out_FIQT_sig <- out_FIQT[2*(pnorm(abs(out_FIQT$beta/out_FIQT$se), lower.tail=FALSE)) < 5e-12,]
out_BR_ss_sig <- out_BR[2*(pnorm(abs(out_BR$beta/out_BR$se), lower.tail=FALSE)) < 5e-12,]
out_CL <- conditional_likelihood(summary_data = summary_stats, alpha=5e-12)
summary_data_sig2 <- summary_data_height_2[2*(pnorm(abs(summary_data_height_2$beta/summary_data_height_2$se), lower.tail=FALSE)) < 5e-12,]
true_beta2 <- summary_data_sig2$beta_rep
# Compute estimated MSE for each method for 5e-12 threshold
mse_height2_5e_34 <- data.frame(naive = mean((true_beta2 - out_EB_sig$beta)^2) - mean(out_EB_sig$se^2), CL1 = mean((true_beta2 - out_CL$beta.cl1)^2) - mean(out_EB_sig$se^2), CL2 = mean((true_beta2 - out_CL$beta.cl2)^2) - mean(out_EB_sig$se^2), CL3 = mean((true_beta2 - out_CL$beta.cl3)^2) - mean(out_EB_sig$se^2), EB = mean((true_beta2 - out_EB_sig$beta_EB)^2) - mean(out_EB_sig$se^2),  EB_df = mean((true_beta2 - out_EB_df_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_scam = mean((true_beta2 - out_EB_scam_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_po = mean((true_beta2 - out_EB_gampo_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_nb = mean((true_beta2 - out_EB_gamnb_sig$beta_EB)^2) - mean(out_EB_sig$se^2), boot = mean((true_beta2 - out_BR_ss_sig$beta_BR_ss)^2) - mean(out_EB_sig$se^2), FIQT = mean((true_beta2 - out_FIQT_sig$beta_FIQT)^2) - mean(out_EB_sig$se^2))

# Use threshold 5e-36 (need to re-apply conditional likelihood with change of threshold)
out_EB_sig <- out_EB[2*(pnorm(abs(out_EB$beta/out_EB$se), lower.tail=FALSE)) < 5e-36,]
out_EB_df_sig <- out_EB_df[2*(pnorm(abs(out_EB_df$beta/out_EB_df$se), lower.tail=FALSE)) < 5e-36,]
out_EB_scam_sig <- out_EB_scam[2*(pnorm(abs(out_EB_scam$beta/out_EB_scam$se), lower.tail=FALSE)) < 5e-36,]
out_EB_gampo_sig <- out_EB_gam_po[2*(pnorm(abs(out_EB_gam_po$beta/out_EB_gam_po$se), lower.tail=FALSE)) < 5e-36,]
out_EB_gamnb_sig <- out_EB_gam_nb[2*(pnorm(abs(out_EB_gam_nb$beta/out_EB_gam_nb$se), lower.tail=FALSE)) < 5e-36,]
out_FIQT_sig <- out_FIQT[2*(pnorm(abs(out_FIQT$beta/out_FIQT$se), lower.tail=FALSE)) < 5e-36,]
out_BR_ss_sig <- out_BR[2*(pnorm(abs(out_BR$beta/out_BR$se), lower.tail=FALSE)) < 5e-36,]
out_CL <- conditional_likelihood(summary_data = summary_stats, alpha=5e-36)
summary_data_sig2 <- summary_data_height_2[2*(pnorm(abs(summary_data_height_2$beta/summary_data_height_2$se), lower.tail=FALSE)) < 5e-36,]
true_beta2 <- summary_data_sig2$beta_rep
# Compute estimated MSE for each method for 5e-36 threshold
mse_height2_5e_36 <- data.frame(naive = mean((true_beta2 - out_EB_sig$beta)^2) - mean(out_EB_sig$se^2), CL1 = mean((true_beta2 - out_CL$beta.cl1)^2) - mean(out_EB_sig$se^2), CL2 = mean((true_beta2 - out_CL$beta.cl2)^2) - mean(out_EB_sig$se^2), CL3 = mean((true_beta2 - out_CL$beta.cl3)^2) - mean(out_EB_sig$se^2), EB = mean((true_beta2 - out_EB_sig$beta_EB)^2) - mean(out_EB_sig$se^2),  EB_df = mean((true_beta2 - out_EB_df_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_scam = mean((true_beta2 - out_EB_scam_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_po = mean((true_beta2 - out_EB_gampo_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_nb = mean((true_beta2 - out_EB_gamnb_sig$beta_EB)^2) - mean(out_EB_sig$se^2), boot = mean((true_beta2 - out_BR_ss_sig$beta_BR_ss)^2) - mean(out_EB_sig$se^2), FIQT = mean((true_beta2 - out_FIQT_sig$beta_FIQT)^2) - mean(out_EB_sig$se^2))

###############################################################################

## Combine all above results and save them in external results files 
mse_5e_32 <- rbind(mse_height1_5e_32, mse_height2_5e_32)
mse_5e_32 <- round(mse_5e_32[,1:11],5)
GWAS <- c("Height 1", "Height 2")
mse_5e_32 <- cbind(GWAS,mse_5e_32)
write.csv(mse_5e_32,"results/mse_5e_32_height.txt", row.names = FALSE)

mse_5e_34 <- rbind(mse_height1_5e_34, mse_height2_5e_34)
mse_5e_34 <- round(mse_5e_34[,1:11],5)
GWAS <- c("Height 1", "Height 2")
mse_5e_34 <- cbind(GWAS,mse_5e_34)
write.csv(mse_5e_34,"results/mse_5e_34_height.txt", row.names = FALSE)

mse_5e_36 <- rbind(mse_height1_5e_36, mse_height2_5e_36)
mse_5e_36 <- round(mse_5e_36[,1:11],5)
GWAS <- c("Height 1", "Height 2")
mse_5e_36 <- cbind(GWAS,mse_5e_36)
write.csv(mse_5e_36,"results/mse_5e_36_height.txt", row.names = FALSE)

###############################################################################

## PART 2) PLOT RESULTS APPROPRIATELY

new_mse_5e_32 <- pivot_longer(mse_5e_32, -c(GWAS), values_to = "MSE", names_to = "Method")
new_mse_5e_34 <- pivot_longer(mse_5e_34, -c(GWAS), values_to = "MSE", names_to = "Method")
new_mse_5e_36 <- pivot_longer(mse_5e_36, -c(GWAS), values_to = "MSE", names_to = "Method")

new_mse_5e_32$Method <- factor(new_mse_5e_32$Method, levels=c("naive", "CL1", "CL2", "CL3", "EB", "EB_df", "EB_gam_nb", "EB_gam_po", "EB_scam", "boot", "FIQT"))
new_mse_5e_32$GWAS <- factor(new_mse_5e_32$GWAS, levels=c("Height 1", "Height 2"))

new_mse_5e_34$Method <- factor(new_mse_5e_34$Method, levels=c("naive", "CL1", "CL2", "CL3", "EB", "EB_df", "EB_gam_nb", "EB_gam_po", "EB_scam", "boot", "FIQT"))
new_mse_5e_34$GWAS <- factor(new_mse_5e_34$GWAS, levels=c("Height 1", "Height 2"))

new_mse_5e_36$Method <- factor(new_mse_5e_36$Method, levels=c("naive", "CL1", "CL2", "CL3", "EB", "EB_df", "EB_gam_nb", "EB_gam_po", "EB_scam", "boot", "FIQT"))
new_mse_5e_36$GWAS <- factor(new_mse_5e_36$GWAS, levels=c("Height 1", "Height 2"))


sum(2*(pnorm(abs(summary_data_height_1$beta/summary_data_height_1$se), lower.tail=FALSE)) < 5e-32)
sum(2*(pnorm(abs(summary_data_height_1$beta/summary_data_height_1$se), lower.tail=FALSE)) < 5e-12)
sum(2*(pnorm(abs(summary_data_height_1$beta/summary_data_height_1$se), lower.tail=FALSE)) < 5e-36)

sum(2*(pnorm(abs(summary_data_height_2$beta/summary_data_height_2$se), lower.tail=FALSE)) < 5e-32)
sum(2*(pnorm(abs(summary_data_height_2$beta/summary_data_height_2$se), lower.tail=FALSE)) < 5e-12)
sum(2*(pnorm(abs(summary_data_height_2$beta/summary_data_height_2$se), lower.tail=FALSE)) < 5e-36)


ann_text1 <- data.frame(MSE = 0.0062,Method = "boot",lab = "label", GWAS = factor("Height 1",levels = c("Height 1", "Height 2")))
ann_text2 <- data.frame(MSE = 0.0062,Method = "boot",lab = "label", GWAS = factor("Height 2",levels = c("Height 1", "Height 2")))
plot_32 <- ggplot(new_mse_5e_32,aes(x=MSE,y=Method,fill=Method)) + geom_bar(stat="identity")+ geom_col(width=1) + facet_wrap(.~GWAS,nrow=6, strip.position = "right") + scale_shape_manual(values=c(6,20,20,20,15,18,18,18,18,10,15)) + scale_fill_manual(values=c(col[1],col[2],col[8],col[4],col[5],col[3],col[6],col[7],col1[1],col1[11],col1[9])) + ylab("Method") + 
  xlab(expression(paste(italic("Estimated"), " MSE at ", 5%*%10^-32))) + theme_bw() + theme(text = element_text(size=12),legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text.y = element_text(angle = 0)) + geom_vline(data = new_mse_5e_32[new_mse_5e_32$Method=="naive",], aes(xintercept = MSE), linetype=2, size=0.5) + geom_vline(xintercept=0) + 
  geom_label(data = ann_text1,label = "3459 sig SNPs",size=3.5,fill="#F4F4F4") +
  geom_label(data = ann_text2,label = "4163 sig SNPs",size=3.5,fill="#F4F4F4") 


ann_text1 <- data.frame(MSE = 0.005,Method = "boot",lab = "label", GWAS = factor("Height 1",levels = c("Height 1", "Height 2")))
ann_text2 <- data.frame(MSE = 0.005,Method = "boot",lab = "label", GWAS = factor("Height 2",levels = c("Height 1", "Height 2")))
plot_34 <- ggplot(new_mse_5e_34,aes(x=MSE,y=Method,fill=Method)) + geom_bar(stat="identity")+ geom_col(width=1) + facet_wrap(GWAS~.,nrow=6, strip.position = "right") + scale_shape_manual(values=c(6,20,20,20,15,18,18,18,18,10,15)) + scale_fill_manual(values=c(col[1],col[2],col[8],col[4],col[5],col[3],col[6],col[7],col1[1],col1[11],col1[9])) + ylab("Method") + 
  xlab(expression(paste(italic("Estimated"), " MSE at ", 5%*%10^-34))) + theme_bw() + theme(text = element_text(size=12),legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text.y = element_text(angle=0)) + geom_vline(data = new_mse_5e_34[new_mse_5e_34$Method=="naive",], aes(xintercept = MSE), linetype=2, size=0.5) + geom_vline(xintercept=0) + 
  geom_label(data = ann_text1,label = "3062 sig SNPs",size=3.5,fill="#F4F4F4") +
  geom_label(data = ann_text2,label = "3358 sig SNPs",size=3.5,fill="#F4F4F4") 

ann_text1 <- data.frame(MSE = 0.007,Method = "boot",lab = "label", GWAS = factor("Height 1",levels = c("Height 1", "Height 2")))
ann_text2 <- data.frame(MSE = 0.007,Method = "boot",lab = "label", GWAS = factor("Height 2",levels = c("Height 1", "Height 2")))
plot_36 <- ggplot(new_mse_5e_30,aes(x=MSE,y=Method,fill=Method)) + geom_bar(stat="identity")+ geom_col(width=1) + facet_wrap(GWAS~.,nrow=6, strip.position = "right") + scale_shape_manual(values=c(6,20,20,20,15,18,18,18,18,10,15)) + scale_fill_manual(values=c(col[1],col[2],col[8],col[4],col[5],col[3],col[6],col[7],col1[1],col1[11],col1[9])) + ylab("Method") + 
  xlab(expression(paste(italic("Estimated"), " MSE at ", 5%*%10^-36))) + theme_bw() + theme(text = element_text(size=12),legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text.y = element_text(angle=0)) + geom_vline(data = new_mse_5e_36[new_mse_5e_36$Method=="naive",], aes(xintercept = MSE), linetype=2, size=0.5) + geom_vline(xintercept=0) + 
  geom_label(data = ann_text1,label = "2748 sig SNPs",size=3.5,fill="#F4F4F4") +
  geom_label(data = ann_text2,label = "2867 sig SNPs",size=3.5,fill="#F4F4F4") 

figure <- plot_32 + plot_34 + plot_36
figure + plot_layout(guides = "collect", ncol=1) + plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = "bold"))

###########################################################################################

