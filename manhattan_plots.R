## creating Manhattan plots 

library(qqman)

summary_data_bmi_1 <-  read.table('data/summary_data_bmi_1.txt',header=TRUE)
summary_data_bmi_2 <-  read.table('data/summary_data_bmi_2.txt',header=TRUE)
summary_data_height_1 <-  read.table('data/summary_data_height_1.txt',header=TRUE)
summary_data_height_2 <-  read.table('data/summary_data_height_2.txt',header=TRUE)

summary_data_bmi_1_sub <- summary_data_bmi_1[2*(pnorm(abs(summary_data_bmi_1$beta/summary_data_bmi_1$se), lower.tail=FALSE)) < 0.05,]
summary_data_bmi_1_new <- data.frame(SNP = summary_data_bmi_1_sub$rsid, CHR = summary_data_bmi_1_sub$chr, BP = summary_data_bmi_1_sub$pos, P = 2*(pnorm(abs(summary_data_bmi_1_sub$beta/summary_data_bmi_1_sub$se), lower.tail=FALSE)))
manhattan(summary_data_bmi_1_new, main = "BMI 1", cex = 0.6, cex.axis = 0.9, suggestiveline = -log10(5e-10))
                                                                                      

summary_data_bmi_2_sub <- summary_data_bmi_2[2*(pnorm(abs(summary_data_bmi_2$beta/summary_data_bmi_2$se), lower.tail=FALSE)) < 0.05,]
summary_data_bmi_2_new <- data.frame(SNP = summary_data_bmi_2_sub$rsid, CHR = summary_data_bmi_2_sub$chr, BP = summary_data_bmi_2_sub$pos, P = 2*(pnorm(abs(summary_data_bmi_2_sub$beta/summary_data_bmi_2_sub$se), lower.tail=FALSE)))
manhattan(summary_data_bmi_2_new, main = "BMI 2", cex = 0.6, cex.axis = 0.9, suggestiveline = -log10(5e-12))


summary_data_height_1_sub <- summary_data_height_1[2*(pnorm(abs(summary_data_height_1$beta/summary_data_height_1$se), lower.tail=FALSE)) < 0.05,]
summary_data_height_1_new <- data.frame(SNP = summary_data_height_1_sub$rsid, CHR = summary_data_height_1_sub$chr, BP = summary_data_height_1_sub$pos, P = 2*(pnorm(abs(summary_data_height_1_sub$beta/summary_data_height_1_sub$se), lower.tail=FALSE)))
manhattan(summary_data_height_1_new, main = "Height 1", cex = 0.6, cex.axis = 0.9, suggestiveline = -log10(5e-32))


summary_data_height_2_sub <- summary_data_height_2[2*(pnorm(abs(summary_data_height_2$beta/summary_data_height_2$se), lower.tail=FALSE)) < 0.05,]
summary_data_height_2_new <- data.frame(SNP = summary_data_height_2_sub$rsid, CHR = summary_data_height_2_sub$chr, BP = summary_data_height_2_sub$pos, P = 2*(pnorm(abs(summary_data_height_2_sub$beta/summary_data_height_2_sub$se), lower.tail=FALSE)))
manhattan(summary_data_height_2_new, main = "Height 2", cex = 0.6, cex.axis = 0.9, suggestiveline = -log10(5e-34))


summary_data_T2D_1 <-  read.table('data/summary_data_T2D_1.txt',header=TRUE)
summary_data_T2D_1$beta <- log(summary_data_T2D_1$beta)
summary_data_T2D_1$beta_rep <- log(summary_data_T2D_1$beta_rep)

summary_data_T2D_2 <-  read.table('data/summary_data_T2D_2.txt',header=TRUE)
summary_data_T2D_2$beta <- log(summary_data_T2D_2$beta)
summary_data_T2D_2$beta_rep <- log(summary_data_T2D_2$beta_rep)


summary_data_T2D_1_sub <- summary_data_T2D_1[2*(pnorm(abs(summary_data_T2D_1$beta/summary_data_T2D_1$se), lower.tail=FALSE)) < 0.05,]
summary_data_T2D_1_new <- data.frame(SNP = summary_data_T2D_1_sub$rsid, CHR = summary_data_T2D_1_sub$chr, BP = summary_data_T2D_1_sub$pos, P = 2*(pnorm(abs(summary_data_T2D_1_sub$beta/summary_data_T2D_1_sub$se), lower.tail=FALSE)))
manhattan(summary_data_T2D_1_new, main = "T2D 1", cex = 0.6, cex.axis = 0.9, suggestiveline = -log10(5e-7))


summary_data_T2D_2_sub <- summary_data_T2D_2[2*(pnorm(abs(summary_data_T2D_2$beta/summary_data_T2D_2$se), lower.tail=FALSE)) < 0.05,]
summary_data_T2D_2_new <- data.frame(SNP = summary_data_T2D_2_sub$rsid, CHR = summary_data_T2D_2_sub$chr, BP = summary_data_T2D_2_sub$pos, P = 2*(pnorm(abs(summary_data_T2D_2_sub$beta/summary_data_T2D_2_sub$se), lower.tail=FALSE)))
manhattan(summary_data_T2D_2_new, main = "T2D 2", cex = 0.6, cex.axis = 0.9, suggestiveline = -log10(5e-7))


## save as 1200x600 