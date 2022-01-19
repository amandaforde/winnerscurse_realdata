# Evaluating and comparing Winner's Curse correctiong methods using real data sets

This repository contains various shell and R scripts which have been used to perform suitable GWASs on two different traits, namely BMI and Type 2 diabetes (T2D). These GWAS results allowed us to investigate methods which have been designed to correct Winner's Curse bias in GWAS SNP-trait association estimates. The methods were applied to our obtained sets of summary statistics using the R package, `winnerscurse`. For more details regarding this package and its functionalities, please see [https://amandaforde.github.io/winnerscurse/](https://amandaforde.github.io/winnerscurse/). 

Details of all current progress to date, together with illustrations of results, can be viewed at [https://amandaforde.github.io/winnerscurse_realdata/](https://amandaforde.github.io/winnerscurse_realdata/). Note that our investigation has only considered evaluating methods which use summary statistics from a discovery GWAS, rather than methods which combine information gained from both discovery and replication GWASs.


Our work began from a point in which appropriate quality control procedures had already been performed on UKBB genotype data and we had access to BMI information for 502,554 individuals and T2D status for 502,524 individuals. Description of the scripts used in our workflow is included below. 

&nbsp;

## Generation of GWAS summary statistics for BMI and T2D

### Step 1: Splitting each data set in two  


**split_bmi.R:** 

This R script takes a text file containing individual ID numbers and corresponding *BMI* values for each and randomly splits the file into two separate data sets. The seed of R's random number generator is set to `1998`.  

- ***Input:*** bmi4plink.txt
- ***Output:*** bmi_gwasA.txt, bmi_gwasB.txt


**split_T2D.R:**

This R script takes a text file containing individual ID numbers and corresponding *T2D* status for each and randomly splits the file into two separate data sets. It ensures that T2D status is encoded in 1|2 format to make it suitable for use with PLINK 2.0. The seed of R's random number generator is set to `1998`.

- ***Input:*** T2D4plink.txt
- ***Output:*** T2D_gwasA.txt, T2D_gwasB.txt



&nbsp;


### Step 2: Performing two GWASs for each trait 


**parallelA_bmi.sh:** 

This shell script allows us to run the commands contained in *bmiA.sh* for each chromosome in parallel. 

- ***Input:*** bmiA.sh


**parallelB_bmi.sh:**

This shell script allows us to run the commands contained in *bmiB.sh* for each chromosome in parallel. 

- ***Input:*** bmiB.sh


**parallelA_T2D.sh:**

This shell script allows us to run the commands contained in *T2DA.sh* for each chromosome in parallel.

- ***Input:*** T2DA.sh

**parallelB_T2D.sh:**

This shell script allows us to run the commands contained in *T2DB.sh* for each chromosome in parallel.

- ***Input:*** T2DB.sh


**bmiA.sh:**

This shell script first uses PLINK 2.0 to produce \*_qcd_bmiA.bim, \*_qcd_bmiA.bed, \*_qcd_bmiA.fam files, with the phenotype data being specified in *bmi_gwasA.txt*. A *linear* model is then fitted for every variant in which 8 principal components are included as well as age.

- ***Input:*** \*_qcd.psam, \*_qcd.pvar, \*_qcd.pgen files for each chromosome, bmi_gwasA.txt, covariates.txt
- ***Output:*** bmiA_res_\*.PHENO1.glm.linear for each chromosome

**bmiB.sh:**

This shell script first uses PLINK 2.0 to produce \*_qcd_bmiB.bim, \*_qcd_bmiB.bed, \*_qcd_bmiB.fam files, with the phenotype data being specified in *bmi_gwasB.txt*. A *linear* model is then fitted for every variant in which 8 principal components are included as well as age.

- ***Input:*** \*_qcd.psam, \*_qcd.pvar, \*_qcd.pgen files for each chromosome, bmi_gwasB.txt, covariate.txt
- ***Output:*** bmiB_res_\*.PHENO1.glm.linear for each chromosome

**T2DA.sh:**

This shell script first uses PLINK 2.0 to produce \*_qcd_T2DA.bim, \*_qcd_T2DA.bed, \*_qcd_T2DA.fam files, with the phenotype data being specified in *T2D_gwasA.txt*. As this is a binary trait, a *Firth regression* model is then fitted for every variant in which 8 principal components are included as well as age.


- ***Input:*** \*_qcd.psam, \*_qcd.pvar, \*_qcd.pgen files for each chromosome, T2D_gwasA.txt, covariates.txt
- ***Output:*** T2DA_res_\*.PHENO1.glm.logistic.hybrid for each chromosome

**T2DB.sh:**

This shell script first uses PLINK 2.0 to produce \*_qcd_T2DB.bim, \*_qcd_T2DB.bed, \*_qcd_T2DB.fam files, with the phenotype data being specified in *T2D_gwasB.txt*. As this is a binary trait, a *Firth regression* model is then fitted for every variant in which 8 principal components are included as well as age.

- ***Input:*** \*_qcd.psam, \*_qcd.pvar, \*_qcd.pgen files for each chromosome, T2D_gwasB.txt, covariates.txt
- ***Output:*** T2DB_res_\*.PHENO1.glm.logistic.hybrid for each chromosome



&nbsp;


### Step 3: Combining results 


**summary_data_bmi.R:**

This R script is designed to combine all summary statistics obtained from performing our two *BMI* GWASs together in a suitable format, in which one GWAS is considered to be the discovery GWAS and the other the replication GWAS. The output is a text file which contains information regarding chromosome, position, rsID, effect size and corresponding standard error from the discovery GWAS as well as effect size from the replication GWAS for each variant. 


**summary_data_T2D.R:**

This R script is designed to combine all summary statistics obtained from performing our two *T2D* GWASs together in a suitable format, in which one GWAS is considered to be the discovery GWAS and the other the replication GWAS. The output is a text file which contains information regarding chromosome, position, rsID, effect size (`OR`) and corresponding standard error (`LOG(OR)_SE`) from the discovery GWAS as well as effect size from the replication GWAS for each variant. 



**summary_data_bmi_1.sh:**

This shell script specifies that the files *bmiA*_res_\*.PHENO1.glm.linear for chromosomes 1 to 22 are to be used as the discovery GWAS results and the files *bmiB*_res_\*.PHENO1.glm.linear as the replication GWAS results for the R script summary_data_bmi.R.

- ***Input:*** summary_data_bmi.R, bmiA_res_\*.PHENO1.glm.linear for each chromosome, bmiB_res_\*.PHENO1.glm.linear for each chromosome
- ***Output:*** summary_data_bmi_1.txt

**summary_data_bmi_2.sh:**

This shell script specifies that the files *bmiB*_res_\*.PHENO1.glm.linear for chromosomes 1 to 22 are to be used as the discovery GWAS results and the files *bmiA*_res_\*.PHENO1.glm.linear as the replication GWAS results for the R script summary_data_bmi.R.

- ***Input:*** summary_data_bmi.R, bmiB_res_\*.PHENO1.glm.linear for each chromosome, bmiA_res_\*.PHENO1.glm.linear for each chromosome
- ***Output:*** summary_data_bmi_2.txt


**summary_data_T2D_1.sh:**

This shell script specifies that the files *T2DA*_res_\*.PHENO1.glm.logistic.hybrid for chromosomes 1 to 22 are to be used as the discovery GWAS results and the files *T2DB*_res_\*.PHENO1.glm.logistic.hybrid as the replication GWAS results for the R script summary_data_bmi.R.

- ***Input:*** summary_data_T2D.R, T2DA_res_\*.PHENO1.glm.logistic.hybrid for each chromosome, T2DB_res_\*.PHENO1.glm.logistic.hybrid for each chromosome
- ***Output:*** summary_data_T2D_1.txt


**summary_data_T2D_2.sh:**

This shell script specifies that the files *T2DB*_res_\*.PHENO1.glm.logistic.hybrid for chromosomes 1 to 22 are to be used as the discovery GWAS results and the files *T2DA*_res_\*.PHENO1.glm.logistic.hybrid as the replication GWAS results for the R script summary_data_bmi.R. 


- ***Input:*** summary_data_T2D.R, T2DB_res_\*.PHENO1.glm.logistic.hybrid for each chromosome, T2DA_res_\*.PHENO1.glm.logistic.hybrid for each chromosome
- ***Output:*** summary_data_T2D_2.txt


