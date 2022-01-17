---
title: "PRISM breast milk EV-microRNAs and maternal asthma"
output:
  html_document:
    toc: true
    toc_float: true
---


## Required Packages and functions

```{r,warning=FALSE,message=FALSE,eval=FALSE}
library(reshape)
library(ggplot2)
library(gridExtra)
library(knitr)
library(dplyr)
library(kableExtra)
library(tidyverse)
library(stringr)
library(UpSetR)
library(reshape)
library(corrplot)
library(factoextra)
library(EnvStats)
library(MASS)
library(lmtest)
library(sandwich)
```

## QC and normalization
```{r,warning=FALSE,message=FALSE,eval=FALSE}
# Raw data
miRNAraw = read.csv('PRISM_Breast Milk_EV_miRNA_raw.csv')

# Cotrol probes
cont = c('ath-miR159a_000338', 'RNU44_001094', 'RNU48_001006', 'U6 rRNA_001973')
# ath-miR159a is used as a control
# RNU33 and RNU 48 are artificual controls
# U6 is a conserved snRNA

# Misannotated probes measuring tRNA
tRNA = c('hsa-miR-1274A_002883', 'hsa-miR-1274B_002884')
length(unique(miRNAraw$Target.Name[!(miRNAraw$Target.Name %in% cont) & !(miRNAraw$Target.Name %in% tRNA)]))
# 752 miRNA probes

# Changing Cq to numeric
miRNAraw$Cq[miRNAraw$Cq == 'Undetermined'] = NA
miRNAraw$Cq = as.numeric(as.character(miRNAraw$Cq))

# Dropping Cq values associated with amplification scores < 1.1
table(miRNAraw$Amp.Score < 1.1)
# FALSE  TRUE 
# 22459 38891 

miRNA = miRNAraw[miRNAraw$Amp.Score >= 1.1,]
# 22459    43

# Dropping Cq values associated with Cq confidence values < 0.8
table(miRNA$Cq.Conf < 0.8)
# FALSE  TRUE 
# 18354  4105 

miRNA = miRNA[miRNA$Cq.Conf >= 0.8,]
# 18354    43

# Dropping Cq values less than the control probe U6 (for this dataset, Cq<10.1)
table(miRNA$Cq < 10.1)
# FALSE  TRUE 
# 18323    31 

miRNA = miRNA[miRNA$Cq >= 10.1,]
# 18323    43

# Dropping control probes
miRNA = miRNA[!(miRNA$Target.Name %in% cont),]
# 15887    43

# Dropping misannotated probes measuring tRNA
miRNA = miRNA[!(miRNA$Target.Name %in% tRNA),]
# 15776    43


# Dropping Cq values > 35
table(miRNA$Cq > 35)
# 15767     9 

miRNA = miRNA[miRNA$Cq <= 35,]
# 15767    43

# Correcting IDs
# Correcting ID as per email with Tessa Bloomquist (12/18/19)
miRNA$Sample.Name[miRNA$Sample.Name == 7070201] = 7074201

# Correcting ID as per email with Mathilda Chiu (2/4/20)
miRNA$Sample.Name[miRNA$Sample.Name == 7006102] = 7006103

# Saving long format QCed data
write.csv(miRNA, '/Users/annebozack/Documents/Lee/miRNA asthma/manuscript/clinical epigenetics submission/PRISM_Breast Milk_EV_miRNA_QCed_20200218.csv')
miRNA = read.csv('/Users/annebozack/Documents/Lee/miRNA asthma/manuscript/clinical epigenetics submission/PRISM_Breast Milk_EV_miRNA_QCed_20200218.csv')[,-1]

# Removing special characters in target name (can't be used in column names)
miRNA$Target.Name = as.character(miRNA$Target.Name)
miRNA$Target.Name2 = str_replace_all(miRNA$Target.Name, '[-#]', '.')

# creating a wide format dataset
miRNAwide = data.frame(Sample.Name = miRNA$Sample.Name, Cq = miRNA$Cq, Target.Name2 = miRNA$Target.Name2) %>% spread(Target.Name2, Cq)

# removing columns/probes that are all NA
miRNAwide = miRNAwide[,colSums(is.na(miRNAwide)) < nrow(miRNAwide)]

dim(miRNAwide)[2] - 1
# 550 miRNAs passed QC

# Saving wide format dataset 
write.csv(miRNAwide, 'PRISM breast milk miRNA datasets for Box/PRISM_Breast Milk_EV_miRNA_QCed_wide_20200218.csv')
miRNAwide = read.csv('/Users/annebozack/Documents/Lee/miRNA asthma/PRISM_Breast Milk_EV_miRNA_QCed_wide_20200218.csv')[,-1]

# Number of samples each miRNA detected in and descriptive statistics
df = data.frame(matrix(nrow = length(colnames(miRNAwide[, 2:ncol(miRNAwide)])), ncol = 7))
df[,1] = colnames(miRNAwide[, 2:ncol(miRNAwide)])
colnames(df) = c('miRNA', 'Detected', 'Median', 'IQR', 'Mean', 'SD', '% missing')

for (i in 1:nrow(df)){
	df[i,2] = sum(!is.na(miRNAwide[[df$miRNA[i]]]))
	df[i,3] = median(miRNAwide[[df$miRNA[i]]], na.rm=T)
	df[i,4] = paste0(round(quantile(miRNAwide[[df$miRNA[i]]], na.rm=T), 2)[2], ', ', round(quantile(miRNAwide[[df$miRNA[i]]], na.rm=T), 2)[4])
	df[i,5] = mean(miRNAwide[[df$miRNA[i]]], na.rm=T)
	df[i,6] = sd(miRNAwide[[df$miRNA[i]]], na.rm=T)
	df[i,7] = sum(is.na(miRNAwide[[df$miRNA[i]]]))/75
}

df_nDet = df[order(-df$Detected),]
write.csv(df_nDet, '/Users/annebozack/Documents/Lee/miRNA asthma/manuscript/clinical epigenetics submission/suppTab1.csv')
```

```{r,warning=FALSE,message=FALSE,eval=FALSE}
# Replacing missing Cq values with 35 and normalizing with global mean.
table(is.na(miRNAwide))
# FALSE  TRUE 
 # 15842 25483 

miRNAwide[,c(2:551)][is.na(miRNAwide[,c(2:551)])] = 35

miRNAnorm1 = data.frame(matrix(ncol = ncol(miRNAwide) - 1, nrow = nrow(miRNAwide)))
colnames(miRNAnorm1) = colnames(miRNAwide)[-1]
rownames(miRNAnorm1) = miRNAwide[,1]

# check alignment
all(colnames(miRNAnorm1) == colnames(miRNAwide[-1]))
# TRUE
all(identical(colnames(miRNAnorm1), colnames(miRNAwide[-1])))
# TRUE

for (i in 1:nrow(miRNAnorm1)){
	gm = geoMean(as.numeric(miRNAwide[,-1][i,]))
	for (j in 1:ncol(miRNAnorm1)){
		miRNAnorm1[i,j] = miRNAwide[,-1][i,j] - gm
	}	
}

dim(miRNAnorm1)
#  75 550

# 80% of samples
75*0.8
# 60

table(df_nDet$Detected >= 60)
# FALSE  TRUE 
  # 420   130 

# Keeping only 130 miRNAs detected in >= 80% of samples
miRNAnorm1 = data.frame(Sample.Name = miRNAwide[,1], miRNAnorm1[,colnames(miRNAnorm1) %in% df$miRNA[df$Detected >= 60]])

dim(miRNAnorm1)
#  75 131
```

## Adding pheno data
```{r,warning=FALSE,message=FALSE,eval=FALSE}
pheno = read.csv('/Users/annebozack/Documents/Lee/transfer/PRISM/breastMilkEVmiRNA/maternal_asthma/AB20200922PRISM.csv')
pheno$Sample.Name = NA
pheno$Sample.Name[pheno$iid == 1] = paste0(pheno$hhid[pheno$iid == 1], '01')
pheno$Sample.Name[pheno$iid == 2] = paste0(pheno$hhid[pheno$iid == 2], '02')
pheno$Sample.Name[pheno$iid == 3] = paste0(pheno$hhid[pheno$iid == 3], '03')
pheno$Sample.Name = as.numeric(pheno$Sample.Name)

pheno = pheno[!(is.na(pheno$Sample.Name)),]
pheno = pheno[(pheno$Sample.Name %in% miRNAnorm1$Sample.Name),]

pheno$edu2[pheno$educ <= 2] = 1 # high school or GED
pheno$edu2[pheno$educ > 2] = 2 # > high school
pheno$edu2 = as.factor(pheno$edu2)

pheno$momrace2 = as.factor(pheno$momrace2)
pheno$momrace3[pheno$momrace2 == 0 | pheno$momrace2 == 3] = 0 # White/other
pheno$momrace3[pheno$momrace2 == 1] = 1 # Black
pheno$momrace3[pheno$momrace2 == 2] = 2 # Hispanic
pheno$momrace3 = as.factor(pheno$momrace3)

pheno$prepreg_ob[pheno$mbmi_prepreg < 30] = 0
pheno$prepreg_ob[pheno$mbmi_prepreg >= 30] = 1
pheno$prepreg_ob = as.factor(pheno$prepreg_ob)

# pheno2 = read.csv('/Users/annebozack/Documents/Lee/transfer/PRISM/450K/AB20191202.csv')
# pheno2$Sample.Name = NA
# pheno2$Sample.Name[pheno2$iid == 1 & !is.na(pheno2$iid)] = paste0(as.character(pheno2$hhid[pheno2$iid == 1 & !is.na(pheno2$iid)]), '01')
# pheno2$Sample.Name[pheno2$iid == 2 & !is.na(pheno2$iid)] = paste0(as.character(pheno2$hhid[pheno2$iid == 2 & !is.na(pheno2$iid)]), '02')
# pheno2$Sample.Name[pheno2$iid == 3 & !is.na(pheno2$iid)] = paste0(as.character(pheno2$hhid[pheno2$iid == 3 & !is.na(pheno2$iid)]), '03')
# pheno2$Sample.Name = as.numeric(pheno2$Sample.Name)
# phenoSm = pheno2[,c(31, 8)]

# pheno = merge(pheno, phenoSm, by = 'Sample.Name')

pheno$chidsex = as.factor(pheno$childsex)

pheno$smoke = NA
pheno$smoke[pheno$ets_preg == 1 | pheno$smk_preg == 1] = 1 # smoking or ets >= 1 hour/week during pregnancy
pheno$smoke[pheno$ets_preg == 0 & pheno$smk_preg == 0] = 0
pheno$smoke = as.factor(pheno$smoke)

pheno$smk_preg = as.factor(pheno$smk_preg)

# adding 3 level variable for asthma and atopy
pheno$masth3[pheno$masthnow == 1] = 2
pheno$masth3[pheno$masthevr == 1 & is.na(pheno$masth3)] = 1
pheno$masth3[pheno$masthevr == 0] = 0
pheno$matop3[pheno$matopynow == 1] = 2
pheno$matop3[pheno$matopy == 1 & is.na(pheno$matop3)] = 1
pheno$matop3[pheno$matopy == 0] = 0

# Sample.Name = data.frame(Sample.Name = rownames(miRNAnorm1))
# miRNAnorm1 = cbind(Sample.Name, miRNAnorm1)

# transforming to negative Ccq values
miRNAnorm1[,c(2:122)] = -miRNAnorm1[,c(2:122)]

# combining microRNA and pheno data
miRNAcomb1 = merge(miRNAnorm1, pheno, on='Sample.Name')

# miRNAwithNAs = read.csv('/Users/annebozack/Documents/Lee/transfer/PRISM/breastMilkEVmiRNA/PRISM breast milk miRNA datasets for Box/PRISM_Breast Milk_EV_miRNA_QCed_wide_20200218.csv')[,-1]
# rownames(miRNAwithNAs) = miRNAwithNAs$Sample.Name
# miRNAwithNAs = miRNAwithNAs[,-1]

# miRNAwithNAsmall = miRNAwithNAs[,(colnames(miRNAwithNAs) %in% colnames(miRNAcomb))]

# miRNAcomb$NmiRNA = NA
# for(i in 1:nrow(miRNAcomb)){
	# miRNAcomb$NmiRNA[i] = table(!(is.na(as.numeric(miRNAwithNAsmall[i,]))))[2]
# }

# adding gestational weeks at breast milk collection 
age = read.csv('/Users/annebozack/Documents/Lee/transfer/PRISM/breastMilkEVmiRNA/AB20200117/AB20200117.csv')
colnames(age)[1] = 'Sample.Name'
age = rbind(age, data.frame(Sample.Name = 7006103, hhid = 70061, iid = 3, chage_breastmilk = 0.057377))

miRNAcomb1 = merge(miRNAcomb1, age, on = 'Sample.Name')

# adding a variable for low RNA used in assay
lowRNA = c(7012001, 7022602, 7023801, 7038101, 7049301, 7049201, 7047101, 7029001, 7068301, 7064001, 7018601)
miRNAcomb1$lowRNA[miRNAcomb1$Sample.Name %in% lowRNA] = 1
miRNAcomb1$lowRNA[!(miRNAcomb1$Sample.Name %in% lowRNA)] = 0
miRNAcomb1$lowRNA = factor(miRNAcomb1$lowRNA)

# removing participant with missing maternal astham/atopy data
miRNAcomb1 = miRNAcomb1[miRNAcomb1$Sample.Name != 7051401,]

# Saving normalized neg Dcq values and pheno data
write.csv(miRNAcomb1, '/Users/annebozack/Documents/Lee/miRNA asthma/manuscript/epigenomics submission/miRNA_normAllAfter35replace_pheno_gt80perc_asthma.csv')
```

## Cleaning pheno data
```{r, eval = F}
# miRNAcomb1 = read.csv('/Users/annebozack/Documents/Lee/miRNA asthma/manuscript/epigenomics submission/miRNA_normAllAfter35replace_pheno_gt80perc_asthma.csv')[,-1]

miRNAcomb1$momAge30[miRNAcomb1$age_birth < 30] = 0
miRNAcomb1$momAge30[miRNAcomb1$age_birth >= 30] = 1

miRNAcomb1$childsex = as.factor(miRNAcomb1$childsex)
miRNAcomb1$momrace3 = as.factor(miRNAcomb1$momrace3)
miRNAcomb1$edu2 = as.factor(miRNAcomb1$edu2)
miRNAcomb1$prepreg_ob = as.factor(miRNAcomb1$prepreg_ob)
miRNAcomb1$smoke = as.factor(miRNAcomb1$smoke)
miRNAcomb1$smk_preg = as.factor(miRNAcomb1$smk_preg)
miRNAcomb1$ets_preg = as.factor(miRNAcomb1$ets_preg)
miRNAcomb1$momAge30 = as.factor(miRNAcomb1$momAge30)
miRNAcomb1$masthevr = as.factor(miRNAcomb1$masthevr)
miRNAcomb1$masthnow = as.factor(miRNAcomb1$masthnow)
miRNAcomb1$matopy = as.factor(miRNAcomb1$matopy)
miRNAcomb1$matopynow = as.factor(miRNAcomb1$matopynow)
miRNAcomb1$masth3 = as.factor(miRNAcomb1$masth3)
miRNAcomb1$matop3 = as.factor(miRNAcomb1$matop3)
```

## Functions 
```{r, eval = F}
# Function for unadjusted analyses 
dfRlmUnadj = function(dfMIRNA, var, NmiRNAs){
	dfMod = data.frame(matrix(ncol = 7, nrow = NmiRNAs))
	colnames(dfMod) = c('miRNA', 'B_preg', 'p_preg', 'p_FDR_preg', 'B_ever', 'p_ever', 'p_FDR_ever')
	dfMod$miRNA = colnames(dfMIRNA[c(1:NmiRNAs)])
	for (i in 1:NmiRNAs){
		set.seed(1234)
		rlm.fit = rlm(dfMIRNA[,i] ~ dfMIRNA[[var]], method = 'MM', psi = psi.bisquare, maxit = 40)
		
		test = coeftest(rlm.fit, vcovHC(rlm.fit, type="HC0"))
		dfMod[i,2] = test[3,1]
		dfMod[i,3] = test[3,4]
		dfMod[i,5] = test[2,1]
		dfMod[i,6] = test[2,4]
	}
	dfMod$p_FDR_preg = p.adjust(dfMod$p_preg, method = 'fdr')
	dfMod$p_FDR_ever = p.adjust(dfMod$p_ever, method = 'fdr')
	return(dfMod)
}







# using rlm with default and coeftest 
# Function for analyses adjusted for sex, race, education, and week of breast milk collection
dfRlmAdjChage = function(dfMIRNA, var, NmiRNAs){
	dfMod = data.frame(matrix(ncol = 7, nrow = NmiRNAs))
	colnames(dfMod) = c('miRNA', 'B_preg', 'p_preg', 'p_FDR_preg', 'B_ever', 'p_ever', 'p_FDR_ever')
	dfMod$miRNA = colnames(dfMIRNA[c(1:NmiRNAs)])
	for (i in 1:NmiRNAs){
		set.seed(1234)
		rlm.fit = rlm(dfMIRNA[,i] ~ dfMIRNA[[var]] + dfMIRNA[['childsex']] + dfMIRNA[['momrace3']] + dfMIRNA[['edu2']] + dfMIRNA[['chage_breastmilk']], maxit = 60)
		test = coeftest(rlm.fit, vcovHC(rlm.fit, type="HC0"))
		dfMod[i,2] = test[3,1]
		dfMod[i,3] = test[3,4]
		dfMod[i,5] = test[2,1]
		dfMod[i,6] = test[2,4]
	}
	dfMod$p_FDR_preg = p.adjust(dfMod$p_preg, method = 'fdr')
	dfMod$p_FDR_ever = p.adjust(dfMod$p_ever, method = 'fdr')
	return(dfMod)
}
write.csv(dfRlmAdjChage(miRNAcomb1[,-c(1:3)], 'masth3', 130), '/Users/annebozack/Documents/Lee/miRNA asthma/manuscript/clinical epigenetics submission/rlmCoeftest.csv')

# using rlm with MM and coeftest 
# Function for analyses adjusted for sex, race, education, and week of breast milk collection
dfRlmAdjChage = function(dfMIRNA, var, NmiRNAs){
	dfMod = data.frame(matrix(ncol = 7, nrow = NmiRNAs))
	colnames(dfMod) = c('miRNA', 'B_preg', 'p_preg', 'p_FDR_preg', 'B_ever', 'p_ever', 'p_FDR_ever')
	dfMod$miRNA = colnames(dfMIRNA[c(1:NmiRNAs)])
	for (i in 1:NmiRNAs){
		set.seed(1234)
		rlm.fit = rlm(dfMIRNA[,i] ~ dfMIRNA[[var]] + dfMIRNA[['childsex']] + dfMIRNA[['momrace3']] + dfMIRNA[['edu2']] + dfMIRNA[['chage_breastmilk']], method = 'MM', maxit = 60)
		test = coeftest(rlm.fit, vcovHC(rlm.fit, type="HC0"))
		dfMod[i,2] = test[3,1]
		dfMod[i,3] = test[3,4]
		dfMod[i,5] = test[2,1]
		dfMod[i,6] = test[2,4]
	}
	dfMod$p_FDR_preg = p.adjust(dfMod$p_preg, method = 'fdr')
	dfMod$p_FDR_ever = p.adjust(dfMod$p_ever, method = 'fdr')
	return(dfMod)
}
write.csv(dfRlmAdjChage(miRNAcomb1[,-c(1:3)], 'masth3', 130), '/Users/annebozack/Documents/Lee/miRNA asthma/manuscript/clinical epigenetics submission/rlmMMCoeftest.csv')

# using rlm with MM and psi.bisquare and coeftest 
# Function for analyses adjusted for sex, race, education, and week of breast milk collection
dfRlmAdjChage = function(dfMIRNA, var, NmiRNAs){
	dfMod = data.frame(matrix(ncol = 7, nrow = NmiRNAs))
	colnames(dfMod) = c('miRNA', 'B_preg', 'p_preg', 'p_FDR_preg', 'B_ever', 'p_ever', 'p_FDR_ever')
	dfMod$miRNA = colnames(dfMIRNA[c(1:NmiRNAs)])
	for (i in 1:NmiRNAs){
		set.seed(1234)
		rlm.fit = rlm(dfMIRNA[,i] ~ dfMIRNA[[var]] + dfMIRNA[['childsex']] + dfMIRNA[['momrace3']] + dfMIRNA[['edu2']] + dfMIRNA[['chage_breastmilk']], method = 'MM', psi = psi.bisquare, maxit = 60)
		test = coeftest(rlm.fit, vcovHC(rlm.fit, type="HC0"))
		dfMod[i,2] = test[3,1]
		dfMod[i,3] = test[3,4]
		dfMod[i,5] = test[2,1]
		dfMod[i,6] = test[2,4]
	}
	dfMod$p_FDR_preg = p.adjust(dfMod$p_preg, method = 'fdr')
	dfMod$p_FDR_ever = p.adjust(dfMod$p_ever, method = 'fdr')
	return(dfMod)
}
write.csv(dfRlmAdjChage(miRNAcomb1[,-c(1:3)], 'masth3', 130), '/Users/annebozack/Documents/Lee/miRNA asthma/manuscript/clinical epigenetics submission/rlmMMpsiBisquareCoeftest.csv')

# using rlm with MM and pt 
# Function for analyses adjusted for sex, race, education, and week of breast milk collection
dfRlmAdjChage = function(dfMIRNA, var, NmiRNAs){
	dfMod = data.frame(matrix(ncol = 7, nrow = NmiRNAs))
	colnames(dfMod) = c('miRNA', 'B_preg', 'p_preg', 'p_FDR_preg', 'B_ever', 'p_ever', 'p_FDR_ever')
	dfMod$miRNA = colnames(dfMIRNA[c(1:NmiRNAs)])
	for (i in 1:NmiRNAs){
		set.seed(1234)
		rlm.fit = rlm(dfMIRNA[,i] ~ dfMIRNA[[var]] + dfMIRNA[['childsex']] + dfMIRNA[['momrace3']] + dfMIRNA[['edu2']] + dfMIRNA[['chage_breastmilk']], method = 'MM', maxit = 60)
		test = data.frame(summary(rlm.fit)$coefficients)
		test$p = 2*pt(abs(test$t.value), summary(rlm.fit)$df[2], lower.tail=F)
		dfMod[i,2] = test[3,1]
		dfMod[i,3] = test[3,4]
		dfMod[i,5] = test[2,1]
		dfMod[i,6] = test[2,4]
	}
	dfMod$p_FDR_preg = p.adjust(dfMod$p_preg, method = 'fdr')
	dfMod$p_FDR_ever = p.adjust(dfMod$p_ever, method = 'fdr')
	return(dfMod)
}
write.csv(dfRlmAdjChage(miRNAcomb1[,-c(1:3)], 'masth3', 130), '/Users/annebozack/Documents/Lee/miRNA asthma/manuscript/clinical epigenetics submission/rlmMMpt.csv')


# using lmrob 
# Function for analyses adjusted for sex, race, education, and week of breast milk collection
library(robustbase)
dfRlmAdjChage = function(dfMIRNA, var, NmiRNAs){
	dfMod = data.frame(matrix(ncol = 7, nrow = NmiRNAs))
	colnames(dfMod) = c('miRNA', 'B_preg', 'p_preg', 'p_FDR_preg', 'B_ever', 'p_ever', 'p_FDR_ever')
	dfMod$miRNA = colnames(dfMIRNA[c(1:NmiRNAs)])
	for (i in 1:NmiRNAs){
		set.seed(1234)
		lmrob.fit = lmrob(dfMIRNA[,i] ~ dfMIRNA[[var]] + dfMIRNA[['childsex']] + dfMIRNA[['momrace3']] + dfMIRNA[['edu2']] + dfMIRNA[['chage_breastmilk']], k.max = 400)
		test = data.frame(summary(lmrob.fit)$coefficients)
		dfMod[i,2] = test[3,1]
		dfMod[i,3] = test[3,4]
		dfMod[i,5] = test[2,1]
		dfMod[i,6] = test[2,4]
	}
	dfMod$p_FDR_preg = p.adjust(dfMod$p_preg, method = 'fdr')
	dfMod$p_FDR_ever = p.adjust(dfMod$p_ever, method = 'fdr')
	return(dfMod)
}
write.csv(dfRlmAdjChage(miRNAcomb1[,-c(1:3)], 'masth3', 130), '/Users/annebozack/Documents/Lee/miRNA asthma/manuscript/clinical epigenetics submission/lmrob.csv')


# using lm 
# Function for analyses adjusted for sex, race, education, and week of breast milk collection
dfRlmAdjChage = function(dfMIRNA, var, NmiRNAs){
	dfMod = data.frame(matrix(ncol = 7, nrow = NmiRNAs))
	colnames(dfMod) = c('miRNA', 'B_preg', 'p_preg', 'p_FDR_preg', 'B_ever', 'p_ever', 'p_FDR_ever')
	dfMod$miRNA = colnames(dfMIRNA[c(1:NmiRNAs)])
	for (i in 1:NmiRNAs){
		set.seed(1234)
		lm.fit = lm(dfMIRNA[,i] ~ dfMIRNA[[var]] + dfMIRNA[['childsex']] + dfMIRNA[['momrace3']] + dfMIRNA[['edu2']] + dfMIRNA[['chage_breastmilk']])
		test = data.frame(summary(lm.fit)$coefficients)
		dfMod[i,2] = test[3,1]
		dfMod[i,3] = test[3,4]
		dfMod[i,5] = test[2,1]
		dfMod[i,6] = test[2,4]
	}
	dfMod$p_FDR_preg = p.adjust(dfMod$p_preg, method = 'fdr')
	dfMod$p_FDR_ever = p.adjust(dfMod$p_ever, method = 'fdr')
	return(dfMod)
}
write.csv(dfRlmAdjChage(miRNAcomb1[,-c(1:3)], 'masth3', 130), '/Users/annebozack/Documents/Lee/miRNA asthma/manuscript/clinical epigenetics submission/lm.csv')




# using coeftest
# Function for analyses adjusted for sex, race, education, and week of breast milk collection
dfRlmAdjChage = function(dfMIRNA, var, NmiRNAs){
	dfMod = data.frame(matrix(ncol = 7, nrow = NmiRNAs))
	colnames(dfMod) = c('miRNA', 'B_preg', 'p_preg', 'p_FDR_preg', 'B_ever', 'p_ever', 'p_FDR_ever')
	dfMod$miRNA = colnames(dfMIRNA[c(1:NmiRNAs)])
	for (i in 1:NmiRNAs){
		set.seed(1234)
		rlm.fit = rlm(dfMIRNA[,i] ~ dfMIRNA[[var]] + dfMIRNA[['childsex']] + dfMIRNA[['momrace3']] + dfMIRNA[['edu2']] + dfMIRNA[['chage_breastmilk']], method = 'MM', psi = psi.bisquare, maxit = 60)
		test = coeftest(rlm.fit, vcovHC(rlm.fit, type="HC0"))
		dfMod[i,2] = test[3,1]
		dfMod[i,3] = test[3,4]
		dfMod[i,5] = test[2,1]
		dfMod[i,6] = test[2,4]
	}
	dfMod$p_FDR_preg = p.adjust(dfMod$p_preg, method = 'fdr')
	dfMod$p_FDR_ever = p.adjust(dfMod$p_ever, method = 'fdr')
	return(dfMod)
}

# using pt
# Function for analyses adjusted for sex, race, education, and week of breast milk collection
dfRlmAdjChage = function(dfMIRNA, var, NmiRNAs){
	dfMod = data.frame(matrix(ncol = 7, nrow = NmiRNAs))
	colnames(dfMod) = c('miRNA', 'B_preg', 'p_preg', 'p_FDR_preg', 'B_ever', 'p_ever', 'p_FDR_ever')
	dfMod$miRNA = colnames(dfMIRNA[c(1:NmiRNAs)])
	for (i in 1:NmiRNAs){
		set.seed(1234)
		rlm.fit = rlm(dfMIRNA[,i] ~ dfMIRNA[[var]] + dfMIRNA[['childsex']] + dfMIRNA[['momrace3']] + dfMIRNA[['edu2']] + dfMIRNA[['chage_breastmilk']], method = 'MM', psi = psi.bisquare, maxit = 60)
		test = data.frame(summary(rlm.fit)$coefficients)
		test = coef$p = 2*pt(abs(coef$t.value), summary(mod)$df[2], lower.tail=F)
		dfMod[i,2] = test[3,1]
		dfMod[i,3] = test[3,4]
		dfMod[i,5] = test[2,1]
		dfMod[i,6] = test[2,4]
	}
	dfMod$p_FDR_preg = p.adjust(dfMod$p_preg, method = 'fdr')
	dfMod$p_FDR_ever = p.adjust(dfMod$p_ever, method = 'fdr')
	return(dfMod)
}




# # Function for analyses adjusted for sex, race, education, and week of breast milk collection with sex interaction term
# dfRlmAdjChageInt = function(dfMIRNA, var, NmiRNAs){
	# dfMod = data.frame(matrix(ncol = 7, nrow = NmiRNAs))
	# colnames(dfMod) = c('miRNA', 'B', 'p', 'B_sex', 'p_sex', 'B_sexInt', 'p_sexInt')
	# dfMod$miRNA = colnames(dfMIRNA[c(1:NmiRNAs)])
	# for (i in 1:NmiRNAs){
		# set.seed(1234)
		# mod = rlm(dfMIRNA[,i] ~ dfMIRNA[[var]]*dfMIRNA[['childsex']] + dfMIRNA[['momrace3']] + dfMIRNA[['edu2']] + dfMIRNA[['chage_breastmilk']], method = 'MM', psi = psi.bisquare, maxit = 60)
		# coef = data.frame(summary(mod)$coefficients)
		# coef$p = 2*pt(abs(coef$t.value), summary(mod)$df[2], lower.tail=F)
		# dfMod[i,2] = coef[2,1]
		# dfMod[i,3] = coef[2,4]
		# dfMod[i,4] = coef[3,1]
		# dfMod[i,5] = coef[3,4]
		# dfMod[i,6] = coef[8,1]
		# dfMod[i,7] = coef[8,4]
	# }
	# return(dfMod)
# }

# Lambda 
lambda <- function(p) median(qchisq(p, df=1, lower.tail=FALSE), na.rm=TRUE) / qchisq(0.5, df=1)

# QQ plot
gg_qqplot = function(pvector) {
	l = round(lambda(pvector), 3)
	o = -log10(sort(pvector, decreasing = FALSE))
	e = -log10(ppoints(length(pvector)))
	df = data.frame(o = o, e = e)
	ggplot(df, aes(e, o)) + geom_point(alpha = 0.5, size = 1) + geom_abline(intercept = 0, slope = 1, color = '#AB3428') + labs(y = expression(Observed ~ ~-log[10](italic(p))), x = expression(Expected ~ ~-log[10](italic(p)))) + theme_classic() + annotate("text", x = 0.5, y = 2.5, label = paste0('lambda = ', l))
}
```

## Matheral asthma
### Unadjusted 
```{r, eval = F}
asthUnadj = dfRlmUnadj(miRNAcomb1[,-c(1:3)], 'masth3', 130)

# Active during pregnancy vs. never
table(asthUnadj$p_preg <0.05)
# FALSE  TRUE 
  # 103    18 
  
table(asthUnadj$p_FDR_preg <0.2)
# FALSE  TRUE 
  # 111    10

table(asthUnadj$p_FDR_preg <0.1)
# FALSE 
# 121 

table(asthUnadj$p_FDR_preg <0.05)
# FALSE 
# 121 

# Inactive vs. never
table(asthUnadj$p_ever <0.05)
# FALSE  TRUE 
  # 114     7 
  
table(asthUnadj$p_FDR_ever <0.2)
# FALSE 
# 121 

table(asthUnadj$p_FDR_ever <0.1)
# FALSE 
# 121 

table(asthUnadj$p_FDR_ever <0.05)
# FALSE 
# 121 
```


### Adjusted for sex, race, education, and week of breast milk collection
```{r, eval = F}
asthAdj = dfRlmAdjChage(miRNAcomb1[,-c(1:3)], 'masth3', 121)

# Active during pregnancy vs. never
table(asthAdj$p_preg < 0.05)
# FALSE  TRUE 
  # 106    15 

table(asthAdj$p_FDR_preg < 0.2)
# FALSE 
  # 121 

table(asthAdj$p_FDR_preg < 0.1)
# FALSE 
# 121 
 
table(asthAdj$p_FDR_preg < 0.05)
# FALSE 
# 121  

# Inactive vs. never
table(asthAdj$p_ever < 0.05)
# FALSE  TRUE 
  # 109    12

table(asthAdj$p_FDR_ever < 0.2)
# FALSE  TRUE 
  # 117     4 
  
table(asthAdj$p_FDR_ever < 0.1)
# FALSE  TRUE 
  # 117     4 
 
table(asthAdj$p_FDR_ever < 0.05)
# FALSE  TRUE 
  # 117     4  

```

```{r, warning=FALSE,message=FALSE,eval = F}
# correcting microRNA names
miRBase = read.csv('/Users/annebozack/Documents/Lee/transfer/PRISM/breastMilkEVmiRNA/lab methods/miRBaseIDs.csv', colClasses = c(rep("character",4)))

miRBase$miRNA = paste0(miRBase$AssayName, '_', miRBase$AssayID)
miRBase$miRNA = str_replace_all(miRBase$miRNA, '[-#]', '.')

asthAdj = merge(asthAdj, miRBase[,c(5,3, 4)], by = 'miRNA')

sigAsthAdj_preg = asthAdj[asthAdj$p_preg <0.05, ]
sigAsthAdj_ever = asthAdj[asthAdj$p_ever <0.05, ]

asthAdj$logp_preg = -log(asthAdj$p_preg, 10)
asthAdj$sig_preg[asthAdj$p_preg < 0.05] = 1
asthAdj$sig_preg[!(asthAdj$p_preg < 0.05)] = 0
asthAdj$logp_ever = -log(asthAdj$p_ever, 10)
asthAdj$sig_ever[asthAdj$p_ever < 0.05] = 1
asthAdj$sig_ever[!(asthAdj$p_ever < 0.05)] = 0

```

#### Active during pregnancy
```{r, echo = F}
sigAsthAdj_preg[order(sigAsthAdj_preg$p_preg),] %>% kable() %>% kable_styling(font_size = 12) 
```

```{r, warning=FALSE,message=FALSE, out.width = '75%', fig.asp = 1}
ggplot(asthAdj, aes(x = B_preg, y = logp_preg, colour = factor(sig_preg))) + geom_point(size = 1.5, alpha = 0.5) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0, size = 0.5) + geom_hline(yintercept=-log(0.05, 10), colour = '#b30000', linetype='dashed', size = 0.25) + geom_vline(xintercept=-0.2, colour = '#b30000', linetype='dashed', size = 0.25) + geom_vline(xintercept=0.2, colour = '#b30000', linetype='dashed', size = 0.25) + ylab(expression('-log'['10']~'('~italic('p')~')')) + xlab(expression(~italic("B"))) + annotate('text', x = 0.75, y = -log(0.045, 10), label = '~italic(p)== 0.05', parse = T, size=3.5) + scale_colour_manual(values = c('black', '#b35806')) + theme_minimal() + theme(text = element_text(size=10)) + theme(legend.position = "none") + geom_hline(yintercept=-log(0.05/130, 10), colour = '#b30000', size = 0.25)
```

#### Inactive
```{r, echo = F}
sigAsthAdj_ever[order(sigAsthAdj_ever$p_ever),] %>% kable() %>% kable_styling(font_size = 12) 
```

```{r, warning=FALSE,message=FALSE, out.width = '75%', fig.asp = 1}
ggplot(asthAdj, aes(x = B_ever, y = logp_ever, colour = factor(sig_ever))) + geom_point(size = 1.5, alpha = 0.5) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0, size = 0.5) + geom_hline(yintercept=-log(0.05, 10), colour = '#b30000', linetype='dashed', size = 0.25) + geom_vline(xintercept=-0.2, colour = '#b30000', linetype='dashed', size = 0.25) + geom_vline(xintercept=0.2, colour = '#b30000', linetype='dashed', size = 0.25) + ylab(expression('-log'['10']~'('~italic('p')~')')) + xlab(expression(~italic("B"))) + annotate('text', x = 0.75, y = -log(0.045, 10), label = '~italic(p)== 0.05', parse = T, size=3.5) + scale_colour_manual(values = c('black', '#b35806')) + theme_minimal() + theme(text = element_text(size=10)) + theme(legend.position = "none") + geom_hline(yintercept=-log(0.05/130, 10), colour = '#b30000', size = 0.25)

```


## Maternal atopy
### Unadjusted 
```{r, eval = F}
atopUnadj = dfRlmUnadj(miRNAcomb1[,-c(1:3)], 'matop3', 121)

# Active during pregnancy vs. never
table(atopUnadj$p_preg <0.05)
# FALSE  TRUE 
  # 119     2 
  
table(atopUnadj$p_FDR_preg <0.2)
# FALSE 
# 121 

table(atopUnadj$p_FDR_preg <0.1)
# FALSE 
# 121 

table(atopUnadj$p_FDR_preg <0.05)
# FALSE 
# 121 

# Inactive vs. never
table(atopUnadj$p_ever <0.05)
# FALSE  TRUE 
  # 114     7 
  
table(atopUnadj$p_FDR_ever <0.2)
# FALSE 
# 121 

table(atopUnadj$p_FDR_ever <0.1)
# FALSE 
# 121 

table(atopUnadj$p_FDR_ever <0.05)
# FALSE 
# 121 
```

### Adjusted for sex, race, education, and week of breast milk collection
```{r, eval = F}
atopAdj = dfRlmAdjChage(miRNAcomb1[,-c(1:3)], 'matop3', 121)

# Active during pregnancy vs. never
table(atopAdj$p_preg < 0.05)
# FALSE  TRUE 
  # 117     4 

table(atopAdj$p_FDR_preg < 0.2)
# FALSE 
  # 121 

table(atopAdj$p_FDR_preg < 0.1)
# FALSE 
# 121 
 
table(atopAdj$p_FDR_preg < 0.05)
# FALSE 
# 121  

# Inactive vs. never
table(atopAdj$p_ever < 0.05)
# FALSE  TRUE 
  # 120     1 

table(atopAdj$p_FDR_ever < 0.2)
# FALSE 
  # 121  
  
table(atopAdj$p_FDR_ever < 0.1)
# FALSE 
  # 121  
 
table(atopAdj$p_FDR_ever < 0.05)
# FALSE 
  # 121  

```

```{r, warning=FALSE,message=FALSE,eval = F}
# correcting microRNA names
atopAdj = merge(atopAdj, miRBase[,c(5,3, 4)], by = 'miRNA')

sigAtopAdj_preg = atopAdj[atopAdj$p_preg <0.05, ]
sigAtopAdj_ever = atopAdj[atopAdj$p_ever <0.05, ]

atopAdj$logp_preg = -log(atopAdj$p_preg, 10)
atopAdj$sig_preg[atopAdj$p_preg < 0.05] = 1
atopAdj$sig_preg[!(atopAdj$p_preg < 0.05)] = 0
atopAdj$logp_ever = -log(atopAdj$p_ever, 10)
atopAdj$sig_ever[atopAdj$p_ever < 0.05] = 1
atopAdj$sig_ever[!(atopAdj$p_ever < 0.05)] = 0

```

#### Active during pregnancy
```{r, echo = F}
sigAtopAdj_preg[order(sigAtopAdj_preg$p_preg),] %>% kable() %>% kable_styling(font_size = 12) 
```

```{r, warning=FALSE,message=FALSE, out.width = '75%', fig.asp = 1}
ggplot(atopAdj, aes(x = B_preg, y = logp_preg, colour = factor(sig_preg))) + geom_point(size = 1.5, alpha = 0.5) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0, size = 0.5) + geom_hline(yintercept=-log(0.05, 10), colour = '#b30000', linetype='dashed', size = 0.25) + geom_vline(xintercept=-0.2, colour = '#b30000', linetype='dashed', size = 0.25) + geom_vline(xintercept=0.2, colour = '#b30000', linetype='dashed', size = 0.25) + ylab(expression('-log'['10']~'('~italic('p')~')')) + xlab(expression(~italic("B"))) + annotate('text', x = 0.75, y = -log(0.045, 10), label = '~italic(p)== 0.05', parse = T, size=3.5) + scale_colour_manual(values = c('black', '#b35806')) + theme_minimal() + theme(text = element_text(size=10)) + theme(legend.position = "none")
```

#### Inactive
```{r, echo = F}
sigAtopAdj_ever[order(sigAtopAdj_ever$p_ever),] %>% kable() %>% kable_styling(font_size = 12) 
```

```{r, warning=FALSE,message=FALSE, out.width = '75%', fig.asp = 1}
ggplot(atopAdj, aes(x = B_ever, y = logp_ever, colour = factor(sig_ever))) + geom_point(size = 1.5, alpha = 0.5) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0, size = 0.5) + geom_hline(yintercept=-log(0.05, 10), colour = '#b30000', linetype='dashed', size = 0.25) + geom_vline(xintercept=-0.2, colour = '#b30000', linetype='dashed', size = 0.25) + geom_vline(xintercept=0.2, colour = '#b30000', linetype='dashed', size = 0.25) + ylab(expression('-log'['10']~'('~italic('p')~')')) + xlab(expression(~italic("B"))) + annotate('text', x = 0.75, y = -log(0.045, 10), label = '~italic(p)== 0.05', parse = T, size=3.5) + scale_colour_manual(values = c('black', '#b35806')) + theme_minimal() + theme(text = element_text(size=10)) + theme(legend.position = "none")
```








