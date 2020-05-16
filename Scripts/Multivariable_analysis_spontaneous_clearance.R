# Created by Marc Niebel Feb 2020
# Multivariable analysis for clinical outcome, spontaneous clearance
#using a cox proportional hazards model

library(survival)
library(survminer)
library(ggplot2)
library(dplyr)
library(tidyr)

#Dataframe which was created from Survivial analysis_v2.R
data_mv <- read.csv("./data/multivariable_df.csv")

#Remember P66 has been identified as being dually infected and thus has two lines.
#Will be removed for downstream analysis
#clean_data <- data %>%
#    filter(record_id != 66)

#STEP 1: PRIMARY ANALYSIS 
#Include variables which were main remit of study
#Use variables for which association is known
# age, gender, HIV status,genotype, peak Bilirubin, peak ALT
# Other ones that are known are PWID and MSM
#Diversity is missing and IL28B

columns_to_convert_to_factor <- c("hiv","Diabetes","alco_excess","PWID","MSM","ARVs")
data_mv[columns_to_convert_to_factor] <- lapply(data_mv[columns_to_convert_to_factor],factor)
data_mv <- rename(data_mv,Gender=gender.factor)
data_mv$Gender <- relevel(data_mv$Gender,ref = "Male")
data_mv$hiv <- relevel(data_mv$hiv,ref = "1")
data_mv$MSM <- relevel(data_mv$MSM,ref = "1")

multivariable_cox <- coxph(Surv(Time,Event)~age+Gender+hiv+peak_Bil_binary+Binary_peakALT+MSM+PWID+Genotype,data=data_mv)
summary(multivariable_cox)
#onefile=FALSE is used to ensure first page not being blank
pdf(file="Forest plot spontaneous clearance.pdf", onefile = FALSE)
#Forest plot to summarise data
ggforest(multivariable_cox,data = data_mv)
dev.off()

#Looking at missing data for spontaneous clearers currently
missing_data_sc <- multivariable_df %>% filter(Event==1)
sapply(missing_data_sc, function(x) sum(is.na(x)))
#Ones of relevance are Genotype n=10 and Il28B n=12 which have missing data
#Problem is that patients who spontaneously clear cannot be typed

#Remove gender
multivariable_cox_remove_gender <- coxph(Surv(Time,Event)~age+hiv+peak_Bil_binary+Binary_peakALT+MSM+PWID+Genotype,data=data_mv)
summary(multivariable_cox_remove_gender)
ggforest(multivariable_cox_remove_gender,data=data_mv)

#Remove gender and genotype
multivariable_cox_remove_gender_genotype <- coxph(Surv(Time,Event)~age+hiv+peak_Bil_binary+Binary_peakALT+MSM+PWID,data=data_mv)
summary(multivariable_cox_remove_gender_genotype)
ggforest(multivariable_cox_remove_gender_genotype,data=data_mv)

#Remove genotype
multivariable_cox_remove_genotype <- coxph(Surv(Time,Event)~age+Gender+hiv+peak_Bil_binary+Binary_peakALT+MSM+PWID,data=data_mv)
summary(multivariable_cox_remove_genotype)
ggforest(multivariable_cox_remove_genotype,data=data_mv)

#MSM and Genotype removed
multivariable_cox_remove_MSM_Genotype <- coxph(Surv(Time,Event)~age+Gender+hiv+peak_Bil_binary+Binary_peakALT+PWID,data=data_mv)
summary(multivariable_cox_remove_MSM_Genotype)
ggforest(multivariable_cox_remove_MSM_Genotype,data=data_mv)

#MSM removed
multivariable_cox_remove_MSM <- coxph(Surv(Time,Event)~age+Gender+hiv+peak_Bil_binary+Binary_peakALT+PWID+Genotype,data=data_mv)
summary(multivariable_cox_remove_MSM_Genotype)
ggforest(multivariable_cox_remove_MSM,data=data_mv)

adding_column_both_risk_factors <- data_mv %>% 
        mutate(MSM_PWID=case_when(MSM == 1 & PWID == 1 ~ "Both_risk",TRUE~"One_risk")) %>%
            filter(MSM_PWID=="One_risk")

#Removed both risk factors together
cox_removed_both_risk_factors <- coxph(Surv(Time,Event)~age+Gender+hiv+peak_Bil_binary+Binary_peakALT+MSM+PWID+Genotype,data=adding_column_both_risk_factors)
ggforest(cox_removed_both_risk_factors,data=adding_column_both_risk_factors)
cox_removed_both_risk_factors_genotype <- coxph(Surv(Time,Event)~age+Gender+hiv+peak_Bil_binary+Binary_peakALT+MSM+PWID,data=adding_column_both_risk_factors)
ggforest(cox_removed_both_risk_factors_genotype,data=adding_column_both_risk_factors)

#HIV postive patients with gender removed cause females are not able to be MSM(mutually exclusive)
hiv_positive_patients <- data_mv %>% filter(hiv==1)
sapply(hiv_positive_patients, function(x) sum(is.na(x)))
multivariable_cox_hiv_positive <- coxph(Surv(Time,Event)~age+peak_Bil_binary+Binary_peakALT+MSM+PWID+Genotype+CD4_count+ARVs,data=hiv_positive_patients)
summary(multivariable_cox_hiv_positive)
ggforest(multivariable_cox_hiv_positive,data=hiv_positive_patients)

#HIV positive patients with ARVs removed
multivariable_cox_hiv_positive_arvs_removed <- coxph(Surv(Time,Event)~age+peak_Bil_binary+Binary_peakALT+MSM+PWID+Genotype+CD4_count,data=hiv_positive_patients)
summary(multivariable_cox_hiv_positive_arvs_removed)
ggforest(multivariable_cox_hiv_positive_arvs_removed,data=hiv_positive_patients)

#DIAGNOSTICS
proportional_assumption <- cox.zph(multivariable_cox)
#Non-significant relationship between residuals and time are refuted by a significant relationship
ggcoxzph(proportional_assumption)
#None were significant
#NEED MORE THOROUGH ANALYSIS
ggcoxdiagnostics(multivariable_cox,type="dfbeta",linear.predictions = FALSE,ggtheme = theme_bw())
ggcoxdiagnostics(multivariable_cox,type="deviance",linear.predictions = FALSE,ggtheme=theme_bw())
#Testing non-linearity for continous variables(IN PROGRESS)





