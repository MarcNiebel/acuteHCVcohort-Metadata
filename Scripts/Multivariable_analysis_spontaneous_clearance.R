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

#STEP 1: PRIMARY ANALYSIS 
#Remit of study: age, ?gender,genotype, HIV-status,peak Bilirubin, peakALT,diversity(missing)
#Use variables(s) for which association is well known
#PWID

columns_to_convert_to_factor <- c("hiv","Diabetes","alco_excess","PWID","MSM","ARVs")
data_mv[columns_to_convert_to_factor] <- lapply(data_mv[columns_to_convert_to_factor],factor)
data_mv <- rename(data_mv,Gender=gender.factor)
data_mv <-rename(data_mv,IL28B=CC_NonCC)
data_mv$Gender <- relevel(data_mv$Gender,ref = "Male")
data_mv$hiv <- relevel(data_mv$hiv,ref = "1")
data_mv$MSM <- relevel(data_mv$MSM,ref = "1")
data_mv$IL28B <- relevel(data_mv$IL28B,ref="CT/TT")

multivariable_cox <- coxph(Surv(Time,Event)~age+Gender+hiv+peak_Bil_binary+Binary_peakALT+PWID+Genotype,data=data_mv)
sink("Output/Multivariable_sc/all_variables_of_interest.txt")
print(multivariable_cox)
sink(file=NULL)

#onefile=FALSE is used to ensure first page not being blank
pdf(file="Output/Multivariable_sc/Forest_plot_sc_all_variables.pdf", onefile = FALSE)
#Forest plot to summarise data
ggforest(multivariable_cox)
dev.off()

#Remove gender cause only 11 females
multivariable_cox_remove_gender <-update(multivariable_cox,~.-Gender)
sink("Output/Multivariable_sc/gender_removed.txt")
print(multivariable_cox_remove_gender)
sink(file=NULL)
pdf(file="Output/Multivariable_sc/Forest_plot_excluded_gender.pdf",onefile = FALSE)
ggforest(multivariable_cox_remove_gender)
dev.off()

#No difference in either model
anova(multivariable_cox,multivariable_cox_remove_gender)

data_combine_genotype <- data_mv %>% 
    mutate(Gt1_nonGt1=case_when(Genotype=="gt1a" ~ "gt1",TRUE~"non-gt1"))
multivariable_cox_combine <- coxph(Surv(Time,Event)~age+Gender+hiv+peak_Bil_binary+Binary_peakALT+PWID+Gt1_nonGt1,data=data_combine_genotype)
pdf(file="Output/Multivariable_sc/Forest_plot_combined_genotype.pdf",onefile = FALSE)
ggforest(multivariable_cox_combine)
dev.off()

multivariable_cox_combine_remove_gender <- update(multivariable_cox_combine,~.-Gender)
summary(multivariable_cox_combine_remove_gender)
pdf(file="Output/Multivariable_sc/Forest_plot_combined_genotype_remove_gender.pdf",onefile = FALSE)
ggforest(multivariable_cox_combine_remove_gender)
dev.off()

#No difference in the two models
anova(multivariable_cox_combine,multivariable_cox_combine_remove_gender)

#Exploratory analysis
#Not including gender(not possible with MSM) but including MSM and IL28B(note only 111 typed)
multivariable_cox_exploratory <- update(multivariable_cox,~. -Gender+MSM+IL28B)
multivariable_cox_exploratory_v2 <- update(multivariable_cox_combine,~. -Gender+MSM+IL28B)

#Concordance index quantifies the level of model fit which could be used for predictions.
#Briefly spoken, the -index can be interpreted as the probability that a patient with a 
#small survival time is associated with a high value of a biomarker combination (and vice versa). 
#Consequently, it measures the concordance between the rankings of the survival times and the biomarker 
#values and therefore the ability of a biomarker to discriminate between patients with small survival times 
#and patients with large survival times (Andreas Mayr. et al Plos One 2014)

#DIAGNOSTICS
proportional_assumption <- cox.zph(multivariable_cox)
#Non-significant relationship between residuals and time are refuted by a significant relationship
ggcoxzph(proportional_assumption)
#None were significant
#NEED MORE THOROUGH ANALYSIS
ggcoxdiagnostics(multivariable_cox,type="dfbeta",linear.predictions = FALSE,ggtheme = theme_bw())
ggcoxdiagnostics(multivariable_cox,type="deviance",linear.predictions = FALSE,ggtheme=theme_bw())
#Testing non-linearity for continous variables(IN PROGRESS)





