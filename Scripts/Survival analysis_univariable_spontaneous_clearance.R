#Created by Marc Niebel November 2019
#Kaplan Meier analysis and cox regression based on endpoint of spontaneous clearance
#Producing dataframe for multivariable analysis

library(survival)
library(survminer)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(broom)

###################################################
#Sourcing the dataframe from the data folder
source("./data/AcuteHCV_R_2020-04-22_1120.r", chdir = TRUE)
#Read in the time to event .csv with the status of event
time_to_event <-read.csv("./data/time_to_event.csv")
#Combining Redcap data and time to event variable
total_data <-merge(data,time_to_event,by="record_id")
#Removing chronic patients
total_data <- total_data %>%
    filter(record_id != 35 & record_id != 213 & record_id != 41)
##################################################
spont_clearance_time_gender <- total_data %>% select(record_id,Event,Time,gender.factor)
spont_clearance_time_gender$gender.factor <-relevel(factor(spont_clearance_time_gender$gender.factor),ref="Male")
gender_KMplot <- ggsurvplot(survfit(Surv(Time,Event)~gender.factor,data=spont_clearance_time_gender),
           xlab="Days",
           ylab="Proportion of HCV persistence",
           legend.title="Sex",
           legend.labs=c("Male","Female"),
           pval = TRUE,
           conf.int = TRUE)
ggsave(file="Output/survival analysis_uni_sc/gender_KM.pdf",print(gender_KMplot),onefile=FALSE)
#Looking at how this has changed at specific points in time
summary(survfit(Surv(Time,Event)~gender.factor,data=spont_clearance_time_gender),times=182.5)
summary(survfit(Surv(Time,Event)~gender.factor,data=spont_clearance_time_gender),times=365.25)
cox_gender_co <-coxph(Surv(Time,Event)~gender.factor,data=spont_clearance_time_gender)
sink("Output/survival analysis_uni_sc/cox_gender.txt")
print(cox_gender_co)
sink(file=NULL)
#Proportional hazards assumption
cox_assumption_gender <- cox.zph(cox_gender_co)
schoenfield_gender <- ggcoxzph(cox_assumption_gender)
ggsave("Output/survival analysis_uni_sc/schoenfield_gender.pdf", arrangeGrob(grobs=schoenfield_gender))
########################################
spont_clearance_time_age <- total_data %>% select(record_id,Time,Event,age)
cox_age_co <-coxph(Surv(Time,Event)~age,data=spont_clearance_time_age)
sink("Output/survival analysis_uni_sc/cox_age.txt")
print(cox_age_co)
sink(file=NULL)
#Proportional hazards assumption
cox_assumption_age <- cox.zph(cox_age_co)
schoenfield_age <- ggcoxzph(cox_assumption_age)
ggsave("Output/survival analysis_uni_sc/schoenfield_age.pdf", arrangeGrob(grobs=schoenfield_age))
########################################
spont_clearance_time_HIVstatus <- total_data %>% select(record_id,Event,Time,hiv)
spont_clearance_time_HIVstatus$hiv <- relevel(factor(spont_clearance_time_HIVstatus$hiv),ref="1")
age_KMplot <- ggsurvplot(survfit(Surv(Time,Event)~hiv,data=spont_clearance_time_HIVstatus),
           xlab="Days",
           ylab="Proportion of HCV persistence",
           legend.title="HIV Status",
           legend.labs=c("HIV negative","HIV positive"),
           pval = TRUE,
           conf.int = TRUE)
ggsave(file="Output/survival analysis_uni_sc/age_KM.pdf",print(age_KMplot),onefile=FALSE)
summary(survfit(Surv(Time,Event)~hiv,data=spont_clearance_time_HIVstatus),times=182.5)
summary(survfit(Surv(Time,Event)~hiv,data=spont_clearance_time_HIVstatus),times=365.25)
cox_hiv_co <- coxph(Surv(Time,Event)~hiv,data=spont_clearance_time_HIVstatus)
sink("Output/survival analysis_uni_sc/cox_hiv.txt")
print(cox_hiv_co)
sink(file=NULL)
#Proportional hazards assumption
cox_assumption_hiv <- cox.zph(cox_hiv_co)
schoenfield_hiv <- ggcoxzph(cox_assumption_hiv)
ggsave("Output/survival analysis_uni_sc/schoenfield_hiv.pdf", arrangeGrob(grobs=schoenfield_hiv))
########################################
spont_clearance_time_peakALT <-total_data %>%select(record_id,Event,Time,alt_peak)
clean_spont_clearance_time_peakALT <-na.omit(spont_clearance_time_peakALT)
peakALT_ClinicalOutcome_binary <- clean_spont_clearance_time_peakALT %>% mutate(Binary_peakALT=case_when(alt_peak >=1000 ~ ">1000",TRUE~"<1000"))
peakALT_ClinicalOutcome_binary <- peakALT_ClinicalOutcome_binary[,-4]
peakALT_KMplot <- ggsurvplot(survfit(Surv(Time,Event)~Binary_peakALT,data=peakALT_ClinicalOutcome_binary),
           xlab="Days",
           ylab="Proportion of HCV persistence",
           legend.title="Peak ALT",
           legend.labs=c("peak ALT<1000","peak ALT>1000"),
           pval = TRUE,
           conf.int = TRUE)
ggsave(file="Output/survival analysis_uni_sc/peakALT_KM.pdf",print(peakALT_KMplot),onefile=FALSE)
cox_peakALT_co <-coxph(Surv(Time,Event)~Binary_peakALT,data=peakALT_ClinicalOutcome_binary)
sink("Output/survival analysis_uni_sc/cox_peakALT.txt")
print(cox_peakALT_co)
sink(file=NULL)
#Proportional hazards assumption
cox_assumption_peakALT <-cox.zph(cox_peakALT_co)
schoenfield_peakALT <- ggcoxzph(cox_assumption_peakALT)
ggsave("Output/survival analysis_uni_sc/schoenfield_peakALT.pdf", arrangeGrob(grobs=schoenfield_hiv))
########################################
spont_clearance_time_drug_use <- total_data %>% select(record_id,Event,Time,drugs___20:drugs___19)
all_Drug_use <- spont_clearance_time_drug_use %>% 
    mutate(Drug_use=case_when(drugs___20 == 1|drugs___21 == 1|
                              drugs___22 == 1|drugs___1 == 1 |
                              drugs___2 ==1 | drugs___3 == 1 |
                              drugs___4 == 1 | drugs___5 == 1|
                              drugs___6 ==1 | drugs___7 ==1|
                              drugs___8 == 1  |drugs___9 == 1| 
                              drugs___10 == 1 |drugs___11 == 1|
                              drugs___12 == 1| drugs___13 ==1|
                              drugs___14 ==1 | drugs___15 == 1|
                              drugs___16 == 1 |drugs___17 == 1|
                              drugs___18 ==1 |drugs___19 == 1
                              ~ "Drug use",TRUE ~ "no use"))
all_Drug_use <-all_Drug_use[c(1:3,26)]
all_Drug_use$Drug_use <-relevel(factor(all_Drug_use$Drug_use),ref = "no use")
druguse_KMplot <- ggsurvplot(survfit(Surv(Time,Event)~Drug_use,data=all_Drug_use),
           xlab="Days",
           ylab="Proportion of HCV persistence",
           legend.title="Drug Use",
           legend.labs=c("No Drug Use","Drug Use"),
           pval = TRUE,
           conf.int = TRUE)
ggsave(file="Output/survival analysis_uni_sc/Drug_Use_KM.pdf",print(druguse_KMplot),onefile=FALSE)
cox_drug_use_co <- coxph(Surv(Time,Event)~Drug_use,data=all_Drug_use)
sink("Output/survival analysis_uni_sc/cox_drug_use.txt")
print(cox_drug_use_co)
sink(file=NULL)
#Proportional hazards assumption
cox_assumption_druguse <- cox.zph(cox_drug_use_co)
schoenfield_druguse <- ggcoxzph(cox_assumption_druguse)
ggsave("Output/survival analysis_uni_sc/schoenfield_druguse.pdf", arrangeGrob(grobs=schoenfield_druguse))
########################################
spont_clearance_time_cocaine_use <-total_data %>%select(record_id,Event,Time,drugs___10,drugs___11,drugs___12,drugs___13)
Cocaine_use <- spont_clearance_time_cocaine_use %>% 
    mutate(Cocaine_use=case_when(drugs___10 ==1|
                                 drugs___11==1|
                                 drugs___12==1|
                                 drugs___13==1 ~ "use", TRUE ~"no use"))
cocaine_use_co <-Cocaine_use[c(1:3,8)]
cocaine_use_co$Cocaine_use <- relevel(factor(cocaine_use_co$Cocaine_use),ref="no use")
cocaineuse_KMplot <- ggsurvplot(survfit(Surv(Time,Event)~Cocaine_use,data=cocaine_use_co),
           xlab="Days",
           ylab="Proportion of HCV persistence",
           legend.title="Cocaine Use",
           legend.labs=c("No Cocaine Use","Cocaine Use"),
           pval = TRUE,
           conf.int = TRUE)
ggsave(file="Output/survival analysis_uni_sc/Cocaine_Use_KM.pdf",print(cocaineuse_KMplot),onefile=FALSE)
cox_cocaine_co <-coxph(Surv(Time,Event)~Cocaine_use,data=cocaine_use_co)
sink("Output/survival analysis_uni_sc/cox_cocaine_use.txt")
print(cox_cocaine_co)
sink(file=NULL)
#Proportional hazards assumption
cox_assumption_cocaine <-cox.zph(cox_cocaine_co)
schoenfield_cocaineuse <- ggcoxzph(cox_assumption_cocaine)
ggsave("Output/survival analysis_uni_sc/schoenfield_cocaineuse.pdf", arrangeGrob(grobs=schoenfield_cocaineuse))
########################################
spont_clearance_time_methuse <-total_data %>% select(record_id,Event,Time,drugs___1:drugs___3)
Meth_use <- spont_clearance_time_methuse %>% mutate(Meth_use=case_when(drugs___1 ==1|
                                                                       drugs___2 ==1|
                                                                       drugs___3 ==1 ~"use",TRUE ~ "no use"))
meth_use_co <- Meth_use[c(1:3,7)]
meth_use_co$Meth_use <- relevel(factor(meth_use_co$Meth_use),ref = "no use")
methuse_KMplot <- ggsurvplot(survfit(Surv(Time,Event)~Meth_use,data=meth_use_co),
           xlab="Days",
           ylab="Proportion of HCV persistence",
           legend.title ="Meth Use",
           legend.labs=c("No Meth Use","Meth Use"),
           pval = TRUE,
           conf.int = TRUE)
ggsave(file="Output/survival analysis_uni_sc/Meth_Use_KM.pdf",print(methuse_KMplot),onefile=FALSE)
cox_meth_co <-coxph(Surv(Time,Event)~Meth_use,data=meth_use_co)
sink("Output/survival analysis_uni_sc/cox_meth_use.txt")
print(cox_meth_co)
sink(file=NULL)
#Proportional hazards assumption
cox_assumption_methuse <- cox.zph(cox_meth_co)
schoenfield_methuse <- ggcoxzph(cox_assumption_methuse)
ggsave("Output/survival analysis_uni_sc/schoenfield_methuse.pdf", arrangeGrob(grobs=schoenfield_methuse))
########################################
spont_clearance_time_heroinuse <-total_data %>% select(record_id,Event,Time,drugs___20)
heroin_use <- spont_clearance_time_heroinuse %>% mutate(Heroin_use=case_when(drugs___20==1 ~ "use",TRUE ~"no use"))
heroin_use_co <- heroin_use[c(1:3,5)]
heroin_use_co$Heroin_use <-relevel(factor(heroin_use_co$Heroin_use),ref="no use")
heroinuse_KMplot <- ggsurvplot(survfit(Surv(Time,Event)~Heroin_use,data=heroin_use_co),
           xlab="Days",
           ylab="Proportion of HCV persistence",
           legend.title="Heroin Use",
           legend.labs=c("Heroin No Use", "Heroin Use"),
           pval = TRUE,
           conf.int = TRUE)
ggsave(file="Output/survival analysis_uni_sc/Heroin_Use_KM.pdf",print(heroinuse_KMplot),onefile=FALSE)
summary(survfit(Surv(Time,Event)~Heroin_use,data=heroin_use_co),times=182.5)
summary(survfit(Surv(Time,Event)~Heroin_use,data=heroin_use_co),times=365.25)
cox_heroin_co <-coxph(Surv(Time,Event)~Heroin_use,data=heroin_use_co)
sink("Output/survival analysis_uni_sc/cox_heroin_use.txt")
print(cox_heroin_co)
sink(file=NULL)
#Proportional hazards assumption
cox_assumption_heroin <-cox.zph(cox_heroin_co)
schoenfield_heroinuse <- ggcoxzph(cox_assumption_heroin)
ggsave("Output/survival analysis_uni_sc/schoenfield_heroinuse.pdf", arrangeGrob(grobs=schoenfield_heroinuse))
########################################
#HIV positive patients only
spont_clearance_time_ARVS <- total_data %>% select(record_id,Event,Time,hiv_tx_at_hcv_diagnosis,hiv)
names(spont_clearance_time_ARVS)[4] <-"ARVs"
spont_clearance_time_ARVS_HIV_pos <- spont_clearance_time_ARVS %>% filter(hiv==1)
spont_clearance_time_ARVS_HIV_pos <- spont_clearance_time_ARVS_HIV_pos[,-5]
clean_spont_clearance_time_ARVS_HIV_pos <-na.omit(spont_clearance_time_ARVS_HIV_pos)
clean_spont_clearance_time_ARVS_HIV_pos$ARVs <- relevel(factor(clean_spont_clearance_time_ARVS_HIV_pos$ARVs),ref="0")
ARVs_KMplot <- ggsurvplot(survfit(Surv(Time,Event)~ARVs,data=clean_spont_clearance_time_ARVS_HIV_pos),
           xlab="Days",
           ylab="Proportion of HCV persistence",
           legend.title ="Patients on ARVs",
           legend.labs=c("On ARVs","Not on ARVs"),
           pval = TRUE,
           conf.int = TRUE)
ggsave(file="Output/survival analysis_uni_sc/ARVs_Use_KM.pdf",print(ARVs_KMplot),onefile=FALSE)
cox_ARVS_co <-coxph(Surv(Time,Event)~ARVs,data=clean_spont_clearance_time_ARVS_HIV_pos)
sink("Output/survival analysis_uni_sc/cox_ARVs_use.txt")
print(cox_ARVs_co)
sink(file=NULL)
#Proportional hazards assumption
cox_assumption_ARVS <- cox.zph(cox_ARVS_co)
schoenfield_ARVs <- ggcoxzph(cox_assumption_ARVS)
ggsave("Output/survival analysis_uni_sc/schoenfield_ARVs.pdf", arrangeGrob(grobs=schoenfield_ARVs))
########################################
spont_clearance_time_diabetes <- total_data %>% select(record_id,Time,Event,comorbidities___2)
names(spont_clearance_time_diabetes)[4] <-"Diabetes"
clean_spont_clearance_time_diabetes <- na.omit(spont_clearance_time_diabetes)
diabetes_KMplot <- ggsurvplot(survfit(Surv(Time,Event) ~ Diabetes,data=clean_spont_clearance_time_diabetes),
           xlab="Days",
           ylab="Proportion of HCV persistence",
           legend.title="Diabetes",
           legend.labs=c("No Diabetes","Diabetes"),
           pval = TRUE,
           conf.int = TRUE)
ggsave(file="Output/survival analysis_uni_sc/Diabetes_KM.pdf",print(diabetes_KMplot),onefile=FALSE)
cox_diabetes_co <-coxph(Surv(Time,Event)~Diabetes,data=clean_spont_clearance_time_diabetes)
sink("Output/survival analysis_uni_sc/cox_diabetes.txt")
print(cox_diabetes_co)
sink(file=NULL)
#Proportional hazards assumption
cox_assumption_diabetes <- cox.zph(cox_diabetes_co)
schoenfield_diabetes <- ggcoxzph(cox_assumption_diabetes)
ggsave("Output/survival analysis_uni_sc/schoenfield_diabetes.pdf", arrangeGrob(grobs=schoenfield_diabetes))
########################################
spont_clearance_time_immuno <- total_data %>% select(record_id,Time,Event,immuno)
clean_spont_clearance_time_immuno <-na.omit(spont_clearance_time_immuno)
immuno_given <-clean_spont_clearance_time_immuno %>% mutate(Immunotherapy=case_when(immuno == 1 ~"Yes", TRUE~"No"))
immuno_given <- immuno_given[,-4]
immuno_KMplot <- ggsurvplot(survfit(Surv(Time,Event)~Immunotherapy,data=immuno_given),
           xlab="Days",
           ylab="Proportion of HCV persistence",
           legend.title="Immunotherapy",
           legend.labs=c("Not given immuno.","Given immuno."),
           pval = TRUE,
           conf.int = TRUE)
ggsave(file="Output/survival analysis_uni_sc/Immunotherapy_KM.pdf",print(immuno_KMplot),onefile=FALSE)
cox_immuno_co <-coxph(Surv(Time,Event)~Immunotherapy,data=immuno_given)
sink("Output/survival analysis_uni_sc/cox_immunotherapy.txt")
print(cox_immuno_co)
sink(file=NULL)
#Proportional hazards assumption
cox_assumption_immuno <- cox.zph(cox_immuno_co)
schoenfield_immuno <- ggcoxzph(cox_assumption_immuno)
ggsave("Output/survival analysis_uni_sc/schoenfield_immuno.pdf", arrangeGrob(grobs=schoenfield_immuno))
########################################
spont_clearance_time_cAb_HBV <- total_data %>% select(record_id,Time,Event,cab.factor)
names(spont_clearance_time_cAb_HBV)[4]<-"Core_HBV_Ab"
clean_spont_clearance_time_cAb_HBV <- na.omit(spont_clearance_time_cAb_HBV)
clean_spont_clearance_time_cAb_HBV$Core_HBV_Ab <-relevel(factor(clean_spont_clearance_time_cAb_HBV$Core_HBV_Ab),ref = "Negative")
cAb_HBV_KMplot <- ggsurvplot(survfit(Surv(Time,Event)~Core_HBV_Ab,data=clean_spont_clearance_time_cAb_HBV),
           xlab="Days",
           ylab="Proportion of HCV persistence",
           legend.title="Core Ab HBV",
           legend.labs=c("Negative","Positive"),
           pval = TRUE,
           conf.int = TRUE)
ggsave(file="Output/survival analysis_uni_sc/cAb_HBV_KM.pdf",print(cAb_HBV_KMplot),onefile=FALSE)
cox_cAb_HBV_co <-coxph(Surv(Time,Event)~Core_HBV_Ab,data=clean_spont_clearance_time_cAb_HBV)
sink("Output/survival analysis_uni_sc/cox_cAb_HBV.txt")
print(cox_cAb_HBV_co)
sink(file=NULL)
#Proportional hazards assumption
cox_assumption_cAb_HBV <- cox.zph(cox_cAb_HBV_co)
schoenfield_cAb_HBV <- ggcoxzph(cox_assumption_cAb_HBV)
ggsave("Output/survival analysis_uni_sc/schoenfield_cAb_HBV.pdf", arrangeGrob(grobs=schoenfield_cAb_HBV))
########################################
#Incomplete dataset(PROBLEMATIC)
spont_clearance_time_chronic_HBV <- total_data %>% select(record_id,Time,Event,hbvsag_pcr.factor)
names(spont_clearance_time_chronic_HBV)[4] <-"Chronic_HBV"
clean_spont_clearance_time_chronic_HBV <- na.omit(spont_clearance_time_chronic_HBV)
clean_spont_clearance_time_chronic_HBV$Chronic_HBV<- relevel(factor(clean_spont_clearance_time_chronic_HBV$Chronic_HBV),ref="Negative")
table(clean_spont_clearance_time_chronic_HBV$Event,clean_spont_clearance_time_chronic_HBV$Chronic_HBV)
#Neither KM or Cox Regression are done due
########################################
#HIV positive patients only
spont_clearance_time_CD4 <- total_data %>% select(record_id,Time,Event,cd4_at_hcv_diagnosis,hiv)
names(spont_clearance_time_CD4)[4]<-"CD4_count"
spont_clearance_time_CD4_hiv_pos <- spont_clearance_time_CD4 %>% filter(hiv==1)
spont_clearance_time_CD4_hiv_pos <- spont_clearance_time_CD4_hiv_pos[,-5]
clean_spont_clear_time_cd4 <- na.omit(spont_clearance_time_CD4_hiv_pos)
cox_CD4_co <-coxph(Surv(Time,Event)~CD4_count,data=clean_spont_clear_time_cd4)
sink("Output/survival analysis_uni_sc/cox_CD4.txt")
print(cox_CD4_co)
sink(file=NULL)
#Proportional hazards assumption
cox_assumption_CD4 <- cox.zph(cox_CD4_co)
schoenfield_CD4 <-ggcoxzph(cox_assumption_CD4)
ggsave("Output/survival analysis_uni_sc/schoenfield_CD4.pdf", arrangeGrob(grobs=schoenfield_CD4))
########################################
spont_clearance_time_peakBil <- total_data %>% select(record_id,Time,Event,bil_peak_1)
clean_spont_clearance_time_peakBil <- na.omit(spont_clearance_time_peakBil)
peakBil_binary <- clean_spont_clearance_time_peakBil  %>% mutate(peak_Bil_binary=case_when(bil_peak_1 >=20 ~">= 20",TRUE ~ "< 20"))
peakBil_binary <-peakBil_binary[,-4]
peakBil_KMplot <- ggsurvplot(survfit(Surv(Time,Event)~peak_Bil_binary,data=peakBil_binary),
           xlab="Days",
           ylab="Proportion of HCV persistence",
           legend.title="peak Bilirubin",
           legend.labs=c("peak Bilirubin < 20","peak Bilirubin > 20"),
           pval = TRUE,
           conf.int = TRUE)
ggsave(file="Output/survival analysis_uni_sc/peakBil_KM.pdf",print(peakBil_KMplot),onefile=FALSE)
summary(survfit(Surv(Time,Event)~peak_Bil_binary,data=peakBil_binary),times=182.5)
summary(survfit(Surv(Time,Event)~peak_Bil_binary,data=peakBil_binary),times=365.25)
cox_peakBil_co <-coxph(Surv(Time,Event)~peak_Bil_binary,data=peakBil_binary)
sink("Output/survival analysis_uni_sc/cox_peakBilirubin.txt")
print(cox_peakBil_co)
sink(file=NULL)
#Proportional hazards assumption
cox_assumption_peakBil <-cox.zph(cox_peakBil_co)
schoenfield_peakBil <- ggcoxzph(cox_assumption_peakBil)
ggsave("Output/survival analysis_uni_sc/schoenfield_peakBil.pdf", arrangeGrob(grobs=schoenfield_peakBil))
########################################
spont_clearance_time_b_viral_load <- total_data %>% select(record_id,Event,Time,vl1)
clean_spont_clearance_time_b_viral_load <- na.omit(spont_clearance_time_b_viral_load)
#Binary interpretation (<800,000 > 800,000)
b_viral_load <-clean_spont_clearance_time_b_viral_load  %>% mutate(Viral_load=case_when(vl1 < 800000 ~ "low",TRUE ~ "high"))
b_viral_load_remove_number <- b_viral_load[,-4]
b_viral_load_remove_number$Viral_load <- relevel(factor(b_viral_load_remove_number$Viral_load),ref="low")
b_viral_load_KMplot <- ggsurvplot(survfit(Surv(Time,Event)~Viral_load,data=b_viral_load_remove_number),
           xlab="Days",
           ylab="Proportion of HCV persistence",
           legend.title="Baseline Viral Load",
           legend.labs=c("low(<8x10^5)","high(>8x10^5)"),
           pval = TRUE,
           conf.int = TRUE)
cox_baseline_viral_load_binary <- coxph(Surv(Time,Event)~Viral_load,data=b_viral_load_remove_number)
sink("Output/survival analysis_uni_sc/cox_viral_load.txt")
print(cox_baseline_viral_load_binary)
sink(file=NULL)
#Proportional hazards assumption
cox_assumption_b_viral_load <- cox.zph(cox_baseline_viral_load_binary)
schoenfield_b_viral_load <- ggcoxzph(cox_assumption_b_viral_load)
ggsave("Output/survival analysis_uni_sc/schoenfield_baseline_viral_load.pdf", arrangeGrob(grobs=schoenfield_b_viral_load))
########################################
#ID 66 is documented as being dually infected.
spont_clearance_time_genotype <- total_data %>% select(record_id,Event,Time,clinical_genotype___1,clinical_genotype___3,clinical_genotype___4)
names(spont_clearance_time_genotype)[4:6]<-(c("gt1a","gt3a","gt4d"))
patients_genotyped <- spont_clearance_time_genotype %>%
    gather(Genotype,ID,gt1a:gt4d) %>%
    filter(ID==1)
#Removed ID 66 dual infection cause no evidence in sequencing data
patients_genotyped <- patients_genotyped %>%
    filter(record_id != 66 | Genotype !="gt4d")
patients_genotyped <- patients_genotyped[,-5]
genotype_KMplot <- ggsurvplot(survfit(Surv(Time,Event)~Genotype,data=patients_genotyped),
           xlab="Days",
           ylab="Proportion of HCV persistence",
           legend.title="Infected Genotype",
           pval = TRUE,
           conf.int = TRUE)
cox_genotype_co <-coxph(Surv(Time,Event)~Genotype,data=patients_genotyped)
summary(cox_genotype_co)
#Showing probability after 365 days
summary(survfit(Surv(Time,Event)~Genotype,data=patients_genotyped),times=182.5)
summary(survfit(Surv(Time,Event)~Genotype,data=patients_genotyped),times=365.25)
#Proportional hazards assumption
cox_assumption_genotype <-cox.zph(cox_genotype_co)
ggcoxzph(cox_assumption_genotype)
########################################
spont_clearance_time_ethnicity <- total_data %>% select(record_id,Time,Event,ethnic.factor)
ethnic_grouping <- spont_clearance_time_ethnicity %>% mutate(Ethnic_grouping=case_when(ethnic.factor =="White British"~"White British",
                                                                                         ethnic.factor =="Any other White background" ~"Any other White background",
                                                                                         TRUE ~ "Other ethnic background"))
ethnic_grouping <- ethnic_grouping[,-4]
ethnic_grouping$Ethnic_grouping <- relevel(factor(ethnic_grouping$Ethnic_grouping),ref="White British")
ggsurvplot(survfit(Surv(Time,Event)~Ethnic_grouping,data=ethnic_grouping),
           xlab="Days",
           ylab="Proportion of HCV persistence",
           risk.table = TRUE,
           pval = TRUE,
           conf.int = TRUE)
cox_ethnicity_co <-coxph(Surv(Time,Event)~Ethnic_grouping,data=ethnic_grouping)
summary(cox_ethnicity_co)
#Proportional hazards assumption
cox_assumption_ethnicity <- cox.zph(cox_ethnicity_co)
ggcoxzph(cox_assumption_ethnicity)
########################################
spont_clearance_time_alcohol_excess <- total_data %>%select(record_id,Time,Event,alco_excess)
clean_spont_clearance_time_alcohol_excess <- na.omit(spont_clearance_time_alcohol_excess)
ggsurvplot(survfit(Surv(Time,Event)~alco_excess,data=clean_spont_clearance_time_alcohol_excess),
           xlab="Days",
           ylab="Proportion of HCV persistence",
           risk.table = TRUE,
           pval = TRUE,
           conf.int = TRUE)
cox_alcohol_co <-coxph(Surv(Time,Event)~alco_excess,data=clean_spont_clearance_time_alcohol_excess)
summary(cox_alcohol_co)
#Proportional hazards assumption
cox_assumption_alcohol <-cox.zph(cox_alcohol_co)
ggcoxzph(cox_assumption_alcohol)
########################################
spont_clearance_time_IDU <- total_data %>% select(record_id,Event,Time,risk___1)
names(spont_clearance_time_IDU)[4]<-"PWID"
ggsurvplot(survfit(Surv(Time,Event)~PWID,data=spont_clearance_time_IDU),
           xlab="Days",
           ylab="Proportion of HCV persistence",
           risk.table = TRUE,
           pval = TRUE,
           conf.int = TRUE)
cox_IDU_co <-coxph(Surv(Time,Event)~PWID,data=spont_clearance_time_IDU)
summary(cox_IDU_co)
#Proportional hazards assumption
cox_assumption_IDU <-cox.zph(cox_IDU_co)
ggcoxzph(cox_assumption_IDU)
########################################
spont_clearance_time_MSM <- total_data %>% select(record_id,risk___11,risk___10,risk___2,Event,Time)
spont_clearance_time_MSM_risk <- spont_clearance_time_MSM %>% mutate(MSM=case_when(risk___11 ==1|risk___10 ==1|risk___2==1 ~ "1",TRUE ~ "0"))
spont_clearance_time_MSM_risk <-spont_clearance_time_MSM_risk[c(1,5:7)]
spont_clearance_time_MSM_risk$MSM <- relevel(factor(spont_clearance_time_MSM_risk$MSM),ref="1")
ggsurvplot(survfit(Surv(Time,Event)~MSM,data=spont_clearance_time_MSM_risk),
           xlab="Days",
           ylab="Proportion of HCV persistence",
           risk.table = TRUE,
           pval = TRUE,
           conf.int = TRUE)
summary(survfit(Surv(Time,Event)~MSM,data=spont_clearance_time_MSM_risk),times=182.5)
summary(survfit(Surv(Time,Event)~MSM,data=spont_clearance_time_MSM_risk),times=365.25)
cox_MSM_co <-coxph(Surv(Time,Event)~MSM,data=spont_clearance_time_MSM_risk)
summary(cox_MSM_co)
#Proportional hazards assumption
cox_assumption_MSM <-cox.zph(cox_MSM_co)
ggcoxzph(cox_assumption_MSM)
########################################
spont_clearance_time_weight <-total_data %>% select(record_id,weight,Event,Time)
clean_spont_clearance_time_weight <-na.omit(spont_clearance_time_weight)
cox_weight_co <- coxph(Surv(Time,Event)~weight,data=clean_spont_clearance_time_weight)
summary(cox_weight_co)
#Proportional hazards assumption
cox_assumption_weight <-cox.zph(cox_weight_co)
ggcoxzph(cox_assumption_weight)
########################################
#n=111
spont_clearance_time_il28b <- total_data %>% select(record_id,ifnl4_860.factor,Event,Time)
names(spont_clearance_time_il28b)[2]<-"il28b"
clean_spont_clearance_time_il28b <-na.omit(spont_clearance_time_il28b)
ggsurvplot(survfit(Surv(Time,Event)~il28b,data=clean_spont_clearance_time_il28b),
           xlab="Days",
           ylab="Proportion of HCV persistence",
           risk.table = TRUE,
           pval = TRUE,
           conf.int = TRUE)
cox_il28b_co <-coxph(Surv(Time,Event)~il28b,data=clean_spont_clearance_time_il28b)
summary(cox_il28b_co)
#Proportional hazards assumption
cox_assumption_il28b <-cox.zph(cox_il28b_co)
ggcoxzph(cox_assumption_il28b)
#Combined non-favourable
combined_il28b <- clean_spont_clearance_time_il28b %>%mutate(CC_NonCC=case_when(il28b == 'CC' ~ 'CC',TRUE ~"CT/TT"))
combined_il28b <- combined_il28b[,-2]
combined_il28b$CC_NonCC <- relevel(factor(combined_il28b$CC_NonCC),ref='CT/TT')
ggsurvplot(survfit(Surv(Time,Event)~CC_NonCC,data=combined_il28b),
           xlab="Days",
          ylab="Proportion of HCV persistence",
           risk.table = TRUE,
           pval = TRUE,
         conf.int = TRUE)
cox_il28b_co_combined <-coxph(Surv(Time,Event)~CC_NonCC,data=combined_il28b)
summary(cox_il28b_co_combined)
#Proportional hazards assumption
cox_assumption_il28b_combined <-cox.zph(cox_il28b_co_combined)
ggcoxzph(cox_assumption_il28b_combined)
########################################
#Making dataframe of all above variables including missing data for multivariable analysis(IN PROGRESS)

#Not included are il28b(Still missing considerable data),CD4(HIV+ve),ARVS(HIV+ve)
dataframe_list <- list(spont_clearance_time_gender,spont_clearance_time_age,spont_clearance_time_HIVstatus,
                       peakALT_ClinicalOutcome_binary,all_Drug_use,cocaine_use_co,meth_use_co,heroin_use_co,
                       clean_spont_clearance_time_diabetes,immuno_given,
                       clean_spont_clearance_time_cAb_HBV,clean_spont_clearance_time_chronic_HBV,
                       peakBil_binary,b_viral_load_remove_number,patients_genotyped,ethnic_grouping,
                       spont_clearance_time_alcohol_excess,spont_clearance_time_IDU,
                       spont_clearance_time_MSM_risk,clean_spont_clearance_time_weight,
                       combined_il28b,clean_spont_clearance_time_ARVS_HIV_pos,
                       clean_spont_clear_time_cd4)
#Note ID 66 has two entries due to being dually infected

#This will recursively merge dataframes from the list above
multivariable_df <- Reduce(function(x,y) merge(x,y,all=TRUE),dataframe_list)

#Look at how many missing values there are:
#Not standardised yet for removing negative times currently.Awaiting Emmas response.
sapply(multivariable_df, function(x) sum(is.na(x)))

#Writing a csv file
write.csv(multivariable_df,"multivariable_df.csv",row.names = FALSE)


