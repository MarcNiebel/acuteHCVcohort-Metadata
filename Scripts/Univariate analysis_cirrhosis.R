#Created by Marc Niebel November 2019
#Univariate analysis on acute HCV metadata vs clinical outcome: cirrhosis

#Libraries required
library(dplyr)
library(tidyr)
library(reshape2)
library(lattice)
library(gridExtra)
library(jtools)

###################################################
#Sourcing the dataframe from the data folder
source("./data/AcuteHCV_R_2020-04-22_1120.r", chdir = TRUE)
###################################################
Gender_Cirrhosis <- data %>% select(gender.factor,cirrhosis,record_id)
#Removing chronic patients
Gender_Cirrhosis <- Gender_Cirrhosis %>%
    filter(record_id != 35 & record_id != 213 & record_id != 41) %>%
    select(-record_id)
#Remove any NA rows
clean_Gender_Cirrhosis <- na.omit(Gender_Cirrhosis)
clean_Gender_Cirrhosis$cirrhosis <- factor(clean_Gender_Cirrhosis$cirrhosis)
pdf("Output/univariable_cirrhosis/Gender_cirrhosis.pdf")
plot1 <- ggplot(clean_Gender_Cirrhosis,aes(cirrhosis, fill=gender.factor))+geom_histogram(stat="count")
plot2 <- histogram(~cirrhosis|gender.factor, data=clean_Gender_Cirrhosis)
grid.arrange(plot1,plot2,nrow=1)
dev.off()
#Shows a contigency table
table_gender_cirrhosis <- table(clean_Gender_Cirrhosis)
#Statistics
fisher_g_cirrhosis <- fisher.test(table_gender_cirrhosis)$p.value
fisher_g_cirrhosis <- round(fisher_g_cirrhosis,digits = 4)
#Odds ratio analysis
clean_Gender_Cirrhosis$gender.factor <- relevel(clean_Gender_Cirrhosis$gender.factor,ref="Female")
logit_gender_cirrhosis <- glm(cirrhosis ~ gender.factor, data=clean_Gender_Cirrhosis,family = "binomial")
logit_gender_summ <- summ(logit_gender_cirrhosis,exp=TRUE,digits = 4)
sink("Output/univariable_cirrhosis/Gender_cirrhosis_regression_summary.txt")
print(logit_gender_summ)
sink(file=NULL)
###################################################
Age_Cirrhosis <- data %>% select(age,cirrhosis,record_id)
Age_Cirrhosis <- Age_Cirrhosis %>%
    filter(record_id != 35 & record_id != 213 & record_id != 41) %>%
    select(-record_id)
clean_Age_Cirrhosis <- na.omit(Age_Cirrhosis)
clean_Age_Cirrhosis$cirrhosis <- factor(clean_Age_Cirrhosis$cirrhosis)
pdf("Output/univariable_cirrhosis/Age_cirrhosis.pdf")
plot1 <- ggplot(clean_Age_Cirrhosis,mapping=aes(x=cirrhosis,y=age))+geom_boxplot()+scale_y_continuous()
plot2 <- ggplot(clean_Age_Cirrhosis,aes(age))+geom_bar()
grid.arrange(plot1,plot2,nrow=1)
dev.off()
#Statistics
wilcox_age_cirrhosis <- wilcox.test(age ~ cirrhosis,clean_Age_Cirrhosis)$p.value
wilcox_age_cirrhosis <- round(wilcox_age_cirrhosis,digits = 4)
#Odds ratio analysis
logit_age_cirrhosis <- glm(cirrhosis ~ age, data=clean_Age_Cirrhosis,family="binomial")
logit_age_summ <- summ(logit_age_cirrhosis,exp=TRUE,digits=4)
sink("Output/univariable_cirrhosis/Age_regression_summary.txt")
print(logit_age_summ)
sink(file=NULL)
###################################################
HIV_Cirrhosis <- data %>% select(hiv.factor, cirrhosis,record_id)
HIV_Cirrhosis <- HIV_Cirrhosis %>%
    filter(record_id != 35 & record_id != 213 & record_id != 41) %>%
    select(-record_id)
clean_HIV_Cirrhosis <- na.omit(HIV_Cirrhosis)
clean_HIV_Cirrhosis$cirrhosis <- factor(clean_HIV_Cirrhosis$cirrhosis)
pdf("Output/univariable_cirrhosis/HIV_co-infection_cirrhosis.pdf")
plot1 <- ggplot(clean_HIV_Cirrhosis,aes(cirrhosis, fill=hiv.factor))+geom_histogram(stat="count")
plot2 <- histogram(~cirrhosis | hiv.factor, data=clean_HIV_Cirrhosis)
grid.arrange(plot1,plot2,ncol=1)
dev.off()
table_HIV_cirrhosis <- table(clean_HIV_Cirrhosis)
#Statistics
fisher_HIV_cirrhosis <- fisher.test(table_HIV_cirrhosis)$p.value
fisher_HIV_cirrhosis <- round(fisher_HIV_cirrhosis,digits = 4)
#Odds ratio analysis
clean_HIV_Cirrhosis <- clean_HIV_Cirrhosis %>% mutate(HIV_status=case_when(hiv.factor == "No" ~ "No HIV", TRUE ~ "HIV detected"))
clean_HIV_Cirrhosis <- clean_HIV_Cirrhosis[,c(2,3)]
clean_HIV_Cirrhosis$HIV_status <- relevel(factor(clean_HIV_Cirrhosis$HIV_status),ref="No HIV")
logit_HIV_cirrhosis <- glm(cirrhosis ~ HIV_status, data=clean_HIV_Cirrhosis,family="binomial")
logit_HIV_summ <- summ(logit_HIV_cirrhosis,exp=TRUE,digits = 4)
sink("Output/univariable_cirrhosis/HIV_co-infection_regression_summary.txt")
print(logit_HIV_summ)
sink(file=NULL)
###################################################
#HIV positive patients only
CD4_Cirrhosis <- data %>% select(cd4_at_hcv_diagnosis, cirrhosis,record_id,hiv)
CD4_Cirrhosis <- CD4_Cirrhosis %>%
    filter(record_id != 35 & record_id != 213 & record_id != 41) %>%
    select(-record_id)
CD4_Cirrhosis_hiv_positive <- CD4_Cirrhosis %>%filter(hiv==1)
CD4_Cirrhosis_hiv_positive <- CD4_Cirrhosis_hiv_positive[,-3]
clean_CD4_Cirrhosis <- na.omit(CD4_Cirrhosis_hiv_positive)
clean_CD4_Cirrhosis$cirrhosis <-factor(clean_CD4_Cirrhosis$cirrhosis)
pdf("Output/univariable_cirrhosis/CD4_count_cirrhosis.pdf")
plot1 <- ggplot(clean_CD4_Cirrhosis)+geom_histogram(mapping=aes(x= cd4_at_hcv_diagnosis),binwidth = 100)
plot2 <- ggplot(clean_CD4_Cirrhosis,mapping=aes(x=cirrhosis,y=cd4_at_hcv_diagnosis))+geom_boxplot()+scale_y_continuous()
grid.arrange(plot1,plot2,nrow=1)
dev.off()
#Statistics
ggplot(clean_CD4_Cirrhosis,aes(sample=cd4_at_hcv_diagnosis))+stat_qq()+stat_qq_line()
wilcox_cd4_cirrhosis <- wilcox.test(cd4_at_hcv_diagnosis ~ cirrhosis,data=clean_CD4_Cirrhosis)$p.value
wilcox_cd4_cirrhosis <-round(wilcox_cd4_cirrhosis,digits=4)
#Odds ratio analysis
logit_cd4_cirrhosis <- glm(cirrhosis ~ cd4_at_hcv_diagnosis , data=clean_CD4_Cirrhosis,family="binomial")
logit_cd4_summ <- summ(logit_cd4_cirrhosis,exp=TRUE,digits = 4)
sink("Output/univariable_cirrhosis/CD4_regression_summary.txt")
print(logit_cd4_summ)
sink(file=NULL)
###################################################
peakALT_Cirrhosis <- data %>% select(alt_peak,cirrhosis,record_id)
peakALT_Cirrhosis <- peakALT_Cirrhosis %>% 
    filter(record_id != 35 & record_id != 213 & record_id != 41) %>%
    select(-record_id)
clean_peakALT_Cirrhosis <- na.omit(peakALT_Cirrhosis)
clean_peakALT_Cirrhosis$cirrhosis <- factor(clean_peakALT_Cirrhosis$cirrhosis)
#Using a binary predcitor variable(<1000 and >1000)
peakALT_cirrhosis_binary <- clean_peakALT_Cirrhosis %>% mutate(Binary_peakALT=case_when(alt_peak >=1000 ~ ">1000",TRUE~"<1000"))
peakALT_cirrhosis_binary <- peakALT_cirrhosis_binary[,-1]
pdf("Output/univariable_cirrhosis/peakALT_cirrhosis.pdf")
plot1 <- ggplot(peakALT_cirrhosis_binary,aes(cirrhosis, fill=Binary_peakALT))+geom_histogram(stat="count")
plot2 <- histogram(~cirrhosis | Binary_peakALT, data=peakALT_cirrhosis_binary)
grid.arrange(plot1,plot2,nrow=1)
dev.off()
#Statistics
table_peakALT_cirrhosis_binary <- table(peakALT_cirrhosis_binary)
fisher_peakALT_cirrhosis_binary <- fisher.test(table_peakALT_cirrhosis_binary)$p.value
fisher_peakALT_cirrhosis_binary <- round(fisher_peakALT_cirrhosis_binary, digits = 4)
#Odds ratio analysis
logit_peakALT_cirrhosis_binary <- glm(cirrhosis ~ Binary_peakALT,data=peakALT_cirrhosis_binary,family="binomial")
logit_peakALT_summ <-summ(logit_peakALT_cirrhosis_binary,exp=TRUE,digits = 4)
sink("Output/univariable_cirrhosis/peakALT_regression_summary.txt")
print(logit_peakALT_summ)
sink(file=NULL)
###################################################
cAb_HBV_Cirrhosis <- data %>% select(cab.factor,cirrhosis,record_id)
cAb_HBV_Cirrhosis <- cAb_HBV_Cirrhosis %>%
    filter(record_id != 35 & record_id != 213 & record_id != 41) %>%
    select(-record_id)
clean_cAb_HBV_Cirrhosis <- na.omit(cAb_HBV_Cirrhosis)
clean_cAb_HBV_Cirrhosis$cirrhosis <-factor(clean_cAb_HBV_Cirrhosis$cirrhosis)
pdf("Output/univariable_cirrhosis/cAb_HBV_cirrhosis.pdf")
plot1 <- ggplot(clean_cAb_HBV_Cirrhosis,aes(cirrhosis, fill=cab.factor))+geom_histogram(stat="count")
plot2 <- histogram(~cirrhosis |cab.factor, data=clean_cAb_HBV_Cirrhosis)
grid.arrange(plot1,plot2,nrow=1)
dev.off()
table_cAb_HBV_cirrhosis <- table(clean_cAb_HBV_Cirrhosis)
#Statistics
fisher_cAb_HBV_cirrhosis <- fisher.test(table_cAb_HBV_cirrhosis)$p.value
fisher_cAb_HBV_cirrhosis <- round(fisher_cAb_HBV_cirrhosis,digits = 4)
#Odds ratio analysis
clean_cAb_HBV_Cirrhosis$cab.factor<- relevel(factor(clean_cAb_HBV_Cirrhosis$cab.factor),ref="Negative")
logit_cAb_HBV_cirrhosis <- glm(cirrhosis ~ cab.factor, data=clean_cAb_HBV_Cirrhosis,family="binomial")
logit_cAb_HBV_summ <- summ(logit_cAb_HBV_cirrhosis,exp=TRUE,digits=4)
sink("Output/univariable_cirrhosis/cAb_HBV_regression_summary.txt")
print(logit_cAb_HBV_summ)
sink(file=NULL)
###################################################
chronic_HBV_Cirrhosis <- data %>% select(hbvsag_pcr.factor,cirrhosis,record_id)
chronic_HBV_Cirrhosis <- chronic_HBV_Cirrhosis %>% 
    filter(record_id != 35 & record_id != 213 & record_id != 41) %>%
    select(-record_id)
clean_chronic_HBV_Cirrhosis <- na.omit(chronic_HBV_Cirrhosis)
clean_chronic_HBV_Cirrhosis$cirrhosis <- factor(clean_chronic_HBV_Cirrhosis$cirrhosis)
pdf("Output/univariable_cirrhosis/chronic_HBV_cirrhosis.pdf")
plot1 <- ggplot(clean_chronic_HBV_Cirrhosis,aes(cirrhosis, fill=hbvsag_pcr.factor))+geom_histogram(stat="count")
plot2 <- histogram(~cirrhosis |hbvsag_pcr.factor, data=clean_chronic_HBV_Cirrhosis)
grid.arrange(plot1,plot2,nrow=1)
dev.off()
table_chronic_HBV_cirrhosis <- table(clean_chronic_HBV_Cirrhosis)
#Statistics
fisher_chronic_HBV_cirrhosis <- fisher.test(table_chronic_HBV_cirrhosis)$p.value
fisher_chronic_HBV_cirrhosis <- round(fisher_chronic_HBV_cirrhosis,digits = 4)
#Odds ratio analysis
clean_chronic_HBV_Cirrhosis$hbvsag_pcr.factor <- relevel(factor(clean_chronic_HBV_Cirrhosis$hbvsag_pcr.factor),ref="Negative")
logit_chronic_HBV_cirrhosis <- glm(cirrhosis ~ hbvsag_pcr.factor, data=clean_chronic_HBV_Cirrhosis,family="binomial")
logit_chronic_HBV_summ <- summ(logit_chronic_HBV_cirrhosis,exp=TRUE,digits = 4)
sink("Output/univariable_cirrhosis/chronic_HBV_regression_summary.txt")
print(logit_chronic_HBV_summ)
sink(file=NULL)
###################################################
#Known drug use
Drug_use_Cirrhosis <- data %>% select(drugs___20:drugs___19,cirrhosis,record_id)
Drug_use_Cirrhosis <- Drug_use_Cirrhosis %>%
    filter(record_id != 35 & record_id != 213 & record_id != 41) %>%
    select(-record_id)
clean_Drug_use_Cirrhosis <-na.omit(Drug_use_Cirrhosis)
all_Drug_use_cirrhosis <- clean_Drug_use_Cirrhosis %>% 
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

all_drug_use_cir <- all_Drug_use_cirrhosis[23:24]
all_drug_use_cir$cirrhosis <- factor(all_drug_use_cir$cirrhosis)
pdf("Output/univariable_cirrhosis/Drug_use_cirrhosis.pdf")
plot1 <- ggplot(all_drug_use_cir,aes(cirrhosis, fill=Drug_use))+geom_histogram(stat="count")
plot2 <- histogram(~cirrhosis |Drug_use, data=all_drug_use_cir)
grid.arrange(plot1,plot2,nrow=1)
dev.off()
table_all_drug_use <- table(all_drug_use_cir)
#Statistics
fisher_all_drug_use_cirrhosis <- fisher.test(table_all_drug_use)$p.value
fisher_all_drug_use_cirrhosis <- round(fisher_all_drug_use_cirrhosis,digits = 4)
#Odds ratio analysis
all_drug_use_cir$Drug_use <- relevel(factor(all_drug_use_cir$Drug_use),ref="no use")
logit_drug_use_cirrhosis <- glm(cirrhosis ~Drug_use, data=all_drug_use_cir,family="binomial")
logit_drug_use_summ <- summ(logit_drug_use_cirrhosis,exp=TRUE,digits = 4)
sink("Output/univariable_cirrhosis/Drug_use_regression_summary.txt")
print(logit_drug_use_summ)
sink(file=NULL)
###################################################
Meth_Drug_use_Cirrhosis <- data %>% select(drugs___1:drugs___3,cirrhosis,record_id)
Meth_Drug_use_Cirrhosis <- Meth_Drug_use_Cirrhosis %>%
    filter(record_id != 35 & record_id != 213 & record_id != 41) %>%
    select(-record_id)
clean_Meth_Drug_use_Cirrhosis <- na.omit(Meth_Drug_use_Cirrhosis)
Meth_use_cirrhosis <- clean_Meth_Drug_use_Cirrhosis %>% 
    mutate(use=case_when(drugs___1 ==1|
                         drugs___2 ==1|
                         drugs___3 ==1 ~"Meth use",TRUE ~ "no use"))
meth_use_cir <- Meth_use_cirrhosis[4:5]
meth_use_cir$cirrhosis <- factor(meth_use_cir$cirrhosis)
pdf("Output/univariable_cirrhosis/Meth_use_cirrhosis.pdf")
plot1 <- ggplot(meth_use_cir,aes(cirrhosis, fill=use))+geom_histogram(stat="count")
plot2 <- histogram(~cirrhosis |use, data=meth_use_cir)
grid.arrange(plot1,plot2,nrow=1)
dev.off()
meth_table_cir <- table(meth_use_cir)
#Statistics
fisher_meth_cir <- fisher.test(meth_table_cir)$p.value
fisher_meth_cir <- round(fisher_meth_cir,digits=4)
#Odds ratio analysis
meth_use_cir$use <- relevel(factor(meth_use_cir$use),ref="no use")
logit_meth_use_cirrhosis <- glm(cirrhosis ~ use, data=meth_use_cir, family="binomial")
logit_meth_use_summ <- summ(logit_meth_use_cirrhosis,exp=TRUE,digits=4)
sink("Output/univariable_cirrhosis/Meth_use_regression_summary.txt")
print(logit_meth_use_summ)
sink(file=NULL)
###################################################
Cocaine_use_Cirrhosis <- data %>% select(drugs___10,drugs___11,drugs___12,drugs___13,cirrhosis,record_id)
Cocaine_use_Cirrhosis <- Cocaine_use_Cirrhosis %>%
    filter(record_id != 35 & record_id != 213 & record_id != 41) %>%
    select(-record_id)
clean_Cocaine_use_Cirrhosis <- na.omit(Cocaine_use_Cirrhosis)
Cocaine_use_cirrhosis <-clean_Cocaine_use_Cirrhosis %>% 
    mutate(use=case_when(drugs___10 ==1|
                         drugs___11==1|
                         drugs___12==1|
                         drugs___13==1 ~ "Cocaine use", TRUE ~"no use"))
cocaine_use_cir <- Cocaine_use_cirrhosis[5:6]
cocaine_use_cir$cirrhosis <- factor(cocaine_use_cir$cirrhosis)
pdf("Output/univariable_cirrhosis/Cocaine_use_cirrhosis.pdf")
plot1 <- ggplot(cocaine_use_cir,aes(cirrhosis, fill=use))+geom_histogram(stat="count")
plot2 <- histogram(~cirrhosis|use, data=cocaine_use_cir)
grid.arrange(plot1,plot2,nrow=1)
dev.off()
cocaine_table_cir <-table(cocaine_use_cir)
#Statistics
fisher_cocaine_cir <-fisher.test(cocaine_table_cir)$p.value
fisher_cocaine_cir <- round(fisher_cocaine_cir,digits= 4)
#Odds ratio analysis
cocaine_use_cir$use <- relevel(factor(cocaine_use_cir$use),ref="no use")
logit_cocaine_use_cirrhosis <- glm(cirrhosis ~ use, data=cocaine_use_cir,family="binomial")
logit_cocaine_use_summ <- summ(logit_cocaine_use_cirrhosis,exp=TRUE,digits = 4)
sink("Output/univariable_cirrhosis/Cocaine_use_regression_summary.txt")
print(logit_cocaine_use_summ)
sink(file=NULL)
###################################################
Heroin_use_Cirrhosis <- data %>% select(drugs___20,cirrhosis,record_id)
Heroin_use_Cirrhosis <- Heroin_use_Cirrhosis %>%
    filter(record_id != 35 & record_id != 213 & record_id != 41) %>%
    select(-record_id)
clean_heroin_use_Cirrhosis <- na.omit(Heroin_use_Cirrhosis)
Heroin_use_Cirrhosis <- clean_heroin_use_Cirrhosis %>% mutate(use=case_when(drugs___20==1 ~ "Heroine use",TRUE ~"No use"))
heroin_use_cir <- Heroin_use_Cirrhosis[2:3]
heroin_use_cir$cirrhosis <- factor(heroin_use_cir$cirrhosis)
pdf("Output/univariable_cirrhosis/Heroin_use_cirrhosis.pdf")
plot1 <- ggplot(heroin_use_cir,aes(cirrhosis,fill=use))+geom_histogram(stat="count")
plot2 <- histogram(~cirrhosis|use,data=heroin_use_cir)
grid.arrange(plot1,plot2,nrow=1)
dev.off()
table_heroin_cir <- table(heroin_use_cir)
#Statistics
fisher_heroin_cir <- fisher.test(table_heroin_cir)$p.value
fisher_heroin_cir <- round(fisher_heroin_cir,digits=4)
#Odds ratio analysis
heroin_use_cir$use <- relevel(factor(heroin_use_cir$use),ref="No use")
logit_heroin_use_cirrhosis <- glm(cirrhosis~use, data=heroin_use_cir,family="binomial")
logit_heroin_use_summ <- summ(logit_heroin_use_cirrhosis,exp=TRUE,digits = 4)
sink("Output/univariable_cirrhosis/heroin_use_regression_summary.txt")
print(logit_heroin_use_summ)
sink(file=NULL)
###################################################
b_viral_load_Cirrhosis <- data %>% select(vl1,cirrhosis,record_id)
b_viral_load_Cirrhosis <- b_viral_load_Cirrhosis %>%
    filter(record_id != 35 & record_id != 213) %>%
    select(-record_id)
clean_viral_load_Cirrhosis <- na.omit(b_viral_load_Cirrhosis)
#Using a binary predictor variable(<800,000 and >800,000)
b_viral_load_Cirrhosis_binary <- clean_viral_load_Cirrhosis %>% mutate(Viral_load=case_when(vl1 < 800000 ~ "low",TRUE ~ "high"))
b_viral_load_Cirrhosis_binary <- b_viral_load_Cirrhosis_binary[,-1]
b_viral_load_Cirrhosis_binary$cirrhosis <- factor(b_viral_load_Cirrhosis_binary$cirrhosis)
ggplot(b_viral_load_Cirrhosis_binary,aes(cirrhosis, fill=Viral_load))+geom_histogram(stat="count")
histogram(~cirrhosis| Viral_load, data=b_viral_load_Cirrhosis_binary)
table_b_viral_load_binary <- table(b_viral_load_Cirrhosis_binary)
fisher_b_viral_load_binary_cirrhosis <- fisher.test(table_b_viral_load_binary)$p.value
fisher_b_viral_load_binary_cirrhosis <- round(fisher_b_viral_load_binary_cirrhosis,digits = 4)
#Odds ratio analysis
b_viral_load_Cirrhosis_binary$Viral_load <- relevel(factor(b_viral_load_Cirrhosis_binary$Viral_load),ref="low")
logit_baselineviralload_cirrhosis <- glm(cirrhosis ~ Viral_load, data=b_viral_load_Cirrhosis_binary,family="binomial")
logit_b_viral_load_summ <- summ(logit_baselineviralload_cirrhosis,exp=TRUE,digits = 4)
###################################################
#HIV patients only!!
ARVS_Cirrhosis <- data %>% select(hiv_tx_at_hcv_diagnosis.factor, cirrhosis,record_id,hiv)
ARVS_Cirrhosis <- ARVS_Cirrhosis %>%
    filter(record_id != 35 & record_id != 213) %>%
    select(-record_id)
ARVS_Cirrhosis_hiv_pos <-ARVS_Cirrhosis %>% filter(hiv==1)
ARVS_Cirrhosis_hiv_pos <- ARVS_Cirrhosis_hiv_pos[,-3]
clean_ARVS_Cirrhosis <- na.omit(ARVS_Cirrhosis_hiv_pos)
clean_ARVS_Cirrhosis$cirrhosis <-factor(clean_ARVS_Cirrhosis$cirrhosis)
ggplot(clean_ARVS_Cirrhosis,aes(cirrhosis, fill=hiv_tx_at_hcv_diagnosis.factor))+geom_histogram(stat="count")
histogram(~cirrhosis |hiv_tx_at_hcv_diagnosis.factor , data=clean_ARVS_Cirrhosis)
table_ARVS_cirrhosis <- table (clean_ARVS_Cirrhosis)
fisher_arvs_cirrhosis <- fisher.test(table_ARVS_cirrhosis)$p.value
fisher_arvs_cirrhosis <- round(fisher_arvs_cirrhosis,digits = 4)
#Odds ratio analysis
clean_ARVS_Cirrhosis$hiv_tx_at_hcv_diagnosis.factor <- relevel(clean_ARVS_Cirrhosis$hiv_tx_at_hcv_diagnosis.factor,ref = "No")
logit_ARVS_cirrhosis <- glm(cirrhosis ~hiv_tx_at_hcv_diagnosis.factor , data=clean_ARVS_Cirrhosis,family="binomial")
logit_ARVS_summ <- summ(logit_ARVS_cirrhosis,exp=TRUE,digits = 4)
###################################################
immunotherapy_Cirrhosis <- data %>% select(immuno.factor,cirrhosis,record_id)
immunotherapy_Cirrhosis <- immunotherapy_Cirrhosis %>%
    filter(record_id != 35 & record_id != 213) %>%
    select(-record_id)
clean_immunotherapy_Cirrhosis <-na.omit(immunotherapy_Cirrhosis)
clean_immunotherapy_Cirrhosis$cirrhosis <- factor(clean_immunotherapy_Cirrhosis$cirrhosis)
ggplot(clean_immunotherapy_Cirrhosis,aes(cirrhosis, fill=immuno.factor))+geom_histogram(stat="count")
histogram(~cirrhosis |immuno.factor, data=clean_immunotherapy_Cirrhosis)
table_immuno_cirrhosis <- table(clean_immunotherapy_Cirrhosis)
fisher_immuno_cirrhosis <- fisher.test(table_immuno_cirrhosis)$p.value
fisher_immuno_cirrhosis <- round(fisher_immuno_cirrhosis,digits = 4)
#Odds ratio analysis
clean_immunotherapy_Cirrhosis$immuno.factor <- relevel(clean_immunotherapy_Cirrhosis$immuno.factor,ref="No")
logit_immuno_cirrhosis <- glm(cirrhosis ~ immuno.factor, data=clean_immunotherapy_Cirrhosis,family="binomial")
logit_immuno_summ <- summ(logit_immuno_cirrhosis,exp=TRUE,digits = 4)
#Caution only 11 patients have immunotherapy!!
###################################################
peakBil_Cirrhosis <- data %>% select(bil_peak_1,cirrhosis,record_id)
peakBil_Cirrhosis <- peakBil_Cirrhosis %>%
    filter(record_id != 35 & record_id != 213) %>%
    select(-record_id)
clean_peakBil_Cirrhosis <- na.omit(peakBil_Cirrhosis)
clean_peakBil_Cirrhosis$cirrhosis <- factor(clean_peakBil_Cirrhosis$cirrhosis)
#Binary predictor outcome(>=20 and <20)
clean_peakBil_Cirrhosis_binary <- clean_peakBil_Cirrhosis %>% mutate(peak_Bil_binary=case_when(bil_peak_1 >=20 ~">= 20",TRUE ~ "< 20"))
clean_peakBil_Cirrhosis_binary <- clean_peakBil_Cirrhosis_binary[,-1]
ggplot(clean_peakBil_Cirrhosis_binary,aes(cirrhosis,fill=peak_Bil_binary))+geom_histogram(stat="count")
histogram(~cirrhosis|peak_Bil_binary,data=clean_peakBil_Cirrhosis_binary)
table_peakBil_binary_Cirrhosis <-table(clean_peakBil_Cirrhosis_binary)
fisher_peakBil_binary_cirrhosis <-fisher.test(table_peakBil_binary_Cirrhosis)$p.value
fisher_peakBil_binary_cirrhosis <- round(fisher_peakBil_binary_cirrhosis,digits = 4)
#Odds ratio analysis
logit_peakBil_cirrhosis <- glm(cirrhosis ~ peak_Bil_binary, data=clean_peakBil_Cirrhosis_binary,family="binomial")
logit_peakBil_summ <- summ(logit_peakBil_cirrhosis,exp=TRUE,digits = 4)
###################################################
IDU_Cirrhosis <- data %>% select(risk___1,cirrhosis,record_id)
IDU_Cirrhosis <- IDU_Cirrhosis %>%
    filter(record_id != 35 & record_id != 213) %>%
    select(-record_id)
clean_IDU_Cirrhosis <-na.omit(IDU_Cirrhosis)
clean_IDU_Cirrhosis$cirrhosis <- factor(clean_IDU_Cirrhosis$cirrhosis)
IDU_risk_cirrhosis <- clean_IDU_Cirrhosis %>% mutate(Risk=case_when(risk___1==1 ~"IDU",TRUE ~ "no IDU"))
IDU_risk_cir <- IDU_risk_cirrhosis[2:3]
ggplot(IDU_risk_cir,aes(cirrhosis,fill=Risk))+geom_histogram(stat="count")
histogram(~cirrhosis|Risk,data=IDU_risk_cir)
table_IDU_cir <- table(IDU_risk_cir)
fisher_IDU_cirrhosis <-fisher.test(table_IDU_cir)$p.value
fisher_IDU_cirrhosis <- round(fisher_IDU_cirrhosis,digits = 4)
#Odds ratio analysis
IDU_risk_cir$Risk <-relevel(factor(IDU_risk_cir$Risk),ref="no IDU")
logit_IDU_cirrhosis <- glm(cirrhosis ~ Risk, data=IDU_risk_cir,family="binomial")
logit_IDU_summ <- summ(logit_IDU_cirrhosis,exp=TRUE,digits = 4)
###################################################
MSM_Cirrhosis <- data %>% select(risk___11,risk___10,risk___2,cirrhosis,record_id)
MSM_Cirrhosis <- MSM_Cirrhosis %>%
    filter(record_id != 35 & record_id != 213) %>%
    select(-record_id)
clean_MSM_Cirrhosis <- na.omit(MSM_Cirrhosis)
MSM_risk_cirrhosis <- clean_MSM_Cirrhosis %>% mutate(Risk=case_when(risk___11 ==1|risk___10 ==1|risk___2==1 ~ "MSM",TRUE ~ "not MSM"))
MSM_risk_cir <-MSM_risk_cirrhosis[4:5]
MSM_risk_cir$cirrhosis <-factor(MSM_risk_cir$cirrhosis)
ggplot(MSM_risk_cir,aes(cirrhosis,fill=Risk))+geom_histogram(stat="count")
histogram(~cirrhosis|Risk,data=MSM_risk_cir)
table_MSM_cir <- table(MSM_risk_cir)
fisher_MSM_cir <- fisher.test(table_MSM_cir)$p.value
fisher_MSM_cir <- round(fisher_MSM_cir,digits = 4)
#Odds ratio analysis
MSM_risk_cir$Risk <- relevel(factor(MSM_risk_cir$Risk),ref="not MSM")
logit_MSM_cirrhosis <- glm(cirrhosis ~ Risk, data=MSM_risk_cir,family="binomial")
logit_MSM_risk_summ <-summ(logit_MSM_cirrhosis,exp=TRUE, digits = 4)
###################################################
Ethnicity_Cirrhosis <- data %>% select(ethnic.factor,cirrhosis,record_id)
Ethnicity_Cirrhosis <- Ethnicity_Cirrhosis %>%
    filter(record_id != 35 & record_id != 213) %>%
    select(-record_id)
clean_ethnicity_cirrhosis <- na.omit(Ethnicity_Cirrhosis)
#Grouped
group_ethnicity_cirrhosis <-clean_ethnicity_cirrhosis %>% mutate(Ethnic_grouping=case_when(ethnic.factor =="White British"~"White British",
                                                                                                 ethnic.factor =="Any other White background" ~"Any other White background",
                                                                                                 TRUE ~ "Other ethnic background"))
group_ethnicity_cirrhosis <- group_ethnicity_cirrhosis[,-1]
group_ethnicity_cirrhosis$cirrhosis <- factor(group_ethnicity_cirrhosis$cirrhosis)
ggplot(group_ethnicity_cirrhosis,aes(cirrhosis, fill=Ethnic_grouping))+geom_histogram(stat="count")
#Alternative
histogram(~cirrhosis | Ethnic_grouping, data=group_ethnicity_cirrhosis)
#contigency table
table_group_ethnicity_cirrhosis <- table(group_ethnicity_cirrhosis)
#Combining further cause 0s in other ethnic background 
group_ethnicity_cirrhosis_combined_ethnicity <- group_ethnicity_cirrhosis %>%
    mutate(New_ethnic_grouping=case_when(Ethnic_grouping == "Other ethnic background"|Ethnic_grouping =="Any other White background" ~"Other",
                                        TRUE ~ "White British"))
group_ethnicity_cirrhosis_combined_ethnicity <- group_ethnicity_cirrhosis_combined_ethnicity[,-2]
ggplot(group_ethnicity_cirrhosis_combined_ethnicity,aes(cirrhosis, fill=New_ethnic_grouping))+geom_histogram(stat="count")
histogram(~cirrhosis | New_ethnic_grouping, data=group_ethnicity_cirrhosis_combined_ethnicity)
table_group_ethnicity_cirrhosis_combined_ethnicity <-table(group_ethnicity_cirrhosis_combined_ethnicity)
fisher_group_eth_cirrhosis <-fisher.test(table_group_ethnicity_cirrhosis_combined_ethnicity)$p.value
fisher_group_eth_cirrhosis <- round(fisher_group_eth_cirrhosis,digits = 4)
#Odds ratio analysis
group_ethnicity_cirrhosis_combined_ethnicity$New_ethnic_grouping <- relevel(factor(group_ethnicity_cirrhosis_combined_ethnicity$New_ethnic_grouping),ref="White British")
logit_ethnicity_cirrhosis <- glm(cirrhosis ~ New_ethnic_grouping, data=group_ethnicity_cirrhosis_combined_ethnicity,family="binomial")
logit_ethnicity_summ <- summ(logit_ethnicity_cirrhosis,exp=TRUE,digits = 4)
##################################################
#Only 1 gt2a so not included(Note that ID 66 is dually infected)
Genotype_Cirrhosis <- data %>% select(clinical_genotype___1,clinical_genotype___3,clinical_genotype___4,cirrhosis,record_id)
Genotype_Cirrhosis <- Genotype_Cirrhosis %>% 
    filter(record_id != 35 & record_id != 213) %>%
    select(-record_id)
Genotype_Cirrhosis$cirrhosis <- factor(Genotype_Cirrhosis$cirrhosis)
clean_Genotype_Cirrhosis <- na.omit(Genotype_Cirrhosis)
patients_genotyped_cirrhosis <- clean_Genotype_Cirrhosis %>%
    gather(Genotype,id,clinical_genotype___1:clinical_genotype___4) %>%
    filter(id==1)
patients_genotyped_cirrhosis <- patients_genotyped_cirrhosis[1:2]
#Looking at the distribution
ggplot(patients_genotyped_cirrhosis,aes(cirrhosis, fill=Genotype))+geom_histogram(stat="count")
histogram(~cirrhosis|Genotype, data=patients_genotyped_cirrhosis)
table_genotype_cirrhosis <- table(patients_genotyped_cirrhosis)
fisher_genotype_cirrhosis <- fisher.test(table_genotype_cirrhosis)$p.value
fisher_genotype_cirrhosis <- round(fisher_genotype_cirrhosis,digits=4)
#Odds ratio analysis
mylogit_genotype_cirrhosis <-glm(cirrhosis ~ Genotype, data=patients_genotyped_cirrhosis,family="binomial")
logit_genotype_summ <- summ(mylogit_genotype_cirrhosis,exp=TRUE,digits = 4)
##################################################
Diabetic_Cirrhosis <- data %>% select(comorbidities___2,cirrhosis,record_id)
Diabetic_Cirrhosis <- Diabetic_Cirrhosis %>%
    filter(record_id != 35 & record_id != 213) %>%
    select(-record_id)
clean_diabetic_Cirrhosis <- na.omit(Diabetic_Cirrhosis)
clean_diabetic_Cirrhosis <- clean_diabetic_Cirrhosis %>% mutate(Diabetes=case_when(comorbidities___2== 0~"No diabetes",TRUE~"Diabetes"))
clean_diabetic_Cirrhosis <- clean_diabetic_Cirrhosis[,-1]
clean_diabetic_Cirrhosis$cirrhosis <- factor(clean_diabetic_Cirrhosis$cirrhosis)
ggplot(clean_diabetic_Cirrhosis,aes(cirrhosis, fill=Diabetes))+geom_histogram(stat="count")
histogram(~cirrhosis |Diabetes , data=clean_diabetic_Cirrhosis)
table_diabetic_cirrhosis <- table (clean_diabetic_Cirrhosis)
fisher_diabetic_cirrhosis <- fisher.test(table_diabetic_cirrhosis)$p.value
fisher_diabetic_cirrhosis <- round(fisher_diabetic_cirrhosis,digits=4)
#Odds ratio analysis
clean_diabetic_Cirrhosis$Diabetes <- relevel(factor(clean_diabetic_Cirrhosis$Diabetes),ref="No diabetes")
logit_diabetes_cirrhosis <- glm(cirrhosis ~ Diabetes, data=clean_diabetic_Cirrhosis,family="binomial")
logit_diabetes_summ <- summ(logit_diabetes_cirrhosis,exp = TRUE,digits = 4)
##################################################
Alcohol_Excess_Cirrhosis <- data %>% select(alco_excess,cirrhosis,record_id)
Alcohol_Excess_Cirrhosis <- Alcohol_Excess_Cirrhosis %>%
    filter(record_id != 35 & record_id != 213) %>%
    select(-record_id)
clean_Alcohol_Excess_Cirrhosis <-na.omit(Alcohol_Excess_Cirrhosis)
clean_Alcohol_Excess_Cirrhosis$cirrhosis <- factor(clean_Alcohol_Excess_Cirrhosis$cirrhosis)
clean_Alcohol_Excess_Cirrhosis$alco_excess <- factor(clean_Alcohol_Excess_Cirrhosis$alco_excess)
ggplot(clean_Alcohol_Excess_Cirrhosis,aes(cirrhosis, fill=alco_excess))+geom_histogram(stat="count")
histogram(~cirrhosis|alco_excess,data=clean_Alcohol_Excess_Cirrhosis)
table_alcohol_excess_cirrhosis <- table(clean_Alcohol_Excess_Cirrhosis)
fisher_alcohol_cirrhosis <- fisher.test(table_alcohol_excess_cirrhosis)$p.value
fisher_alcohol_cirrhosis <- round(fisher_alcohol_cirrhosis,digits = 4)
#Odds ratio analysis
clean_Alcohol_Excess_Cirrhosis$alco_excess <- relevel(clean_Alcohol_Excess_Cirrhosis$alco_excess,ref="0")
logit_alcohol_cirrhosis <- glm(cirrhosis ~ alco_excess, data=clean_Alcohol_Excess_Cirrhosis,family="binomial")
logit_alcohol_excess_summ <- summ(logit_alcohol_cirrhosis,exp = TRUE,digits=4)
##################################################
weight_cirrhosis <- data %>% select(weight,cirrhosis,record_id)
weight_cirrhosis <- weight_cirrhosis %>% 
    filter(record_id != 35 & record_id != 213) %>%
    select(-record_id)
clean_weight_cirrhosis <- na.omit(weight_cirrhosis)
clean_weight_cirrhosis$cirrhosis <- factor(clean_weight_cirrhosis$cirrhosis)
ggplot(clean_weight_cirrhosis,mapping=aes(x=cirrhosis,y=weight))+geom_boxplot()+scale_y_continuous()
ggplot(clean_weight_cirrhosis)+geom_histogram(mapping=aes(x= weight),binwidth = 10)
ggplot(clean_weight_cirrhosis,aes(sample=weight))+stat_qq()+stat_qq_line()
wilcox_weight_cirrhosis <- wilcox.test(weight ~ cirrhosis,clean_weight_cirrhosis)$p.value
wilcox_weight_cirrhosis <- round(wilcox_weight_cirrhosis,digits=4)
#Odds ratio analysis
logit_weight_cirrhosis <- glm(cirrhosis ~ weight, data=clean_weight_cirrhosis,family="binomial")
logit_weight_summ <- summ(logit_weight_cirrhosis,exp=TRUE)
##################################################
#INCOMPLETE(VERY LOW NUMBERS CURRENTLY IN REGARDS TO CIRRHOTIC PATIENTS BEING GENOTYPED)
#il28b_Cirrhosis <- data %>% select(ifnl4_860.factor,cirrhosis,record_id)
#il28b_Cirrhosis <- il28b_Cirrhosis %>% 
#    filter(record_id != 35 & record_id != 213) %>%
#    select(-record_id)
#clean_il28b_Cirrhosis <-na.omit(il28b_Cirrhosis)
#clean_il28b_Cirrhosis$cirrhosis <- factor(clean_il28b_Cirrhosis$cirrhosis)
#ggplot(clean_il28b_Cirrhosis,aes(cirrhosis,fill=ifnl4_860.factor))+geom_histogram(stat="count")
#histogram(~cirrhosis|ifnl4_860.factor, data=clean_il28b_Cirrhosis)
#table_il28b_cirrhosis <- table(clean_il28b_Cirrhosis)
#fisher_il28b_cirrhosis <-fisher.test(table_il28b_cirrhosis)$p.value
#Odds ratio analysis
#logit_il28b_cirrhosis <- glm(cirrhosis ~ ifnl4_860.factor, data=,family="binomial")
#logit_il28b_cirrhosis_summ <-summ(logit_il28b_cirrhosis,exp=TRUE)
#Combine CT and TT due to low numbers of TT and not wanting to use CC as reference
#combined_CT_TT_levels_cirrhosis <- clean_il28b_Cirrhosis %>%
#    mutate(Genotype_binary=case_when(ifnl4_860.factor == "CT"| ifnl4_860.factor == "TT" ~"non-CC",
#                                     TRUE ~ "CC"))
#combined_CT_TT_levels_cirrhosis <- combined_CT_TT_levels_cirrhosis[,-1]
#table_combined_il28b <-table(combined_CT_TT_levels_cirrhosis)
#fisher_combined_il28b_co <- fisher.test(table_combined_il28b)$p.value
#Odds ratio analysis
#combined_CT_TT_levels_cirrhosis$Genotype_binary <- relevel(factor(combined_CT_TT_levels_cirrhosis$Genotype_binary),ref = "non-CC")
#logit_combined_il28b_cirrhosis <- glm(cirrhosis ~ Genotype_binary,data=combined_CT_TT_levels_cirrhosis,family = "binomial")
#logit_combined_il28b_summ <- summ(logit_combined_il28b_cirrhosis,exp=TRUE)
###################################################
#IL28b missing
stats_table <- data.frame(Variable=c("Gender","Age","Ethnicity","HIV status","CD4(HIV +ve patients only)","ARVS(HIV+ve patients only)","peak ALT(>1000)","cAb HBV","chronic HBV",
                                     "Drug use","Cocaine use","Methamphetamine use","Heroin use", "baseline viral load(>800,000 IU/ml)",
                                     "peak Bilirubin(>20)","Risk factor:PWID", "Risk factor:MSM","Genotype","Diabetes",
                                     "Alcohol excess","Weight","Immunotherapy"),
                          p_value=c(fisher_g_cirrhosis,wilcox_age_cirrhosis,fisher_group_eth_cirrhosis,fisher_HIV_cirrhosis,wilcox_cd4_cirrhosis,
                                    fisher_arvs_cirrhosis,fisher_peakALT_cirrhosis_binary,fisher_cAb_HBV_cirrhosis,fisher_chronic_HBV_cirrhosis,
                                    fisher_all_drug_use_cirrhosis,fisher_cocaine_cir,fisher_meth_cir,fisher_heroin_cir,
                                    fisher_b_viral_load_binary_cirrhosis,fisher_peakBil_binary_cirrhosis,fisher_IDU_cirrhosis,fisher_MSM_cir,
                                    fisher_genotype_cirrhosis,fisher_diabetic_cirrhosis,fisher_alcohol_cirrhosis,wilcox_weight_cirrhosis,
                                    fisher_immuno_cirrhosis))
pdf("stats_table_cirrhosis.pdf")
grid.table(stats_table)
dev.off()

