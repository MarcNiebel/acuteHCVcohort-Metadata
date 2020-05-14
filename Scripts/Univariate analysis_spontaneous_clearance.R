#Created by Marc Niebel November 2019
#Univariate analysis on acute HCV metadata for clinical outcomes
#progression/spontaneous clearance
#Libraries required
library(dplyr)
library(tidyr)
library(lattice)
library(gridExtra)
library(jtools)

###################################################
#Sourcing the dataframe from the data folder
source("./data/AcuteHCV_R_2020-04-22_1120.r", chdir = TRUE)
###################################################
Gender_ClinicalOutcome <- data %>% select(gender.factor,sc,record_id)
#Removing chronic patients
Gender_ClinicalOutcome <- Gender_ClinicalOutcome %>% 
    filter(record_id != 35 & record_id != 213 & record_id != 41) %>%
    select(-record_id)
Gender_ClinicalOutcome$sc <- factor(Gender_ClinicalOutcome$sc)
pdf("Output/univariable_sc/Gender_sc.pdf")
plot1 <- ggplot(Gender_ClinicalOutcome,aes(sc, fill=gender.factor))+geom_histogram(stat="count")
#Alternative
plot2 <- histogram(~sc | gender.factor, data=Gender_ClinicalOutcome)
grid.arrange(plot1,plot2,nrow=1)
dev.off()
#Shows a contigency table
table_gender <- table(Gender_ClinicalOutcome)
#Statistics
fisher_go_co <-fisher.test(table_gender)$p.value
fisher_go_co <-round(fisher_go_co,digits = 4)
#Odds ratio analysis
Gender_ClinicalOutcome$gender.factor <- relevel(Gender_ClinicalOutcome$gender.factor,ref="Male")
logit_gender_clinicaloutcome <-glm(sc ~ gender.factor, data=Gender_ClinicalOutcome,family="binomial")
logit_gender_co <- summ(logit_gender_clinicaloutcome,exp=TRUE,digits = 4)
sink("Output/univariable_sc/Gender_regression_summary.txt")
print(logit_gender_co)
sink(file=NULL)
###################################################
Age_ClinicalOutcome <- data %>% select(age,sc,record_id)
Age_ClinicalOutcome <- Age_ClinicalOutcome %>% 
    filter(record_id != 35 & record_id != 213 & record_id != 41) %>%
    select(-record_id)
Age_ClinicalOutcome$sc <- factor(Age_ClinicalOutcome$sc)
pdf("Output/univariable_sc/Age_sc.pdf")
plot1 <- ggplot(Age_ClinicalOutcome,mapping=aes(x=sc,y=age))+geom_boxplot()+scale_y_continuous()
plot2 <- ggplot(Age_ClinicalOutcome,aes(age))+geom_bar()
grid.arrange(plot1,plot2,nrow=1)
dev.off()
#Statistics
shapiro.test(Age_ClinicalOutcome$age)
wilcox_age_co <- wilcox.test(age ~ sc,Age_ClinicalOutcome)$p.value
wilcox_age_co <-round(wilcox_age_co,digits = 4)
#Odds ratio analysis
logit_age_clinicaloutcome <- glm(sc ~ age, data=Age_ClinicalOutcome,family="binomial")
logit_age_co <-summ(logit_age_clinicaloutcome,exp=TRUE,digits=4)
sink("Output/univariable_sc/Age_regression_summary.txt")
print(logit_age_co)
sink(file=NULL)
###################################################
Ethnicity_ClinicalOutcome <- data %>% select(ethnic.factor, sc,record_id)
Ethnicity_ClinicalOutcome <-Ethnicity_ClinicalOutcome %>% 
    filter(record_id != 35 & record_id != 213 & record_id != 41) %>%
    select(-record_id)
#Numbers are low so grouping into White british, any other white background, other ethnic background
group_ethnicity_clinicaloutcome <- Ethnicity_ClinicalOutcome %>% 
    mutate(Ethnic_grouping=case_when(ethnic.factor =="White British"~"White British",
                                     ethnic.factor =="Any other White background" ~"Any other White background",
                                     TRUE ~ "Other ethnic background"))
group_ethnicity_clinicaloutcome <- group_ethnicity_clinicaloutcome[,-1]
group_ethnicity_clinicaloutcome$sc <-factor(group_ethnicity_clinicaloutcome$sc)
pdf("Output/univariable_sc/Ethnicity_sc.pdf",paper="USr")
plot1 <- ggplot(group_ethnicity_clinicaloutcome,aes(sc, fill=Ethnic_grouping))+geom_histogram(stat="count")
#Alternative
plot2 <- histogram(~sc | Ethnic_grouping, data=group_ethnicity_clinicaloutcome)
grid.arrange(plot1,plot2,ncol=1)
dev.off()
table_group_ethnicity <- table(group_ethnicity_clinicaloutcome)
#Statistics
fisher_eth_group_co <- fisher.test(table_group_ethnicity)$p.value
fisher_eth_group_co <- round(fisher_eth_group_co,digits = 4)
#Odds ratio analysis
group_ethnicity_clinicaloutcome$Ethnic_grouping <-factor(group_ethnicity_clinicaloutcome$Ethnic_grouping)
group_ethnicity_clinicaloutcome$Ethnic_grouping <-relevel(group_ethnicity_clinicaloutcome$Ethnic_grouping,ref="White British")
logit_ethnicity_clinicaloutcome <-glm(sc ~ Ethnic_grouping, data=group_ethnicity_clinicaloutcome,family="binomial")
logit_ethnicty_co <- summ(logit_ethnicity_clinicaloutcome,exp=TRUE,digits=4)
sink("Output/univariable_sc/Ethnicity_regression_summary.txt")
print(logit_ethnicty_co)
sink(file=NULL)
###################################################
HIV_ClinicalOutcome <- data %>% select(hiv.factor,sc,record_id)
HIV_ClinicalOutcome <- HIV_ClinicalOutcome %>% 
    filter(record_id != 35 & record_id != 213 & record_id != 41) %>%
    select(-record_id)
clean_HIV_clinicaloutcome <- na.omit(HIV_ClinicalOutcome)
clean_HIV_clinicaloutcome$sc <- factor(clean_HIV_clinicaloutcome$sc)
#Looking at distribution
pdf("Output/univariable_sc/HIV_co-infection_sc.pdf")
plot1 <- ggplot(clean_HIV_clinicaloutcome,aes(sc, fill=hiv.factor))+geom_histogram(stat="count")
plot2 <- histogram(~sc | hiv.factor, data=clean_HIV_clinicaloutcome,col=c("red","seagreen"))
grid.arrange(plot1,plot2,nrow=1)
dev.off()
table_HIV <- table(clean_HIV_clinicaloutcome)
#Statistics
fisher_hiv_co <- fisher.test(table_HIV)$p.value
fisher_hiv_co <- round(fisher_hiv_co,digits = 4)
#Odds ratio analysis
clean_HIV_clinicaloutcome <- clean_HIV_clinicaloutcome %>% 
    mutate(HIV_status=case_when(hiv.factor == "No" ~ "No HIV", TRUE ~ "HIV detected"))
clean_HIV_clinicaloutcome <- clean_HIV_clinicaloutcome[,c(2:3)]
clean_HIV_clinicaloutcome$HIV_status <-factor(clean_HIV_clinicaloutcome$HIV_status)
logit_HIV_clinicaloutcome <- glm(sc ~ HIV_status,data=clean_HIV_clinicaloutcome,family="binomial")
logit_HIV_co <- summ(logit_HIV_clinicaloutcome,exp=TRUE, digits=4)
sink("Output/univariable_sc/HIV_co-infection_regression_summary.txt")
print(logit_HIV_co)
sink(file=NULL)
###################################################
#HIV positive patients only
CD4_ClinicalOutcome <- data %>% select(cd4_at_hcv_diagnosis, sc,hiv,record_id)
CD4_ClinicalOutcome <- CD4_ClinicalOutcome %>% 
    filter(record_id != 35 & record_id != 213 & record_id != 41) %>%
    select(-record_id)
CD4_ClinicalOutcome_hiv_positive <- CD4_ClinicalOutcome %>%filter(hiv==1)
CD4_ClinicalOutcome_hiv_positive <- CD4_ClinicalOutcome_hiv_positive[,-3]
clean_CD4_ClinicalOutcome <- na.omit(CD4_ClinicalOutcome_hiv_positive)
clean_CD4_ClinicalOutcome$sc <- factor(clean_CD4_ClinicalOutcome$sc)
pdf("Output/univariable_sc/CD4_count_sc.pdf")
plot1 <- ggplot(clean_CD4_ClinicalOutcome,mapping=aes(x=sc,y=cd4_at_hcv_diagnosis))+geom_boxplot()+scale_y_continuous()
plot2 <- ggplot(clean_CD4_ClinicalOutcome)+geom_histogram(mapping=aes(x=cd4_at_hcv_diagnosis),binwidth = 100)
grid.arrange(plot1,plot2,nrow=1)
dev.off()
#Statistics
ggplot(clean_CD4_ClinicalOutcome,aes(sample=cd4_at_hcv_diagnosis))+stat_qq()+stat_qq_line()
shapiro.test(clean_CD4_ClinicalOutcome$cd4_at_hcv_diagnosis)
wilcox_cd4_co <- wilcox.test(cd4_at_hcv_diagnosis ~ sc,data=clean_CD4_ClinicalOutcome)$p.value
wilcox_cd4_co <- round(wilcox_cd4_co,digits = 4)
#Odds ratio analysis
logit_CD4_clinicaloutcome <- glm(sc ~ cd4_at_hcv_diagnosis, data=clean_CD4_ClinicalOutcome,family="binomial")
logit_cd4_co <- summ(logit_CD4_clinicaloutcome,exp=TRUE,digits = 4)
sink("Output/univariable_sc/CD4_regression_summary.txt")
print(logit_cd4_co)
sink(file=NULL)
###################################################
peakALT_ClinicalOutcome <- data %>% select(alt_peak,sc,record_id)
peakALT_ClinicalOutcome <- peakALT_ClinicalOutcome %>% 
    filter(record_id != 35 & record_id != 213 & record_id != 41) %>%
    select(-record_id)
clean_peakALT_ClinicalOutcome <- na.omit(peakALT_ClinicalOutcome)
clean_peakALT_ClinicalOutcome$sc <- factor(clean_peakALT_ClinicalOutcome$sc)
#Using a binary predictor variable(<1000 and >1000)
peakALT_ClinicalOutcome_binary <- clean_peakALT_ClinicalOutcome %>% 
    mutate(Binary_peakALT=case_when(alt_peak >=1000 ~ ">1000",TRUE~"<1000"))
peakALT_ClinicalOutcome_binary <- peakALT_ClinicalOutcome_binary[,-1]
pdf("Output/univariable_sc/peakALT_sc.pdf")
plot1 <- ggplot(peakALT_ClinicalOutcome_binary,aes(sc, fill=Binary_peakALT))+geom_histogram(stat="count")
plot2 <- histogram(~sc | Binary_peakALT, data=peakALT_ClinicalOutcome_binary)
grid.arrange(plot1,plot2,nrow=1)
dev.off()
table_peakALT_binary <- table(peakALT_ClinicalOutcome_binary)
#Statistics
fisher_peakALT_binary_co <- fisher.test(table_peakALT_binary)$p.value
fisher_peakALT_binary_co <- round(fisher_peakALT_binary_co,digits = 4)
#Odds ratio analysis
logit_peakALT_clincaloutcome_binary <- glm(sc ~ Binary_peakALT,data=peakALT_ClinicalOutcome_binary,family="binomial")
logit_peakALT_binary_co <-summ(logit_peakALT_clincaloutcome_binary,exp=TRUE,digits = 4)
sink("Output/univariable_sc/peakALT_regression_summary.txt")
print(logit_peakALT_binary_co)
sink(file=NULL)
###################################################
cAb_HBV_ClinicalOutcome <- data %>% select(cab.factor,sc,record_id)
cAb_HBV_ClinicalOutcome <- cAb_HBV_ClinicalOutcome %>% 
    filter(record_id != 35 & record_id != 213 & record_id != 41) %>%
    select(-record_id)
clean_cAb_HBV_ClinicalOutcome <- na.omit(cAb_HBV_ClinicalOutcome)
clean_cAb_HBV_ClinicalOutcome$sc <-factor(clean_cAb_HBV_ClinicalOutcome$sc)
pdf("Output/univariable_sc/cAb_HBV_sc.pdf")
plot1 <- ggplot(clean_cAb_HBV_ClinicalOutcome,aes(sc, fill=cab.factor))+geom_histogram(stat="count")
plot2 <- histogram(~sc |cab.factor, data=clean_cAb_HBV_ClinicalOutcome)
grid.arrange(plot1,plot2,nrow=1)
dev.off()
table_cAb_HBV <- table(clean_cAb_HBV_ClinicalOutcome)
#Statistics
fisher_cAb_HBV_co <- fisher.test(table_cAb_HBV)$p.value
fisher_cAb_HBV_co <- round(fisher_cAb_HBV_co,digits = 4)
#Odds ratio analysis
clean_cAb_HBV_ClinicalOutcome$cab.factor <-relevel(clean_cAb_HBV_ClinicalOutcome$cab.factor,ref="Negative")
logit_cAb_HBV_clinicaloutcome <- glm(sc ~ cab.factor, data=clean_cAb_HBV_ClinicalOutcome,family="binomial")
logit_cAb_HBV_co <- summ(logit_cAb_HBV_clinicaloutcome,exp=TRUE,digits = 4)
sink("Output/univariable_sc/cAb_HBV_regression_summary.txt")
print(logit_cAb_HBV_co)
sink(file=NULL)
###################################################
chronic_HBV_ClinicalOutcome <- data %>% select(hbvsag_pcr.factor,sc,record_id)
chronic_HBV_ClinicalOutcome <- chronic_HBV_ClinicalOutcome %>% 
    filter(record_id != 35 & record_id != 213 & record_id != 41) %>%
    select(-record_id)
clean_chronic_HBV_ClinicalOutcome <- na.omit(chronic_HBV_ClinicalOutcome)
clean_chronic_HBV_ClinicalOutcome$sc <- factor(clean_chronic_HBV_ClinicalOutcome$sc)
pdf("Output/univariable_sc/chronic_HBV_sc.pdf")
plot1 <- ggplot(clean_chronic_HBV_ClinicalOutcome,aes(sc, fill=hbvsag_pcr.factor))+geom_histogram(stat="count")
plot2 <- histogram(~sc |hbvsag_pcr.factor, data=clean_chronic_HBV_ClinicalOutcome)
grid.arrange(plot1,plot2,nrow=1)
dev.off()
table_chronic_HBV <- table(clean_chronic_HBV_ClinicalOutcome)
#No statistics due to 0 value cells
###################################################
Drug_use_ClinicalOutcome <- data %>% select(drugs___20:drugs___19,sc,record_id)
Drug_use_ClinicalOutcome <- Drug_use_ClinicalOutcome %>% 
    filter(record_id != 35 & record_id != 213 & record_id != 41) %>%
    select(-record_id)
all_Drug_use <- Drug_use_ClinicalOutcome %>% 
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

all_drug_use_co <- all_Drug_use[23:24]
all_drug_use_co$sc <-factor(all_drug_use_co$sc)
pdf("Output/univariable_sc/Drug_use_sc.pdf")
plot1 <- ggplot(all_drug_use_co,aes(sc, fill=Drug_use))+geom_histogram(stat="count")
plot2 <- histogram(~sc |Drug_use, data=all_drug_use_co)
grid.arrange(plot1,plot2,nrow=1)
dev.off()
table_all_drug_use <- table(all_drug_use_co)
#Statistics
fisher_all_drug_use_co <- fisher.test(table_all_drug_use)$p.value
fisher_all_drug_use_co <-round(fisher_all_drug_use_co,digits = 4)
#Odds ratio analysis
all_drug_use_co$Drug_use <- relevel(factor(all_drug_use_co$Drug_use),ref="no use")
logit_druguse_clinicaloutcome <- glm(sc ~Drug_use, data=all_drug_use_co,family="binomial")
logit_drug_use_co <- summ(logit_druguse_clinicaloutcome,exp=TRUE,digits=4)
sink("Output/univariable_sc/Drug_use_regression_summary.txt")
print(logit_drug_use_co)
sink(file=NULL)
###################################################
#Number of drugs taken by patients
Drugs_used_record_id <- data %>% select(record_id,sc,drugs___20:drugs___19)
Drugs_taken <- Drugs_used_record_id %>%
    gather(Drugs,Recorded,drugs___20:drugs___19) %>%
    filter(!Recorded==0) %>%
    group_by(record_id) %>%
    arrange(record_id) %>%
    mutate(count=n())
#Patients ranging from 1 drug use to 10 drugs!! Cannot disentangle them easily
###################################################
Cocaine_use_ClinicalOutcome <- data %>% select(drugs___10,drugs___11,drugs___12,drugs___13,sc,record_id)
Cocaine_use_ClinicalOutcome <- Cocaine_use_ClinicalOutcome %>% 
    filter(record_id != 35 & record_id != 213 & record_id != 41) %>%
    select(-record_id)
Cocaine_use <-Cocaine_use_ClinicalOutcome %>% 
    mutate(use=case_when(drugs___10 ==1|
                         drugs___11==1|
                         drugs___12==1|
                         drugs___13==1 ~ "Cocaine use", TRUE ~"no use"))
cocaine_use_co <- Cocaine_use[5:6]
cocaine_use_co$sc <- factor(cocaine_use_co$sc)
pdf("Output/univariable_sc/Cocaine_use_sc.pdf")
plot1 <- ggplot(cocaine_use_co,aes(sc, fill=use))+geom_histogram(stat="count")
plot2 <- histogram(~sc|use, data=cocaine_use_co)
grid.arrange(plot1,plot2,nrow=1)
dev.off()
cocaine_table <-table(cocaine_use_co)
#Statistics
fisher_cocaine_use_co <- fisher.test(cocaine_table)$p.value
fisher_cocaine_use_co <- round(fisher_cocaine_use_co,digits=4)
#Odds ratio analysis
cocaine_use_co$use <- relevel(factor(cocaine_use_co$use), ref="no use")
logit_cocaineuse_clinicaloutcome <- glm(sc ~ use, data=cocaine_use_co,family="binomial")
logit_cocaine_co <-summ(logit_cocaineuse_clinicaloutcome,exp=TRUE,digits = 4)
sink("Output/univariable_sc/Cocaine_use_regression_summary.txt")
print(logit_cocaine_co)
sink(file=NULL)
###################################################
Meth_Drug_use_ClinicalOutcome <- data %>% select(drugs___1:drugs___3,sc,record_id)
Meth_Drug_use_ClinicalOutcome <- Meth_Drug_use_ClinicalOutcome %>% 
    filter(record_id != 35 & record_id != 213 & record_id != 41) %>%
    select(-record_id)
Meth_use <- Meth_Drug_use_ClinicalOutcome %>% 
    mutate(use=case_when(drugs___1 ==1|
                         drugs___2 ==1|
                         drugs___3 ==1 ~"Meth use",TRUE ~ "no use"))
meth_use_co <- Meth_use[4:5]
meth_use_co$sc <-factor(meth_use_co$sc)
pdf("Output/univariable_sc/Meth_use_sc.pdf")
plot1 <- ggplot(meth_use_co,aes(sc, fill=use))+geom_histogram(stat="count")
plot2 <- histogram(~sc |use, data=meth_use_co)
grid.arrange(plot1,plot2,nrow=1)
dev.off()
meth_table <- table(meth_use_co)
#Statistics
fisher_meth_co <- fisher.test(meth_table)$p.value
fisher_meth_co <- round(fisher_meth_co,digits =4)
#Odds ratio analysis
meth_use_co$use <- relevel(factor(meth_use_co$use),ref="no use")
logit_meth_use_clinicaloutcome <- glm(sc ~ use, data=meth_use_co, family="binomial")
logit_meth_use_co <-summ(logit_meth_use_clinicaloutcome,exp=TRUE,digits = 4)
sink("Output/univariable_sc/Meth_use_regression_summary.txt")
print(logit_meth_use_co)
sink(file=NULL)
###################################################
Heroin_use_ClinicalOutcome <- data %>% select(drugs___20,sc,record_id)
Heroin_use_ClinicalOutcome <- Heroin_use_ClinicalOutcome %>% 
    filter(record_id != 35 & record_id != 213 & record_id != 41 ) %>%
    select(-record_id)
Heroin_use <- Heroin_use_ClinicalOutcome %>% 
    mutate(use=case_when(drugs___20==1 ~ "Heroin use",TRUE ~"No use"))
heroin_use_co <- Heroin_use[,-1]
heroin_use_co$sc <- factor(heroin_use_co$sc)
pdf("Output/univariable_sc/Heroin_use_sc.pdf")
plot1 <- ggplot(heroin_use_co,aes(sc,fill=use))+geom_histogram(stat="count")
plot2 <- histogram(~sc|use,data=heroin_use_co)
grid.arrange(plot1,plot2,nrow=1)
dev.off()
table_heroin <- table(heroin_use_co)
#Statistics
fisher_heroin_co <- fisher.test(table_heroin)$p.value
fisher_heroin_co <- round(fisher_heroin_co,digits=4)
#Odds ratio analysis
heroin_use_co$use <-factor(heroin_use_co$use)
heroin_use_co$use <-relevel(factor(heroin_use_co$use), ref="No use")
logit_heroinuse_clinicaloutcome <- glm(sc~use, data=heroin_use_co,family="binomial")
logit_heroin_co <- summ(logit_heroinuse_clinicaloutcome,exp=TRUE,digits = 4)
sink("Output/univariable_sc/heroin_use_regression_summary.txt")
print(logit_heroin_co)
sink(file=NULL)
###################################################
b_viral_load_ClinicalOutcome <- data %>% select(vl1,sc,record_id)
b_viral_load_ClinicalOutcome <- b_viral_load_ClinicalOutcome %>% 
    filter(record_id != 35 & record_id != 213 & record_id != 41) %>%
    select(-record_id)
clean_viral_load_ClinicalOutcome <- na.omit(b_viral_load_ClinicalOutcome)
#Using a binary predictor variable(<800,000 and >800,000)
b_viral_load_ClinicalOutcome_binary <- clean_viral_load_ClinicalOutcome %>% 
    mutate(Viral_load=case_when(vl1 < 800000 ~ "low",TRUE ~ "high"))
b_viral_load_ClinicalOutcome_binary <- b_viral_load_ClinicalOutcome_binary[,-1]
b_viral_load_ClinicalOutcome_binary$sc <- factor(b_viral_load_ClinicalOutcome_binary$sc)
pdf("Output/univariable_sc/basline_viral_load_sc.pdf")
plot1 <- ggplot(b_viral_load_ClinicalOutcome_binary,aes(sc, fill=Viral_load))+geom_histogram(stat="count")
plot2 <- histogram(~sc | Viral_load, data=b_viral_load_ClinicalOutcome_binary)
grid.arrange(plot1,plot2,nrow=1)
dev.off()
table_b_viral_load_binary <- table(b_viral_load_ClinicalOutcome_binary)
#Statistics
fisher_b_viral_load_binary_co <- fisher.test(table_b_viral_load_binary)$p.value
fisher_b_viral_load_binary_co <- round(fisher_b_viral_load_binary_co,digits = 4)
#odds ratio analysis
b_viral_load_ClinicalOutcome_binary$Viral_load <- relevel(factor(b_viral_load_ClinicalOutcome_binary$Viral_load),ref="low")
logit_baselineviralload_clinicaloutcome <- glm(sc ~ Viral_load, data=b_viral_load_ClinicalOutcome_binary,family="binomial")
logit_b_viral_load_co <-summ(logit_baselineviralload_clinicaloutcome,exp=TRUE,digits=4)
sink("Output/univariable_sc/baseline_viral_load_regresssion_summary.txt")
print(logit_b_viral_load_co)
sink(file=NULL)
###################################################
#HIV_positive patients only
ARVS_ClinicalOutcome <- data %>% select(hiv_tx_at_hcv_diagnosis.factor, sc,hiv,record_id)
ARVS_ClinicalOutcome <- ARVS_ClinicalOutcome %>% 
    filter(record_id != 35 & record_id != 213 & record_id != 41) %>%
    select(-record_id)
ARVS_ClinicalOutcome_hiv_pos <- ARVS_ClinicalOutcome %>% filter(hiv==1)
ARVS_ClinicalOutcome_hiv_pos <- ARVS_ClinicalOutcome_hiv_pos[,-3]
clean_ARVS_ClincalOutcome <- na.omit(ARVS_ClinicalOutcome_hiv_pos)
clean_ARVS_ClincalOutcome$sc <-factor(clean_ARVS_ClincalOutcome$sc)
pdf("Output/univariable_sc/ARVS_sc.pdf")
plot1 <- ggplot(clean_ARVS_ClincalOutcome,aes(sc, fill=hiv_tx_at_hcv_diagnosis.factor))+geom_histogram(stat="count")
plot2 <- histogram(~sc |hiv_tx_at_hcv_diagnosis.factor , data=clean_ARVS_ClincalOutcome)
grid.arrange(plot1,plot2,nrow=1)
dev.off()
table_ARVS_co <- table(clean_ARVS_ClincalOutcome)
#Statistics
fisher_arvs_co <- fisher.test(table_ARVS_co)$p.value
fisher_arvs_co <- round(fisher_arvs_co,digits=4)
#Odds ratio analysis
clean_ARVS_ClincalOutcome$hiv_tx_at_hcv_diagnosis.factor <- relevel(clean_ARVS_ClincalOutcome$hiv_tx_at_hcv_diagnosis.factor,ref="No")
logit_ARVS_clinicaloutcome <- glm(sc ~ hiv_tx_at_hcv_diagnosis.factor, data=clean_ARVS_ClincalOutcome,family="binomial")
logit_ARVS_co <-summ(logit_ARVS_clinicaloutcome,exp=TRUE,digits=4)
sink("Output/univariable_sc/ARVS_sc.txt")
print(logit_ARVS_co)
sink(file=NULL)
###################################################
immunotherapy_ClinicalOutcome <- data %>% select(immuno.factor,sc,record_id)
immunotherapy_ClinicalOutcome <- immunotherapy_ClinicalOutcome %>% 
    filter(record_id != 35 & record_id != 213 & record_id != 41) %>%
    select(-record_id)
clean_immunotherapy_ClinicalOutcome <-na.omit(immunotherapy_ClinicalOutcome)
clean_immunotherapy_ClinicalOutcome$sc <- factor(clean_immunotherapy_ClinicalOutcome$sc)
pdf("Output/univariable_sc/immunotherapy_sc.pdf")
plot1 <- ggplot(clean_immunotherapy_ClinicalOutcome,aes(sc, fill=immuno.factor))+geom_histogram(stat="count")
plot2 <- histogram(~sc |immuno.factor, data=clean_immunotherapy_ClinicalOutcome)
grid.arrange(plot1,plot2,nrow=1)
dev.off()
table_immuno <- table(clean_immunotherapy_ClinicalOutcome)
#Statistic
fisher_immuno_co <- fisher.test(table_immuno)$p.value
fisher_immuno_co <- round(fisher_immuno_co,digits=4)
#Odds ratio analysis
clean_immunotherapy_ClinicalOutcome$immuno.factor <- relevel(clean_immunotherapy_ClinicalOutcome$immuno.factor,ref="No")
logit_immuno_clinicaloutcome <- glm(sc ~ immuno.factor, data=clean_immunotherapy_ClinicalOutcome,family="binomial")
logit_immuno_co <- summ(logit_immuno_clinicaloutcome,exp=TRUE,digits=4)
sink("Output/univariable_sc/Immunotherapy_regression_summary.txt")
print(logit_immuno_co)
sink(file=NULL)
###################################################
peakBil_ClinicalOutcome <- data %>% select(bil_peak_1,sc,record_id)
peakBil_ClinicalOutcome <- peakBil_ClinicalOutcome %>% 
    filter(record_id != 35 & record_id != 213 & record_id != 41) %>%
    select(-record_id)
clean_peakBil_ClinicalOutcome <- na.omit(peakBil_ClinicalOutcome)
clean_peakBil_ClinicalOutcome$sc <- factor(clean_peakBil_ClinicalOutcome$sc)
#Binary predictor outcome(>=20 and <20)
clean_peakBil_ClinicalOutcome_binary <- clean_peakBil_ClinicalOutcome %>% 
    mutate(peak_Bil_binary=case_when(bil_peak_1 >=20 ~">= 20",TRUE ~ "< 20"))
clean_peakBil_ClinicalOutcome_binary <- clean_peakBil_ClinicalOutcome_binary[,-1]
pdf("Output/univariable_sc/peakBilirubin_sc.pdf")
plot1 <- ggplot(clean_peakBil_ClinicalOutcome_binary,aes(sc,fill=peak_Bil_binary))+geom_histogram(stat="count")
plot2 <- histogram(~sc|peak_Bil_binary,data=clean_peakBil_ClinicalOutcome_binary)
grid.arrange(plot1,plot2,nrow=1)
dev.off()
table_peakBil_binary_clinicaloutcome <-table(clean_peakBil_ClinicalOutcome_binary)
#Statistics
fisher_peakBil_binary_co <-fisher.test(table_peakBil_binary_clinicaloutcome)$p.value
fisher_peakBil_binary_co <- round(fisher_peakBil_binary_co,digits=4)
#Odds ratio analysis
logit_peakBil_clinicaloutcome_binary <- glm(sc ~ peak_Bil_binary, data=clean_peakBil_ClinicalOutcome_binary,family="binomial")
logit_peakBil_co <-summ(logit_peakBil_clinicaloutcome_binary,exp=TRUE,digits=4)
sink("Output/univariable_sc/peakBilirubin_regression_summary.txt")
print(logit_peakBil_co)
sink(file=NULL)
###################################################
IDU_ClinicalOutcome <- data %>% select(risk___1,sc,record_id)
IDU_ClinicalOutcome <- IDU_ClinicalOutcome %>% 
    filter(record_id != 35 & record_id != 213) %>%
    select(-record_id)
clean_IDU_ClinicalOutcome <-na.omit(IDU_ClinicalOutcome)
IDU_risk <- clean_IDU_ClinicalOutcome %>% 
    mutate(Risk=case_when(risk___1==1 ~"IDU",TRUE ~ "no IDU"))
IDU_risk_co <- IDU_risk[2:3]
IDU_risk_co$sc <- factor(IDU_risk_co$sc)
pdf("Output/univariable_sc/IDU_risk_sc.pdf")
plot1 <- ggplot(IDU_risk_co,aes(sc,fill=Risk))+geom_histogram(stat="count")
plot2 <- histogram(~sc|Risk,data=IDU_risk_co)
grid.arrange(plot1,plot2,nrow=1)
dev.off()
table_IDU <- table(IDU_risk_co)
#Statistics
fisher_IDU_co <-fisher.test(table_IDU)$p.value
fisher_IDU_co <- round(fisher_IDU_co,digits = 4)
#Odds ratio analysis
IDU_risk_co$Risk <- relevel(factor(IDU_risk_co$Risk),ref="no IDU")
logit_IDU_clinicaloutcome <- glm(sc ~ Risk, data=IDU_risk_co,family="binomial")
logit_IDU_co <-summ(logit_IDU_clinicaloutcome,exp=TRUE,digits = 4)
sink("Output/univariable_sc/IDU_risk_regression.txt")
print(logit_IDU_co)
sink(file=NULL)
###################################################
MSM_ClincalOutcome <- data %>% select(risk___11,risk___10,risk___2,sc,record_id)
MSM_ClincalOutcome <- MSM_ClincalOutcome %>% 
    filter(record_id != 35 & record_id != 213) %>%
    select(-record_id)
clean_MSM_ClinicalOutcome <- na.omit(MSM_ClincalOutcome)
MSM_risk <- clean_MSM_ClinicalOutcome %>% mutate(Risk=case_when(risk___11 ==1|risk___10 ==1|risk___2==1 ~ "MSM",TRUE ~ "not MSM"))
MSM_risk_co <-MSM_risk[4:5]
MSM_risk_co$sc <- factor(MSM_risk_co$sc)
pdf("Output/univariable_sc/MSM_sc.pdf")
plot1 <- ggplot(MSM_risk_co,aes(sc,fill=Risk))+geom_histogram(stat="count")
plot2 <- histogram(~sc|Risk,data=MSM_risk_co)
grid.arrange(plot1,plot2,nrow=1)
dev.off()
table_MSM <- table(MSM_risk_co)
#Statistics
fisher_MSM_co <- fisher.test(table_MSM)$p.value
fisher_MSM_co <- round(fisher_MSM_co,digits = 4)
#Odds ratio analysis
logit_MSM_clinicaloutcome <- glm(sc ~ Risk, data=MSM_risk_co,family="binomial")
logit_MSM_co <- summ(logit_MSM_clinicaloutcome,exp=TRUE,digits = 4)
sink("Output/univariable_sc/MSM_risk_regression_summary.txt")
print(logit_MSM_co)
sink(file=NULL)
###################################################
#Selecting genotype and clinical outcome(excluding gt2 cause only 1 patient
#ID 66 is documented as being dually infected.
Genotype_ClinicalOutcome <- data %>% select(clinical_genotype___1,clinical_genotype___3,clinical_genotype___4,sc,record_id)
Genotype_ClinicalOutcome <- Genotype_ClinicalOutcome %>% 
    filter(record_id != 35 & record_id != 213 & record_id != 41) %>%
    select(-record_id)
Genotype_ClinicalOutcome$sc <- factor(Genotype_ClinicalOutcome$sc)
clean_Genotype_clinicalOutcome <- na.omit(Genotype_ClinicalOutcome)
patients_genotyped <-clean_Genotype_clinicalOutcome %>%
    gather(Genotype,id,clinical_genotype___1:clinical_genotype___4) %>%
    filter(id==1)
patients_genotyped <- patients_genotyped[1:2]
pdf("Output/univariable_sc/infected_genotype_sc.pdf")
plot1 <- ggplot(patients_genotyped,aes(sc, fill=Genotype))+geom_histogram(stat="count")
plot2 <- histogram(~sc|Genotype, data=patients_genotyped)
grid.arrange(plot1,plot2,ncol=1)
dev.off()
table_genotype <- table(patients_genotyped)
#Statistics
fisher_genotype_co <-fisher.test(table_genotype)$p.value
fisher_genotype_co <- round(fisher_genotype_co,digits=4)
#Odds ratio analysis
logit_genotype_clinical_outcome <-glm(sc ~ Genotype, data=patients_genotyped,family="binomial")
logit_genotype_co <-summ(logit_genotype_clinical_outcome,exp=TRUE,digits = 4)
sink("Output/univariable_sc/Infected_genotype_regression_summary.txt")
print(logit_genotype_co)
sink(file=NULL)
###################################################
diabetic_ClinicalOutcome <- data %>% select(comorbidities___2,sc,record_id)
diabetic_ClinicalOutcome <- diabetic_ClinicalOutcome %>% 
    filter(record_id != 35 & record_id != 213 & record_id != 41) %>%
    select(-record_id)
clean_diabetic_ClinicalOutcome <-na.omit(diabetic_ClinicalOutcome)
clean_diabetic_ClinicalOutcome <- clean_diabetic_ClinicalOutcome %>% mutate(Diabetes=case_when(comorbidities___2 == 0 ~ "No diabetes",TRUE~"diabetes"))
clean_diabetic_ClinicalOutcome <- clean_diabetic_ClinicalOutcome[,-1]
clean_diabetic_ClinicalOutcome$sc <- factor(clean_diabetic_ClinicalOutcome$sc)
pdf("Output/univariable_sc/Diabetes_sc.pdf")
plot1 <- ggplot(clean_diabetic_ClinicalOutcome,aes(sc, fill=Diabetes))+geom_histogram(stat="count")
plot2 <- histogram(~sc |Diabetes, data=clean_diabetic_ClinicalOutcome)
grid.arrange(plot1,plot2,nrow=1)
dev.off()
table_diabetic <- table(clean_diabetic_ClinicalOutcome)
#Statistics
fisher_diabetic_co <- fisher.test(table_diabetic)$p.value
fisher_diabetic_co <-round(fisher_diabetic_co,digits=4)
#Odds ratio analysis
clean_diabetic_ClinicalOutcome$Diabetes <- relevel(factor(clean_diabetic_ClinicalOutcome$Diabetes),ref="No diabetes")
logit_diabetic_clinicaloutcome <- glm(sc ~ Diabetes, data=clean_diabetic_ClinicalOutcome,family="binomial")
logit_diabetic_co <-summ(logit_diabetic_clinicaloutcome,exp=TRUE,digits = 4)
sink("Output/univariable_sc/Diabetes_regression_summary.txt")
print(logit_diabetic_co)
sink(file=NULL)
###################################################
il28b_ClinicalOutcome <- data %>% select(ifnl4_860.factor,sc,record_id)
il28b_ClinicalOutcome <-il28b_ClinicalOutcome %>% 
filter(record_id != 35 & record_id != 213 & record_id != 41) %>%
   select(-record_id)
clean_il28b_ClinicalOutcome <-na.omit(il28b_ClinicalOutcome)
clean_il28b_ClinicalOutcome$sc <- factor(clean_il28b_ClinicalOutcome$sc)
pdf("Output/univariable_sc/IL28b_sc.pdf")
plot1 <- ggplot(clean_il28b_ClinicalOutcome,aes(sc,fill=ifnl4_860.factor))+geom_histogram(stat="count")
plot2 <- histogram(~sc|ifnl4_860.factor, data=clean_il28b_ClinicalOutcome)
grid.arrange(plot1,plot2,nrow=1)
dev.off()
table_il28b <- table(clean_il28b_ClinicalOutcome)
#Statistic
fisher_il28b_co <-fisher.test(table_il28b)$p.value
fisher_il28b_co <- round(fisher_il28b_co,digits=4)
#Odds ratio analysis
logit_il28b_clinicaloutcome <- glm(sc ~ ifnl4_860.factor, data=clean_il28b_ClinicalOutcome,family="binomial")
logit_il28b_summ <-summ(logit_il28b_clinicaloutcome,exp=TRUE)
#Combine CT and TT due to low numbers of TT and not wanting to use CC(associated with SC) as reference
combined_CT_TT_levels <- clean_il28b_ClinicalOutcome %>% 
    mutate(Genotype_binary=case_when(ifnl4_860.factor == "CT"| 
                                     ifnl4_860.factor == "TT" ~"non-CC",
                                     TRUE ~ "CC"))
combined_CT_TT_levels <- combined_CT_TT_levels[,-1]
table_combined_ifnl4_860 <-table(combined_CT_TT_levels)
#Statistics
fisher_combined_il28b_co <- fisher.test(table_combined_ifnl4_860)$p.value
fisher_combined_il28b_co <- round(fisher_combined_il28b_co,digits=4)
#Odds ratio analysis
combined_CT_TT_levels$Genotype_binary <- relevel(factor(combined_CT_TT_levels$Genotype_binary),ref = "non-CC")
logit_combined_il28b_clinicaloutcome <- glm(sc ~ Genotype_binary,data=combined_CT_TT_levels,family = "binomial")
logit_combined_il28b_co <- summ(logit_combined_il28b_clinicaloutcome,exp=TRUE)
sink("Output/univariable_sc/IL28B_regression_summary.txt")
print(logit_combined_il28b_co)
sink(file=NULL)
###################################################
alcohol_excess_Clinical_Outcome <- data %>% select(alco_excess,sc,record_id) 
alcohol_excess_Clinical_Outcome <- alcohol_excess_Clinical_Outcome %>% 
    filter(record_id != 35 & record_id != 213 & record_id != 41 ) %>%
    select(-record_id)
clean_alcohol_excess_Clinical_Outcome <- na.omit(alcohol_excess_Clinical_Outcome)
clean_alcohol_excess_Clinical_Outcome <- clean_alcohol_excess_Clinical_Outcome %>% 
    mutate(Alcohol_excess=case_when(alco_excess=="1" ~ "alcohol excess",TRUE ~ "no excess"))
clean_alcohol_excess_Clinical_Outcome <- clean_alcohol_excess_Clinical_Outcome[,-1]
clean_alcohol_excess_Clinical_Outcome$sc <- factor(clean_alcohol_excess_Clinical_Outcome$sc)
pdf("Output/univariable_sc/alcohol_excess_sc.pdf")
plot1 <- ggplot(clean_alcohol_excess_Clinical_Outcome,aes(sc,fill=Alcohol_excess))+geom_histogram(stat="count")
plot2 <- histogram(~sc|Alcohol_excess,data=clean_alcohol_excess_Clinical_Outcome)
grid.arrange(plot1,plot2,nrow=1)
dev.off()
table_alcohol_excess_clinicaloutcome <-table(clean_alcohol_excess_Clinical_Outcome)
fisher_alcohol_excess_co <- fisher.test(table_alcohol_excess_clinicaloutcome)$p.value
fisher_alcohol_excess_co <-round(fisher_alcohol_excess_co,digits=4)
#Odds ratio analysis
clean_alcohol_excess_Clinical_Outcome$Alcohol_excess <-relevel(factor(clean_alcohol_excess_Clinical_Outcome$Alcohol_excess),ref="no excess")
logit_alcohol_excess_clinicaloutcome <- glm(sc ~ Alcohol_excess,data=clean_alcohol_excess_Clinical_Outcome,family="binomial")
logit_alcohol_excess_co <- summ(logit_alcohol_excess_clinicaloutcome,exp=TRUE,digits = 4)
sink("Output/univariable_sc/Alcohol_excess_regression_summary.txt")
print(logit_alcohol_excess_co)
sink(file=NULL)
###################################################
weight_Clinical_Outcome <- data %>% select(weight,sc,record_id)
weight_Clinical_Outcome <- weight_Clinical_Outcome %>% 
    filter(record_id != 35 & record_id != 213 & record_id != 41) %>%
    select(-record_id)
clean_weight_Clinical_Outcome <- na.omit(weight_Clinical_Outcome)
clean_weight_Clinical_Outcome$sc <- factor(clean_weight_Clinical_Outcome$sc)
pdf("Output/univariable_sc/Weight_sc.pdf")
plot1 <- ggplot(clean_weight_Clinical_Outcome,mapping=aes(x=sc,y=weight))+geom_boxplot()+scale_y_continuous()
plot2 <- ggplot(clean_weight_Clinical_Outcome)+geom_histogram(mapping=aes(x=weight),binwidth = 10)
grid.arrange(plot1,plot2,nrow=1)
dev.off()
ggplot(clean_weight_Clinical_Outcome,aes(sample=weight))+stat_qq()+stat_qq_line()
#Statistics
shapiro.test(clean_weight_Clinical_Outcome$weight)
wilcox_weight_co <- wilcox.test(weight ~ sc,clean_weight_Clinical_Outcome)$p.value
wilcox_weight_co <- round(wilcox_weight_co,digits=4)
#Odds ratio analysis
logit_weight_clinicaloutcome <- glm(sc ~ weight, data=clean_weight_Clinical_Outcome,family="binomial")
logit_weight_co <- summ(logit_weight_clinicaloutcome,exp=TRUE,digits = 4)
sink("Output/univariable_sc/Weight_regression_summary.txt")
print(logit_weight_co)
sink(file=NULL)
##################################################
stats_table <- data.frame(
    Variable=c("Gender","Age","Ethnicity","HIV status","CD4(HIV +ve patients only)","ARVS(HIV+ve patients only)",
                                     "peak ALT(>1000)","cAb HBV","chronic HBV","Drug use","Cocaine use","Methamphetamine use","Heroin use", 
                                     "baseline viral load(>800,000 IU/ml)","peak Bilirubin(>20)","Risk factor:PWID", 
                                     "Risk factor:MSM","Genotype","Diabetes","Alcohol excess","Weight","Immunotherapy","IL28B(CC vs CT/TT)"),
                          p_value=c(fisher_go_co,wilcox_age_co,fisher_eth_group_co,fisher_hiv_co,wilcox_cd4_co,fisher_arvs_co,
                                    fisher_peakALT_binary_co,fisher_cAb_HBV_co,"NA",fisher_all_drug_use_co,fisher_cocaine_use_co,fisher_meth_co,
                                    fisher_heroin_co,fisher_b_viral_load_binary_co,fisher_peakBil_binary_co,fisher_IDU_co,fisher_MSM_co,
                                    fisher_genotype_co,fisher_diabetic_co,fisher_alcohol_excess_co,wilcox_weight_co,fisher_immuno_co, 
                                    fisher_combined_il28b_co))
pdf(file="Output/univariable_sc/Assoications of variables with spontaneous clearance.pdf")
grid.table(stats_table)
dev.off()


