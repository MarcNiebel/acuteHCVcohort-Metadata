#Created by Marc Niebel Jan 2021

#Purpose of the script is to carry out exploratory data analysis
#and clean up of variables prior to modelling

#Libraries needed####
library(dplyr)
library(tidyr)
library(summarytools)


#Data for analysis####(WILL NEED TO BE CHANGED)
source("./data/AcuteHCV_R_2020-04-22_1120.r", chdir = TRUE)

#Choosing variables of interest
variables_of_interest <- data %>% 
    select(record_id,acute.factor,sc.factor,gender.factor,age,hiv.factor,bil_peak_1,
           alt_peak,risk___1.factor,hbvsag_pcr.factor,cab.factor,hiv_tx_at_hcv_diagnosis.factor,
           cd4_at_hcv_diagnosis,ifnl4_860.factor,risk___11.factor,risk___10.factor,risk___2.factor,
           drugs___10.factor:drugs___13.factor,drugs___1.factor:drugs___3.factor,drugs___20.factor,
           clinical_genotype___1.factor:clinical_genotype___4.factor)

#Remove chronic patients####
removed_chronic_patients <- variables_of_interest %>% 
    filter(acute.factor == "Yes")

#Change names of columns####
colnames(removed_chronic_patients) <- c("record_id","acute_infection","spont_clearance",
                                        "Gender","Age","hiv_pos","peak_bilirubin",
                                        "peak_ALT","PWID","chronic_HBV","core_hbv_Ab",
                                        "hiv_treatment_at_hcv_diag.","CD4_count_at_hcv_diag.",
                                        "IL28B_genotype","MSM_unknown","MSM_insertive",
                                        "MSM_receptive","Cocaine_smoked","Cocaine_injected",
                                        "Cocaine_intranasal","Cocaine_rectal","Methamphetamine_smoked",
                                        "Methamphetamine_iv","Methamphetamine_other","Heroin_iv",
                                        "clinical_gt1","clinical_gt2","clinical_gt3","clinical_gt4")

#Removing the labels on integers
removed_chronic_patients[,c(5,7:8,13)] <- sapply(removed_chronic_patients[,c(5,7:8,13)],as.numeric)


#Step 1: Visualisations of categorical and numerical variables####
#Pulling out categorical variables only
categorical_variables <- select_if(removed_chronic_patients,is.factor)

bar_plots <- function(data_fr){
    var_names <- names(data_fr)[-c(1:3)]
    for(i in seq_along(var_names)){
        plots <- ggplot(data_fr,aes_string(x="spont_clearance",fill=var_names[i]))+
            xlab("Spontaneous Clearance")+
            ylim(0,205)+
            geom_bar(position = "dodge")+
            geom_text(aes(label = ..count..), 
                      stat = "count", vjust = -0.5, 
                      colour = "black",
                      position = position_dodge(width = 1)) +
            theme_classic()
        ggsave(plots,filename=paste("Barplots_",var_names[i],".pdf"),path = "Output/univariable_sc/")
    }
}

#Generate barplots for all categorical variables seperated by clinical outcome
bar_plots(categorical_variables) 

numerical_variables <- select_if(removed_chronic_patients,is.numeric)

histogram_plots <- function(num_df,na.rm=TRUE){
    num_var <- names(num_df)
    for(i in seq_along(num_var)){
        histo_plots <- ggplot(num_df,aes_string(x=num_var[i]))+
            geom_histogram(color="white")+
            theme_classic()
        ggsave(histo_plots,filename=paste("Histograms_",num_var[i],".pdf"),path = "Output/univariable_sc/")
    }
}

#Age characteristics by clinical outcome
age_by_clinical_outcome <- removed_chronic_patients %>%
    select(Age,spont_clearance)
median_age_by_clinical_outcome <- age_by_clinical_outcome %>%
    group_by(spont_clearance) %>%
    summarise(median=median(Age))
#95% CI
progressors_age <- age_by_clinical_outcome %>%
    filter(spont_clearance =="No")
confidence_intervals_progressor <- wilcox.test(progressors_age$Age,
                                               alternative = "two.sided",
                                               conf.int = TRUE) 
spontaneous_clearers_age <- age_by_clinical_outcome %>%
    filter(spont_clearance =="Yes")
confidence_intervals_spont_clearers <- wilcox.test(spontaneous_clearers_age$Age,
                                                   alternative = "two.sided",
                                                   conf.int = TRUE)


#Transformations of data####
#Peak_bilirubin
#Using clinical cut-off of 20 mg/dL
assign_peak_bil_group <-removed_chronic_patients %>% 
    mutate(peak_bilirubin_binary=case_when(peak_bilirubin >=20 ~">= 20",
                                           peak_bilirubin < 20 ~"< 20",
                                           TRUE ~ NA_character_),
                                           peak_bilirubin_binary=as.factor(peak_bilirubin_binary)) %>%
    select(-peak_bilirubin)

#Peak ALT
#Using clinical cut-off of 1000 IU/ml
assign_peak_alt_group <- assign_peak_bil_group %>%
    mutate(Binary_peakALT=case_when(peak_ALT >=1000 ~ ">=1000",
                                    peak_ALT < 1000 ~"<1000",
                                    TRUE ~ NA_character_),
                                    Binary_peakALT =as.factor(Binary_peakALT)) %>%
    select(-peak_ALT)

#CD4 count
#Using clinical cut-off of 200 cells/mm^3(severally immunocompromised)
assign_cd4_group <- assign_peak_alt_group %>%
    mutate(Group_CD4_count=case_when(CD4_count_at_hcv_diag. < 200 ~"<200",
                           CD4_count_at_hcv_diag. >=200 ~">=200",
                           TRUE ~ NA_character_),
                            Group_CD4_count=as.factor(Group_CD4_count)) %>%
    select(-CD4_count_at_hcv_diag.)

#MSM
#Due to low numbers of the different transmission routes MSM will be amalgamated
MSM_risk <- assign_cd4_group %>% 
    mutate(MSM=case_when(MSM_unknown =="Checked"|
                         MSM_insertive =="Checked"|
                         MSM_receptive=="Checked" ~ "Yes",
                         TRUE ~ "No"),
                         MSM=as.factor(MSM)) %>%
    select(-MSM_unknown,-MSM_receptive,-MSM_insertive)

#Cocaine Usage
#Low numbers in each therefore amalgamate
Cocaine_usage <- MSM_risk %>% 
    mutate(Cocaine_use=case_when(Cocaine_smoked=="Checked"|
                                 Cocaine_injected=="Checked"|
                                 Cocaine_intranasal =="Checked"|
                                 Cocaine_smoked == "Checked" ~ "Yes",
                                 TRUE ~"No"),
                                 Cocaine_use=as.factor(Cocaine_use))%>%
    select(-Cocaine_injected,-Cocaine_smoked,-Cocaine_rectal,-Cocaine_intranasal)

#Heroin Usage
#Changing designation from heroin iv to heroin usage
Heroin_usage <- Cocaine_usage %>% 
    mutate(Heroin_use=case_when(Heroin_iv=="Checked" ~ "Yes",
                                TRUE ~"No"),
                                Heroin_use=as.factor(Heroin_use)) %>%
    select(-Heroin_iv)

#Methamphetamine Usage
#Low numbers lead to amalgamation of data
Meth_usage <- Heroin_usage %>% 
    mutate(Methamphetamine_use=case_when(Methamphetamine_smoked =="Checked"|
                             Methamphetamine_iv =="Checked"|
                             Methamphetamine_other =="Checked" ~"Yes",
                             TRUE ~ "No"),
           Methamphetamine_use=as.factor(Methamphetamine_use)) %>%
    select(-Methamphetamine_iv,-Methamphetamine_smoked,-Methamphetamine_other)

#Chronic HBV
#Chronic HBV positive has no spontaneous clearers so it will be removed
removal_chronic_hbv_variable <- Meth_usage %>% select(-chronic_HBV)

#PWID
#Changing checked/unchecked to Yes/No
removal_chronic_hbv_variable <-removal_chronic_hbv_variable %>% 
    mutate(PWID=case_when(PWID == "Checked" ~ "Yes",
                          TRUE ~ "No"))

#Infected genotype(WILL NEED TO BE REVISED DEPENDING ON NGS and CLINICAL GENOTYPE)
#Need to remove genotye 2 result as only 1 patient genotyped and also for 
#ID 66 there are two results i.e. mixed but no evidence from sequencing data
infected_genotype <- removal_chronic_hbv_variable %>% 
    select(record_id,clinical_gt1:clinical_gt4)

#Change into long format
infected_genotype_long_format <- infected_genotype %>% 
    pivot_longer(cols = clinical_gt1:clinical_gt4,names_to = "Infected_genotype")
#Only one gt 2 so remove but need for summary table
infected_genotype_long_format %>% 
    filter(value =="Checked") %>% 
    group_by(Infected_genotype) %>% 
    summarise(n=n())
#Only the infected genotype
infected_genotype_only <- infected_genotype_long_format %>% 
    filter(value =="Checked") %>%
    select(-value)
#Checking for duplicate IDs indicating mixed infection
infected_genotype_only$record_id[duplicated(infected_genotype_only$record_id)] 
#Shows ID 66 has evidence of dual infection. Not shown in sequencing data. Only 1a
removal_id66_genotype_4 <- infected_genotype_only %>%
    filter(!record_id == 66 | !Infected_genotype =="clinical_gt4")
#Join the infected data to main dataframe
transformed_data <- left_join(removal_chronic_hbv_variable,
                              removal_id66_genotype_4,by="record_id") %>%
    select(-c(clinical_gt1:clinical_gt4))

#Output of barplots for transformed data####
manipulated_data <- transformed_data %>%
    select(record_id,acute_infection,spont_clearance,peak_bilirubin_binary,Binary_peakALT,
           Group_CD4_count,MSM,Cocaine_use,Heroin_use,Methamphetamine_use,PWID,
           Infected_genotype)
bar_plots(manipulated_data)

#Summary of data separated by clinical outcome
#Removed both record_id and whether acute infection
transformed_data_removed_variables <- transformed_data[-c(1:2)]

#Removing graph column from output
st_options(dfSummary.graph.col=FALSE)

cross_tabulation <- stby(transformed_data_removed_variables,
                         transformed_data_removed_variables$spont_clearance,dfSummary)
print(cross_tabulation,file = "Output/univariable_sc/cross_tabulation.txt")

#Writing the data for downstream analysis####
#Replace genotype 2 result with NA
transformed_data_removed_gt2 <- transformed_data %>%
    mutate(Infected_genotype=case_when(Infected_genotype =="clinical_gt2" ~ NA_character_,
                                       TRUE ~ Infected_genotype))
write.csv(transformed_data_removed_gt2,file="Output/univariable_sc/transformed_acute_data.csv",row.names = FALSE)


