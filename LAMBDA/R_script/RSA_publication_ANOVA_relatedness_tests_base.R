#RSA: ANOVA for relatedness test results [for base models]===========================================
library(dplyr)
library(tidyr)
library(ez)
library(R.matlab)
library(stringr)
library(gtools)
library(WRS2)
##################################################################################
#Z-transform the tau:
#so that it becomes unbounded and tackle the unequal amount of noise between age group x brain region.
#(2018/6/10)
#Constants-------------------------------------------------------------------------
#Parameter
NUM_TRIAL=16
NUN_LAYERS=3 #Base model type, age group, brain region
LIST_NAMES_BRAIN_REGION=c("R_IPS","L_IPS","R_V1","L_V1")

#Path
PATH_RSA_ROOT="D:\\Yun-Shiuan_LAMBDA\\RSA"
PATH_RSA_CURRENT_TRIAL=file.path(PATH_RSA_ROOT,"trial_22")
PATH_RSA_ANOVA_RESULT_PARAM_BASE=file.path(PATH_RSA_CURRENT_TRIAL,"Tables",
                                     "relatedness_test_Anova","base_models")
PATH_RSA_ANOVA_RESULT_PARAM_RAW_TAU=file.path(PATH_RSA_ANOVA_RESULT_PARAM_BASE,"parametric","raw_tau")
PATH_RSA_ANOVA_RESULT_ROBUST_RAW_TAU=file.path(PATH_RSA_ANOVA_RESULT_PARAM_BASE,"robust","raw_tau")

PATH_RSA_ANOVA_RESULT_PARAM_Z_TAU=file.path(PATH_RSA_ANOVA_RESULT_PARAM_BASE,"parametric","z_tau")
PATH_RSA_ANOVA_RESULT_ROBUST_Z_TAU=file.path(PATH_RSA_ANOVA_RESULT_PARAM_BASE,"robust","z_tau")

#Could be adjusted manually
tau_type="z_tau" # or "cand_r" (the raw tau)

if(tau_type=="z_tau"){
  path_RSA_ANOVA_result_param=PATH_RSA_ANOVA_RESULT_PARAM_Z_TAU
  path_RSA_ANOVA_result_robust=PATH_RSA_ANOVA_RESULT_ROBUST_Z_TAU
}else if(tau_type=="cand_r"){
  path_RSA_ANOVA_result_param=PATH_RSA_ANOVA_RESULT_PARAM_RAW_TAU
  path_RSA_ANOVA_result_robust=PATH_RSA_ANOVA_RESULT_ROBUST_RAW_TAU
}
dir.create(path_RSA_ANOVA_result_param,recursive = T)
dir.create(path_RSA_ANOVA_result_robust,recursive = T)

#File
FILE_RSA_RESULT_BASE_MODEL_MAT=file.path(PATH_RSA_ROOT,
                                         "trial_22","Part3_2nd_order_analysis","base_models",
                                         "without_ACC","relatedness_test",
                                         "Relatedness_Test_Results_22_10-Jun-2018.mat")
FILE_RELATEDNESS_TEST_HELPER_FUNCTIONS="D:\\GoogleDrive\\Lambda_code\\R_script\\R_functions\\RSA_relatedness_test_visualization_functions.R"
#Helper functions---------------------------------------------------------------
source(FILE_RELATEDNESS_TEST_HELPER_FUNCTIONS)

#Read in relatedness results====================================================
#Read in the mat file and flatten it into a list-------------------------------
list_rsa_result=RSA_relatedness_nested_mat_to_single_mats(FILE_RSA_RESULT_BASE_MODEL_MAT,NUN_LAYERS)

#Convert the flatten list of list to a list of data frame-----------------------
list_df_rsa_result=lapply(list_rsa_result,FUN = RSA_relatedness_single_mat_to_df)

#Collect the interested subject-wise data across age group and brain region-----
df_names=names(list_df_rsa_result)
df_collect_subject_cand_r=lapply(1:length(list_df_rsa_result),
                                 FUN = function(df_index){
                                   #data to collect
                                   df_data=list_df_rsa_result[[df_index]]$df_relatedness_subjects_r
                                   df_name=df_names[df_index]
                                   df_long=
                                     df_data%>%
                                     mutate(model_type=str_extract(string=df_name,
                                                                   pattern=".*_base"),
                                            age_group=str_extract(string=df_name,
                                                                  pattern="(?<=base_).*(?=_R_|_L_)"),
                                            brain_region=str_extract(string=df_name,
                                                                     pattern="(R|L)_.*$"))%>%
                                     gather(key="subject_num",value="cand_r",matches("cand_r_\\d+"))%>%
                                     mutate(subject_num=paste0(age_group,"_",
                                                               str_extract(string=subject_num,
                                                                           pattern="(?<=_)\\d+$")))
                                   return(df_long)
                                 })
df_collect_subject_cand_r=do.call("rbind.data.frame",df_collect_subject_cand_r)
#Remove the uninterested R_DLPFC
df_collect_subject_cand_r=
  df_collect_subject_cand_r%>%
  filter(brain_region!="R_DLPFC",
         age_group!="all")
#Z-transform per age group x brain region--------------------------------------------
df_collect_subject_cand_r=
  df_collect_subject_cand_r%>%
  group_by(age_group,brain_region)%>%
  mutate(z_tau=scale(cand_r))%>%
  as.data.frame()
#Split data frame into format and magnitude models-----------------------------------------------------------------------
df_format_models=
  df_collect_subject_cand_r%>%
  filter(model_type=="format_base",
         cand_RDM_names %in% c("g","n","null"))
df_mag_abs_models=
  df_collect_subject_cand_r%>%
  filter(model_type=="mag_base",
         cand_RDM_names %in% c("agv2","a"))
df_mag_sign_models=
  df_collect_subject_cand_r%>%
  filter(model_type=="mag_base",
         cand_RDM_names %in% c("sgv2","s"))
#Parametric ANOVA=======================================================================================)============================
#Format models------------------------------------------)-------------------------------------------------
#Format models (3 way: Age group x Brain Region x Model)-----------------------------------------------------------------------
#Anova
result_anova_format_3_way=
  ezANOVA(
    data = df_format_models,
    dv = .(z_tau),
    wid = .(subject_num),
    type = 3,
    within = .(brain_region,cand_RDM_names),
    between = .(age_group),
    detailed=T
  )
# #Post-hoc t test
# result_pairwise_format=
#   pairwise.t.test(df_format_models$cand_r,df_format_models$cand_RDM_names,
#                   paired=T, p.adjust.method="bonferroni")

#Format: Collapse age group (2 way: Brain Region x Model)----------------------------------------
#Anova
result_anova_format_2_way_collapse_age_group=
  ezANOVA(
    data = df_format_models,
    dv = .(z_tau),
    wid = .(subject_num),
    type = 3,
    within = .(brain_region,cand_RDM_names),
    detailed=T
  )
#Format: Collapse age split brain region (1 way: Brain Region x Model)----------------------------------------
list_format_anova_result_collapse_age_group_split_brain_region=c()
for (brain_region_index in 1:length(LIST_NAMES_BRAIN_REGION)){
  name_brain_region=LIST_NAMES_BRAIN_REGION[brain_region_index]
  result_anova=
    ezANOVA(
      data = df_format_models%>%
        filter(brain_region==name_brain_region),
      dv = .(z_tau),
      wid = .(subject_num),
      type = 3,
      within = .(cand_RDM_names),
      detailed=T)
  list_format_anova_result_collapse_age_group_split_brain_region[[brain_region_index]]=
    result_anova
}
names(list_format_anova_result_collapse_age_group_split_brain_region)=LIST_NAMES_BRAIN_REGION
#Format: Collapse age: Pairwise comparison for model effect------------------------------------------------
#All pairs-----------------------------------------------
list_format_pariwise_comparison_collapse_age_group_split_brain_region=c()

for (brain_region_index in 1:length(LIST_NAMES_BRAIN_REGION)){
  name_brain_region=LIST_NAMES_BRAIN_REGION[brain_region_index]
  list_pariwise_comparison_split_brain_region=c()#Collect the age group results per brain region
  
  #Order the model by the cand_r value
  df_model_ordered_level=
    df_format_models%>%
    filter(brain_region==name_brain_region)%>%
    group_by(cand_RDM_names)%>%
    rename(dv=tau_type)%>%
    summarise(mean_dv=mean(dv))%>%
    arrange(desc(mean_dv))%>%
    mutate(cand_RDM_names=as.character(cand_RDM_names),
           cand_full_RDM_names=
             case_when(cand_RDM_names=="g" ~ "G-Format",
                       cand_RDM_names=="n" ~ "N-format",
                       cand_RDM_names=="null" ~ "Null",
                       TRUE~cand_RDM_names))
  result_pairwise_format=
    pairwise.t.test(df_format_models%>%
                      rename(dv=tau_type)%>%
                      filter(brain_region==name_brain_region)%>%
                      pull(dv),
                    df_format_models%>%
                      filter(brain_region==name_brain_region)%>%
                      pull(cand_RDM_names)%>%
                      factor(levels=df_model_ordered_level$cand_RDM_names,
                             labels = df_model_ordered_level$cand_full_RDM_names),
                    paired=T, p.adjust.method="bonferroni")
  
  list_format_pariwise_comparison_collapse_age_group_split_brain_region[[brain_region_index]]=
    result_pairwise_format
}
names(list_format_pariwise_comparison_collapse_age_group_split_brain_region)=LIST_NAMES_BRAIN_REGION

#Only focus on the best model (compare it to the others)--------------------------------------
list_format_best_model_comparison_collapse_age_group_split_brain_region=c()

for (brain_region_index in 1:length(LIST_NAMES_BRAIN_REGION)){
  name_brain_region=LIST_NAMES_BRAIN_REGION[brain_region_index]
  #Order the model by the cand_r value
  df_model_ordered_level=
    df_format_models%>%
    filter(brain_region==name_brain_region)%>%
    group_by(cand_RDM_names)%>%
    rename(dv=tau_type)%>%
    summarise(mean_dv=mean(dv))%>%
    arrange(desc(mean_dv))%>%
    mutate(cand_RDM_names=as.character(cand_RDM_names),
           cand_full_RDM_names=
             case_when(cand_RDM_names=="g" ~ "G-Format",
                       cand_RDM_names=="n" ~ "N-format",
                       cand_RDM_names=="null" ~ "Null",
                       TRUE~cand_RDM_names))
  result_pairwise_format=
    pairwise.t.test(df_format_models%>%
                      filter(brain_region==name_brain_region)%>%
                      rename(dv=tau_type)%>%
                      pull(dv),
                    df_format_models%>%
                      filter(brain_region==name_brain_region)%>%
                      pull(cand_RDM_names)%>%
                      factor(levels=df_model_ordered_level$cand_RDM_names,
                             labels = df_model_ordered_level$cand_full_RDM_names),
                    paired=T, p.adjust.method="none",alternative = "less")#Remember to correct the p value afterwards
  list_format_best_model_comparison_collapse_age_group_split_brain_region[[brain_region_index]]=
    result_pairwise_format
}
names(list_format_best_model_comparison_collapse_age_group_split_brain_region)=LIST_NAMES_BRAIN_REGION
#Absolute magnitude models------------------------------------------)-------------------------------------------------
#Magnitude models - Absolute magnitude-----------------------------------------------------------------------
#Abs models (3 way: Age group x Brain Region x Model)-----------------------------------------------------------------------
#Anova
result_anova_mag_abs_3_way=
  ezANOVA(
    data = df_mag_abs_models,
    dv = .(z_tau),
    wid = .(subject_num),
    type = 3,
    within = .(brain_region,cand_RDM_names),
    between = .(age_group),
    detailed=T
  )

#Abs: Collapse age group (2 way: Brain Region x Model)----------------------------------------
#Anova
result_anova_mag_abs_2_way_collapse_age_group=
  ezANOVA(
    data = df_mag_abs_models,
    dv = .(z_tau),
    wid = .(subject_num),
    type = 3,
    within = .(brain_region,cand_RDM_names),
    detailed=T
  )
#Abs: Collapse age split brain region (1 way: Brain Region x Model)----------------------------------------
list_mag_abs_anova_result_collapse_age_group_split_brain_region=c()
for (brain_region_index in 1:length(LIST_NAMES_BRAIN_REGION)){
  name_brain_region=LIST_NAMES_BRAIN_REGION[brain_region_index]
  result_anova=
    ezANOVA(
      data = df_mag_abs_models%>%
        filter(brain_region==name_brain_region),
      dv = .(z_tau),
      wid = .(subject_num),
      type = 3,
      within = .(cand_RDM_names),
      detailed=T)
  list_mag_abs_anova_result_collapse_age_group_split_brain_region[[brain_region_index]]=
    result_anova
}
names(list_mag_abs_anova_result_collapse_age_group_split_brain_region)=LIST_NAMES_BRAIN_REGION
#Abs: Collapse age: Pairwise comparison for model effect------------------------------------------------
#All pairs-----------------------------------------------
list_mag_abs_pariwise_comparison_collapse_age_group_split_brain_region=c()

for (brain_region_index in 1:length(LIST_NAMES_BRAIN_REGION)){
  name_brain_region=LIST_NAMES_BRAIN_REGION[brain_region_index]
  list_pariwise_comparison_split_brain_region=c()#Collect the age group results per brain region
  
  #Order the model by the cand_r value
  df_model_ordered_level=
    df_mag_abs_models%>%
    filter(brain_region==name_brain_region)%>%
    group_by(cand_RDM_names)%>%
    rename(dv=tau_type)%>%
    summarise(mean_dv=mean(dv))%>%
    arrange(desc(mean_dv))%>%
    mutate(cand_RDM_names=as.character(cand_RDM_names),
           cand_full_RDM_names=
             case_when(cand_RDM_names=="g" ~ "G-Format",
                       cand_RDM_names=="n" ~ "N-format",
                       cand_RDM_names=="null" ~ "Null",
                       TRUE~cand_RDM_names))
  result_pairwise_format=
    pairwise.t.test(df_mag_abs_models%>%
                      rename(dv=tau_type)%>%
                      filter(brain_region==name_brain_region)%>%
                      pull(dv),
                    df_mag_abs_models%>%
                      filter(brain_region==name_brain_region)%>%
                      pull(cand_RDM_names)%>%
                      factor(levels=df_model_ordered_level$cand_RDM_names,
                             labels = df_model_ordered_level$cand_full_RDM_names),
                    paired=T, p.adjust.method="bonferroni")
  
  list_mag_abs_pariwise_comparison_collapse_age_group_split_brain_region[[brain_region_index]]=
    result_pairwise_format
}
names(list_mag_abs_pariwise_comparison_collapse_age_group_split_brain_region)=LIST_NAMES_BRAIN_REGION

#Only focus on the best model (compare it to the others)--------------------------------------
list_mag_abs_best_model_comparison_collapse_age_group_split_brain_region=c()

for (brain_region_index in 1:length(LIST_NAMES_BRAIN_REGION)){
  name_brain_region=LIST_NAMES_BRAIN_REGION[brain_region_index]
  #Order the model by the cand_r value
  df_model_ordered_level=
    df_mag_abs_models%>%
    filter(brain_region==name_brain_region)%>%
    group_by(cand_RDM_names)%>%
    rename(dv=tau_type)%>%
    summarise(mean_dv=mean(dv))%>%
    arrange(desc(mean_dv))%>%
    mutate(cand_RDM_names=as.character(cand_RDM_names),
           cand_full_RDM_names=
             case_when(cand_RDM_names=="g" ~ "G-Format",
                       cand_RDM_names=="n" ~ "N-format",
                       cand_RDM_names=="null" ~ "Null",
                       TRUE~cand_RDM_names))
  result_pairwise_format=
    pairwise.t.test(df_mag_abs_models%>%
                      filter(brain_region==name_brain_region)%>%
                      rename(dv=tau_type)%>%
                      pull(dv),
                    df_mag_abs_models%>%
                      filter(brain_region==name_brain_region)%>%
                      pull(cand_RDM_names)%>%
                      factor(levels=df_model_ordered_level$cand_RDM_names,
                             labels = df_model_ordered_level$cand_full_RDM_names),
                    paired=T, p.adjust.method="none",alternative = "less")#Remember to correct the p value afterwards
  list_mag_abs_best_model_comparison_collapse_age_group_split_brain_region[[brain_region_index]]=
    result_pairwise_format
}
names(list_mag_abs_best_model_comparison_collapse_age_group_split_brain_region)=LIST_NAMES_BRAIN_REGION
#Signed magnitude models---------------------------------------------)-------------------------------------------------
#Magnitude models - Signed magnitude-----------------------------------------------------------------------
#Anova
result_anova_mag_sign=
  ezANOVA(
    data = df_mag_sign_models,
    dv = .(z_tau),
    wid = .(subject_num),
    type = 3,
    within = .(brain_region,cand_RDM_names),
    between = .(age_group),
    detailed=T
  )
#Print out results as csv==================================================================================)----------------------------------------------------------------------------------
#Format models------------------------------------------)-------------------------------------------------
#Format models(3 way: Age group x Brain Region x Model)------------------------------------------------------------------------------
do.call("smartbind",c(result_anova_format_3_way,fill=""))%>%
  tibble::rownames_to_column()%>%
  mutate(Effect=gsub(x=Effect,pattern="_",replacement=" "),
         Effect=gsub(x=Effect,pattern=":",replacement=" x "),
         Effect=gsub(x=Effect,pattern="age group",replacement="Age Group"),
         Effect=gsub(x=Effect,pattern="brain region",replacement="Brain Region"),
         Effect=gsub(x=Effect,pattern="cand RDM names",replacement="Model"))%>%
  mutate_if(is.numeric,round,digits=3)%>%
  mutate_at(vars(matches("p|p[GG]|p[HF]")),
            funs(gsub(x=.,pattern="^0$",replacement="<.001")))%>%
  # mutate_at(vars(matches("p|p[GG]|p[HF]")),
  #           funs(gsub(x=.,pattern="^0(?=\\.)",replacement="",perl = T)))%>%
  write.csv(file = file.path(path_RSA_ANOVA_result_param,"Format_Way_3_R~agexbrainxmodel_anova_table.csv"))
#Format: Collapse age group (2 way: Brain Region x Model)-----------------------------
do.call("smartbind",c(result_anova_format_2_way_collapse_age_group,fill=""))%>%
  tibble::rownames_to_column()%>%
  mutate(Effect=gsub(x=Effect,pattern="_",replacement=" "),
         Effect=gsub(x=Effect,pattern=":",replacement=" x "),
         Effect=gsub(x=Effect,pattern="age group",replacement="Age Group"),
         Effect=gsub(x=Effect,pattern="brain region",replacement="Brain Region"),
         Effect=gsub(x=Effect,pattern="cand RDM names",replacement="Model"))%>%
  mutate_if(is.numeric,round,digits=3)%>%
  mutate_at(vars(matches("p|p[GG]|p[HF]")),
            funs(gsub(x=.,pattern="^0$",replacement="<.001")))%>%
  # mutate_at(vars(matches("p|p[GG]|p[HF]")),
  #           funs(gsub(x=.,pattern="^0(?=\\.)",replacement="",perl = T)))%>%
  write.csv(file = file.path(path_RSA_ANOVA_result_param,"Format_Collapse_age_group_Way_2_R~formatxbrainxmodel_anova_table.csv"))
#Format: Collapse age split brain region (1 way: Brain Region x Model)----------------
for (brain_region_index in 1:length(LIST_NAMES_BRAIN_REGION)){
  name_brain_region=LIST_NAMES_BRAIN_REGION[brain_region_index]
  anova_result=list_format_anova_result_collapse_age_group_split_brain_region[[brain_region_index]]
  
  do.call("smartbind",c(anova_result,fill=""))%>%
    tibble::rownames_to_column()%>%
    mutate(Effect=gsub(x=Effect,pattern="_",replacement=" "),
           Effect=gsub(x=Effect,pattern=":",replacement=" x "),
           Effect=gsub(x=Effect,pattern="age group",replacement="Age Group"),
           Effect=gsub(x=Effect,pattern="brain region",replacement="Brain Region"),
           Effect=gsub(x=Effect,pattern="cand RDM names",replacement="Model"))%>%
    mutate_if(is.numeric,round,digits=3)%>%
    mutate_at(vars(matches("p|p[GG]|p[HF]")),
              funs(gsub(x=.,pattern="^0$",replacement="<.001")))%>%
    # mutate_at(vars(matches("p|p[GG]|p[HF]")),
    #           funs(gsub(x=.,pattern="^0(?=\\.)",replacement="",perl = T)))%>%
    write.csv(file = file.path(path_RSA_ANOVA_result_param,
                               paste0("Format_Collapse_age_group_Way_1_R~model_anova_table_",
                                      name_brain_region,".csv")))
  
}
#Format: Collapse age: Pairwise comparison for model effect)--------------------------
#All pairs----------------------------------------------------------------------------
for (brain_region_index in 1:length(LIST_NAMES_BRAIN_REGION)){
  name_brain_region=LIST_NAMES_BRAIN_REGION[brain_region_index]
  list_format_pariwise_comparison_collapse_age_group_split_brain_region[[name_brain_region]]$p.value%>%
    as.data.frame()%>%
    tibble::rownames_to_column(var="Model")%>%
    mutate_if(is.numeric,round,digits=3)%>%
    mutate_if(is.numeric,gsub,pattern="^0$",replacement="<.001")%>%
    write.csv(file = file.path(path_RSA_ANOVA_result_param,
                               paste0("Format_Collapse_age_group_Pairwise_comparison_",
                                      name_brain_region,".csv")))
  
}
#Only focus on the best model (compare it to the other)-------------------------------
for (brain_region_index in 1:length(LIST_NAMES_BRAIN_REGION)){
  name_brain_region=LIST_NAMES_BRAIN_REGION[brain_region_index]
  list_format_best_model_comparison_collapse_age_group_split_brain_region[[name_brain_region]]$p.value%>%
    as.data.frame()%>%
    tibble::rownames_to_column(var="Model")%>%
    select(1:2)%>%#The first two columns (model name and the best model)
    mutate_if(is.numeric,p.adjust,"bonferroni")%>%#Correct the raw p-value
    mutate_if(is.numeric,round,digits=3)%>%
    mutate_if(is.numeric,gsub,pattern="^0$",replacement="<.001")%>%
    write.csv(file = file.path(path_RSA_ANOVA_result_param,
                               paste0("Format_Collapse_age_group_Pair_to_best_model_",
                                      name_brain_region,".csv")))
}     
#Absolute magnitude models------------------------------------------)-------------------------------------------------
#Abs mag models(3 way: Age group x Brain Region x Model)------------------------------------------------------------------------------
do.call("smartbind",c(result_anova_mag_abs_3_way,fill=""))%>%
  tibble::rownames_to_column()%>%
  mutate(Effect=gsub(x=Effect,pattern="_",replacement=" "),
         Effect=gsub(x=Effect,pattern=":",replacement=" x "),
         Effect=gsub(x=Effect,pattern="age group",replacement="Age Group"),
         Effect=gsub(x=Effect,pattern="brain region",replacement="Brain Region"),
         Effect=gsub(x=Effect,pattern="cand RDM names",replacement="Model"))%>%
  mutate_if(is.numeric,round,digits=3)%>%
  mutate_at(vars(matches("p|p[GG]|p[HF]")),
            funs(gsub(x=.,pattern="^0$",replacement="<.001")))%>%
  # mutate_at(vars(matches("p|p[GG]|p[HF]")),
  #           funs(gsub(x=.,pattern="^0(?=\\.)",replacement="",perl = T)))%>%
  write.csv(file = file.path(path_RSA_ANOVA_result_param,"Abs_mag_Way_3_R~agexbrainxmodel_anova_table.csv"))
#Abs mag: Collapse age group (2 way: Brain Region x Model)-----------------------------
do.call("smartbind",c(result_anova_mag_abs_2_way_collapse_age_group,fill=""))%>%
  tibble::rownames_to_column()%>%
  mutate(Effect=gsub(x=Effect,pattern="_",replacement=" "),
         Effect=gsub(x=Effect,pattern=":",replacement=" x "),
         Effect=gsub(x=Effect,pattern="age group",replacement="Age Group"),
         Effect=gsub(x=Effect,pattern="brain region",replacement="Brain Region"),
         Effect=gsub(x=Effect,pattern="cand RDM names",replacement="Model"))%>%
  mutate_if(is.numeric,round,digits=3)%>%
  mutate_at(vars(matches("p|p[GG]|p[HF]")),
            funs(gsub(x=.,pattern="^0$",replacement="<.001")))%>%
  # mutate_at(vars(matches("p|p[GG]|p[HF]")),
  #           funs(gsub(x=.,pattern="^0(?=\\.)",replacement="",perl = T)))%>%
  write.csv(file = file.path(path_RSA_ANOVA_result_param,"Abs_mag_Collapse_age_group_Way_2_R~formatxbrainxmodel_anova_table.csv"))
#Abs mag: Collapse age split brain region (1 way: Brain Region x Model)----------------
for (brain_region_index in 1:length(LIST_NAMES_BRAIN_REGION)){
  name_brain_region=LIST_NAMES_BRAIN_REGION[brain_region_index]
  anova_result=list_mag_abs_anova_result_collapse_age_group_split_brain_region[[brain_region_index]]
  
  do.call("smartbind",c(anova_result,fill=""))%>%
    tibble::rownames_to_column()%>%
    mutate(Effect=gsub(x=Effect,pattern="_",replacement=" "),
           Effect=gsub(x=Effect,pattern=":",replacement=" x "),
           Effect=gsub(x=Effect,pattern="age group",replacement="Age Group"),
           Effect=gsub(x=Effect,pattern="brain region",replacement="Brain Region"),
           Effect=gsub(x=Effect,pattern="cand RDM names",replacement="Model"))%>%
    mutate_if(is.numeric,round,digits=3)%>%
    mutate_at(vars(matches("p|p[GG]|p[HF]")),
              funs(gsub(x=.,pattern="^0$",replacement="<.001")))%>%
    # mutate_at(vars(matches("p|p[GG]|p[HF]")),
    #           funs(gsub(x=.,pattern="^0(?=\\.)",replacement="",perl = T)))%>%
    write.csv(file = file.path(path_RSA_ANOVA_result_param,
                               paste0("Abs_mag_Collapse_age_group_Way_1_R~model_anova_table_",
                                      name_brain_region,".csv")))
  
}
#Abs mag: Collapse age: Pairwise comparison for model effect)--------------------------
#All pairs----------------------------------------------------------------------------
for (brain_region_index in 1:length(LIST_NAMES_BRAIN_REGION)){
  name_brain_region=LIST_NAMES_BRAIN_REGION[brain_region_index]
  list_mag_abs_pariwise_comparison_collapse_age_group_split_brain_region[[name_brain_region]]$p.value%>%
    as.data.frame()%>%
    tibble::rownames_to_column(var="Model")%>%
    mutate_if(is.numeric,round,digits=3)%>%
    mutate_if(is.numeric,gsub,pattern="^0$",replacement="<.001")%>%
    write.csv(file = file.path(path_RSA_ANOVA_result_param,
                               paste0("Abs_mag_Collapse_age_group_Pairwise_comparison_",
                                      name_brain_region,".csv")))
  
}
#Only focus on the best model (compare it to the other)-------------------------------
for (brain_region_index in 1:length(LIST_NAMES_BRAIN_REGION)){
  name_brain_region=LIST_NAMES_BRAIN_REGION[brain_region_index]
  list_mag_abs_best_model_comparison_collapse_age_group_split_brain_region[[name_brain_region]]$p.value%>%
    as.data.frame()%>%
    tibble::rownames_to_column(var="Model")%>%
    select(1:2)%>%#The first two columns (model name and the best model)
    mutate_if(is.numeric,p.adjust,"bonferroni")%>%#Correct the raw p-value
    mutate_if(is.numeric,round,digits=3)%>%
    mutate_if(is.numeric,gsub,pattern="^0$",replacement="<.001")%>%
    write.csv(file = file.path(path_RSA_ANOVA_result_param,
                               paste0("Abs_mag_Collapse_age_group_Pair_to_best_model_",
                                      name_brain_region,".csv")))
}     
#Magnitude models - abs-------------------------------------------------------------------
do.call("smartbind",c(result_anova_mag_abs,fill=""))%>%
  tibble::rownames_to_column()%>%
  mutate(Effect=gsub(x=Effect,pattern="_",replacement=" "),
         Effect=gsub(x=Effect,pattern=":",replacement=" x "),
         Effect=gsub(x=Effect,pattern="age group",replacement="Age Group"),
         Effect=gsub(x=Effect,pattern="brain region",replacement="Brain Region"),
         Effect=gsub(x=Effect,pattern="cand RDM names",replacement="Model"))%>%
  mutate_if(is.numeric,round,digits=3)%>%
  mutate_at(vars(matches("p|p[GG]|p[HF]")),
            funs(gsub(x=.,pattern="^0$",replacement="<.001")))%>%
  write.csv(file = file.path(path_RSA_ANOVA_result_param,"mag_abs_models_anova_table.csv"))
#Signed magnitude models------------------------------------------)-------------------------------------------------
#Magnitude models - sign------------------------------------------------------------------
do.call("smartbind",c(result_anova_mag_sign,fill=""))%>%
  tibble::rownames_to_column()%>%
  mutate(Effect=gsub(x=Effect,pattern="_",replacement=" "),
         Effect=gsub(x=Effect,pattern=":",replacement=" x "),
         Effect=gsub(x=Effect,pattern="age group",replacement="Age Group"),
         Effect=gsub(x=Effect,pattern="brain region",replacement="Brain Region"),
         Effect=gsub(x=Effect,pattern="cand RDM names",replacement="Model"))%>%
  mutate_if(is.numeric,round,digits=3)%>%
  mutate_at(vars(matches("p|p[GG]|p[HF]")),
            funs(gsub(x=.,pattern="^0$",replacement="<.001")))%>%
  write.csv(file = file.path(path_RSA_ANOVA_result_param,"mag_sign_models_anova_table.csv"))

#Robust version ANOVA=========================================================================================)===========================
#Format: Collapse age group split brain region (1 way: Brain Region x Model) (1w: trimmed)----------------------------------------------
list_format_anova_result_collapse_age_group_split_brain_region_trimmed=c()

for (brain_region_index in 1:length(LIST_NAMES_BRAIN_REGION)){
  name_brain_region=LIST_NAMES_BRAIN_REGION[brain_region_index]
  df=
    df_format_models%>%
    filter(brain_region==name_brain_region)%>%
    rename(dv=tau_type)
  result_anov_trimmed=
    rmanova(y = df$dv,
            groups = df$cand_RDM_names,
            blocks = df$subject_num)
  list_format_anova_result_collapse_age_group_split_brain_region_trimmed[[brain_region_index]]=result_anov_trimmed
}

names(list_format_anova_result_collapse_age_group_split_brain_region_trimmed)=LIST_NAMES_BRAIN_REGION
#Format: Collapse age group: Post-hoc: RDM correlation ~ Model (post-hoc: trimmed)--------------------------------------------
#Format Pairwise: Only focus on the best model (compare it to the others)--------------------------------------
list_format_best_model_comparison_collapse_age_group_split_brain_region_trimmed=c()
#The p criteria for Rom's method (used by the package and recommended by Wilcox, R. (2012))
list_rom_correction=c(0.05000,0.02500,0.01690,0.01270,0.01020,
                      0.00851,0.00730,0.00639,0.00568,0.00511)
list_ordered_RDM=c("G-Format","N-Format","Null","log-A","log-S","linear-A","linear-S")
for (brain_region_index in 1:length(LIST_NAMES_BRAIN_REGION)){
  name_brain_region=LIST_NAMES_BRAIN_REGION[brain_region_index]
    df=
      df_format_models%>%
      filter(brain_region==name_brain_region)%>%
      rename(dv=tau_type)
    #Only retain the comparison between the best model and the rest
    result_anova=
      rmmcp(y = df$dv,
            groups = df$cand_RDM_names,
            blocks = df$subject_num)
    
    dimnames(result_anova$comp)[[2]][1]="Group1"
    dimnames(result_anova$comp)[[2]][2]="Group2"
    
    result_anova$comp=
      result_anova$comp%>%
      as.data.frame()%>%
      filter(Group1==1)%>%
      arrange(desc(p.value))%>%
      mutate(p.crit=list_rom_correction[1:n()])%>%
      arrange(Group2)%>%
      mutate(Group1=result_anova$fnames[1],
             Group2=result_anova$fnames[2:length(result_anova$fnames)],
             Significant_Rom=ifelse(p.value<p.crit,yes = "*",no = ""),
             Group1=case_when(Group1=="g" ~ "G-Format",
                              Group1=="n" ~ "N-Format",
                              Group1=="null" ~ "Null",
                              TRUE~Group1
             ),
             Group2=case_when(Group2=="g" ~ "G-Format",
                              Group2=="n" ~ "N-Format",
                              Group2=="null" ~ "Null",
                              TRUE~Group2))%>%
      rename(p.value.raw=p.value)%>%
      mutate(p.bonferroni=p.value.raw*(length(result_anova$fnames)-1),
             p.bonferroni=ifelse(p.bonferroni>1,yes = 1,no = p.bonferroni))%>%
      mutate_if(is.numeric,round,3)%>%
      mutate_at(vars(matches("^p\\..*")),
                funs(gsub(x=.,pattern="^0$",replacement="<.001")))%>%
      mutate(Group2=factor(Group2,levels = list_ordered_RDM))%>%
      arrange(Group2)
    list_format_best_model_comparison_collapse_age_group_split_brain_region_trimmed[[brain_region_index]]=result_anova
}
names(list_format_best_model_comparison_collapse_age_group_split_brain_region_trimmed)=LIST_NAMES_BRAIN_REGION


#Abs mag: Collapse age group split brain region (1 way: Brain Region x Model) (1w: trimmed)----------------------------------------------
list_abs_mag_anova_result_collapse_age_group_split_brain_region_trimmed=c()

for (brain_region_index in 1:length(LIST_NAMES_BRAIN_REGION)){
  name_brain_region=LIST_NAMES_BRAIN_REGION[brain_region_index]
  df=
    df_mag_abs_models%>%
    filter(brain_region==name_brain_region)%>%
    rename(dv=tau_type)
  result_anov_trimmed=
    rmanova(y = df$dv,
            groups = df$cand_RDM_names,
            blocks = df$subject_num)
  list_abs_mag_anova_result_collapse_age_group_split_brain_region_trimmed[[brain_region_index]]=result_anov_trimmed
}
names(list_abs_mag_anova_result_collapse_age_group_split_brain_region_trimmed)=LIST_NAMES_BRAIN_REGION
#Abs mag: Collapse age group: Post-hoc: RDM correlation ~ Model (post-hoc: trimmed)--------------------------------------------
#Abs mag Pairwise: Only focus on the best model (compare it to the others)--------------------------------------
list_mag_abs_best_model_comparison_collapse_age_group_split_brain_region_trimmed=c()
#The p criteria for Rom's method (used by the package and recommended by Wilcox, R. (2012))
list_rom_correction=c(0.05000,0.02500,0.01690,0.01270,0.01020,
                      0.00851,0.00730,0.00639,0.00568,0.00511)
list_ordered_RDM=c("G-Format","N-Format","Null","log-A","log-S","linear-A","linear-S")
for (brain_region_index in 1:length(LIST_NAMES_BRAIN_REGION)){
  name_brain_region=LIST_NAMES_BRAIN_REGION[brain_region_index]
  df=
    df_mag_abs_models%>%
    filter(brain_region==name_brain_region)%>%
    rename(dv=tau_type)
  #Only retain the comparison between the best model and the rest
  result_anova=
    rmmcp(y = df$dv,
          groups = df$cand_RDM_names,
          blocks = df$subject_num)
  #Fix the bug (when there's only two groups being compared)
  dimnames(result_anova$comp)[[2]][6]="p.value"
  result_anova$comp=
    cbind(result_anova$comp,list_rom_correction[nrow(result_anova$comp)])
  dimnames(result_anova$comp)[[2]][7]="p.crit"
  
  
  dimnames(result_anova$comp)[[2]][1]="Group1"
  dimnames(result_anova$comp)[[2]][2]="Group2"
  
  result_anova$comp=
    result_anova$comp%>%
    as.data.frame()%>%
    filter(Group1==1)%>%
    arrange(desc(p.value))%>%
    mutate(p.crit=list_rom_correction[1:n()])%>%
    arrange(Group2)%>%
    mutate(Group1=result_anova$fnames[1],
           Group2=result_anova$fnames[2:length(result_anova$fnames)],
           Significant_Rom=ifelse(p.value<p.crit,yes = "*",no = ""),
           Group1=case_when(Group1=="agv2" ~ "log-A",
                            Group1=="a" ~ "linear-A",
                            TRUE~Group1
           ),
           Group2=case_when(Group2=="agv2" ~ "log-A",
                            Group2=="a" ~ "linear-A",
                            TRUE~Group2))%>%
    rename(p.value.raw=p.value)%>%
    mutate(p.bonferroni=p.value.raw*(length(result_anova$fnames)-1),
           p.bonferroni=ifelse(p.bonferroni>1,yes = 1,no = p.bonferroni))%>%
    mutate_if(is.numeric,round,3)%>%
    mutate_at(vars(matches("^p\\..*")),
              funs(gsub(x=.,pattern="^0$",replacement="<.001")))%>%
    mutate(Group2=factor(Group2,levels = list_ordered_RDM))%>%
    arrange(Group2)
  list_mag_abs_best_model_comparison_collapse_age_group_split_brain_region_trimmed[[brain_region_index]]=result_anova
}
names(list_mag_abs_best_model_comparison_collapse_age_group_split_brain_region_trimmed)=LIST_NAMES_BRAIN_REGION


#Print out results as csv=====================================================================================)=======================================================
if(tau_type=="cand_r"){
  PATH_RSA_ANOVA_RESULT_ROBUST=PATH_RSA_ANOVA_RESULT_ROBUST_RAW_TAU
}else if(tau_type=="z_tau"){
  PATH_RSA_ANOVA_RESULT_ROBUST=PATH_RSA_ANOVA_RESULT_ROBUST_Z_TAU
}
#Format: Collapse age group split brain region (1 way: Brain Region x Model) (1w: trimmed)-------------------------------------
for (brain_region_index in 1:length(LIST_NAMES_BRAIN_REGION)){
  name_brain_region=LIST_NAMES_BRAIN_REGION[brain_region_index]
  anova_result=list_format_anova_result_collapse_age_group_split_brain_region_trimmed[[brain_region_index]]
  
  data.frame(Effect="Model",
             DFn=anova_result$df1,
             DFd=anova_result$df2,
             `F`=anova_result$test,
             p=anova_result$p.value)%>%
    mutate_if(is.numeric,round,3)%>%
    mutate_at(vars(matches("p")),
              funs(gsub(x=.,pattern="^0$",replacement="<.001")))%>%
    write.csv(file = file.path(path_RSA_ANOVA_result_robust,
                               paste0("Format_Collapse_age_group_Way_1_R~model_anova_table_",
                                      name_brain_region,".csv")))
  
}
#Format: Collapse age group: Post-hoc: RDM correlation ~ Model (post-hoc: trimmed)--------------------------------------------------------------------------
#Format: Pairwise: Only focus on the best model (compare it to the others)-----------------------------------------------
for (brain_region_index in 1:length(LIST_NAMES_BRAIN_REGION)){
  name_brain_region=LIST_NAMES_BRAIN_REGION[brain_region_index]
  anova_result=list_format_best_model_comparison_collapse_age_group_split_brain_region_trimmed[[name_brain_region]]
    
  anova_result$comp%>%
    write.csv(file = file.path(PATH_RSA_ANOVA_RESULT_ROBUST,
                               paste0("Format_Collapse_age_group_Pair_to_best_model_",
                                      name_brain_region,".csv")))
  
}


#Abs mag: Collapse age group split brain region (1 way: Brain Region x Model) (1w: trimmed)-------------------------------------
for (brain_region_index in 1:length(LIST_NAMES_BRAIN_REGION)){
  name_brain_region=LIST_NAMES_BRAIN_REGION[brain_region_index]
  anova_result=list_abs_mag_anova_result_collapse_age_group_split_brain_region_trimmed[[brain_region_index]]
  
  data.frame(Effect="Model",
             DFn=anova_result$df1,
             DFd=anova_result$df2,
             `F`=anova_result$test,
             p=anova_result$p.value)%>%
    mutate_if(is.numeric,round,3)%>%
    mutate_at(vars(matches("p")),
              funs(gsub(x=.,pattern="^0$",replacement="<.001")))%>%
    write.csv(file = file.path(path_RSA_ANOVA_result_robust,
                               paste0("Abs_mag_Collapse_age_group_Way_1_R~model_anova_table_",
                                      name_brain_region,".csv")))
}
#Abs mag: Collapse age group: Post-hoc: RDM correlation ~ Model (post-hoc: trimmed)--------------------------------------------------------------------------
#Abs mag: Pairwise: Only focus on the best model (compare it to the others)-----------------------------------------------
for (brain_region_index in 1:length(LIST_NAMES_BRAIN_REGION)){
  name_brain_region=LIST_NAMES_BRAIN_REGION[brain_region_index]
  anova_result=list_mag_abs_best_model_comparison_collapse_age_group_split_brain_region_trimmed[[name_brain_region]]
  
  anova_result$comp%>%
    write.csv(file = file.path(PATH_RSA_ANOVA_RESULT_ROBUST,
                               paste0("Abs_mag_Collapse_age_group_Pair_to_best_model_",
                                      name_brain_region,".csv")))
  
}

