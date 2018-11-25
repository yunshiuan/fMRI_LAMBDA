#RSA: ANOVA for relatedness test results[for base models]===========================================
library(dplyr)
library(tidyr)
library(ez)
library(R.matlab)
library(stringr)
library(gtools)
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
PATH_RSA_ANOVA_RESULT_BASE=file.path(PATH_RSA_CURRENT_TRIAL,"Tables",
                                "relatedness_test_Anova","base_models")
PATH_RSA_ANOVA_RESULT_PARAM_RAW_TAU=file.path(PATH_RSA_ANOVA_RESULT_BASE,"parametric","raw_tau")
PATH_RSA_ANOVA_RESULT_ROBUST_RAW_TAU=file.path(PATH_RSA_ANOVA_RESULT_BASE,"robust","raw_tau")
PATH_RSA_ANOVA_RESULT_PARAM_Z_TAU=file.path(PATH_RSA_ANOVA_RESULT_BASE,"parametric","z_tau")
PATH_RSA_ANOVA_RESULT_ROBUST_Z_TAU=file.path(PATH_RSA_ANOVA_RESULT_BASE,"robust","z_tau")

dir.create(PATH_RSA_ANOVA_RESULT,recursive = T)
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
#Raw tau: Split into format and magnitude models========================================================)=========================
#Format models-----------------------------------------------------------------------
df_format_models=
  df_collect_subject_cand_r%>%
    filter(model_type=="format_base",
           cand_RDM_names %in% c("g","n","null"))
#Anova
result_anova_format=
  ezANOVA(
  data = df_format_models,
  dv = .(cand_r),
  wid = .(subject_num),
  type = 3,
  within = .(brain_region,cand_RDM_names),
  between = .(age_group),
  detailed=T
)
#Post-hoc t test
result_pairwise_format=
  pairwise.t.test(df_format_models$cand_r,df_format_models$cand_RDM_names,
                  paired=T, p.adjust.method="bonferroni")

#Magnitude models - Absolute magnitude-----------------------------------------------------------------------
df_mag_abs_models=
  df_collect_subject_cand_r%>%
  filter(model_type=="mag_base",
         cand_RDM_names %in% c("agv2","a"))
#Anova
result_anova_mag_abs=
  ezANOVA(
    data = df_mag_abs_models,
    dv = .(cand_r),
    wid = .(subject_num),
    type = 3,
    within = .(brain_region,cand_RDM_names),
    between = .(age_group),
    detailed=T
  )
#Magnitude models - Signed magnitude-----------------------------------------------------------------------
df_mag_sign_models=
  df_collect_subject_cand_r%>%
  filter(model_type=="mag_base",
         cand_RDM_names %in% c("sgv2","s"))
#Anova
result_anova_mag_sign=
  ezANOVA(
    data = df_mag_sign_models,
    dv = .(cand_r),
    wid = .(subject_num),
    type = 3,
    within = .(brain_region,cand_RDM_names),
    between = .(age_group),
    detailed=T
  )
#Print out results as csv---------------------------------------------------)----------------------------------------------------------------------------------
#Format models
do.call("smartbind",c(result_anova_format,fill=""))%>%
tibble::rownames_to_column()%>%
mutate(Effect=gsub(x=Effect,pattern="_",replacement=" "),
       Effect=gsub(x=Effect,pattern=":",replacement=" x "),
       Effect=gsub(x=Effect,pattern="age group",replacement="Age Group"),
       Effect=gsub(x=Effect,pattern="brain region",replacement="Brain Region"),
       Effect=gsub(x=Effect,pattern="cand RDM names",replacement="Format"))%>%
mutate_if(is.numeric,round,digits=3)%>%
mutate_at(vars(matches("p|p[GG]|p[HF]")),
          funs(gsub(x=.,pattern="^0$",replacement="<.001")))%>%
# mutate_at(vars(matches("p|p[GG]|p[HF]")),
#           funs(gsub(x=.,pattern="^0(?=\\.)",replacement="",perl = T)))%>%
write.csv(file = file.path(PATH_RSA_ANOVA_RESULT,"format_models_anova_table.csv"))
#Pairwise comparison
write.csv(x=result_pairwise_format$p.value,
          file = file.path(PATH_RSA_ANOVA_RESULT,"format_pairwise.csv"))

#Magnitude models - abs
do.call("smartbind",c(result_anova_mag_abs,fill=""))%>%
tibble::rownames_to_column()%>%
mutate(Effect=gsub(x=Effect,pattern="_",replacement=" "),
       Effect=gsub(x=Effect,pattern=":",replacement=" x "),
       Effect=gsub(x=Effect,pattern="age group",replacement="Age Group"),
       Effect=gsub(x=Effect,pattern="brain region",replacement="Brain Region"),
       Effect=gsub(x=Effect,pattern="cand RDM names",replacement="Log"))%>%
mutate_if(is.numeric,round,digits=3)%>%
mutate_at(vars(matches("p|p[GG]|p[HF]")),
          funs(gsub(x=.,pattern="^0$",replacement="<.001")))%>%
write.csv(file = file.path(PATH_RSA_ANOVA_RESULT,"mag_abs_models_anova_table.csv"))
#Magnitude models - sign
do.call("smartbind",c(result_anova_mag_sign,fill=""))%>%
tibble::rownames_to_column()%>%
mutate(Effect=gsub(x=Effect,pattern="_",replacement=" "),
       Effect=gsub(x=Effect,pattern=":",replacement=" x "),
       Effect=gsub(x=Effect,pattern="age group",replacement="Age Group"),
       Effect=gsub(x=Effect,pattern="brain region",replacement="Brain Region"),
       Effect=gsub(x=Effect,pattern="cand RDM names",replacement="Log"))%>%
mutate_if(is.numeric,round,digits=3)%>%
mutate_at(vars(matches("p|p[GG]|p[HF]")),
          funs(gsub(x=.,pattern="^0$",replacement="<.001")))%>%
write.csv(file = file.path(PATH_RSA_ANOVA_RESULT,"mag_sign_models_anova_table.csv"))
#Z-transformed tau: Split into format and magnitude models==============================================)=========================
tau_type="z_tau"
PATH_RSA_ANOVA_RESULT=PATH_RSA_ANOVA_RESULT_PARAM_Z_TAU
dir.create(PATH_RSA_ANOVA_RESULT,recursive = T)
#Format models------------------------------------------)-------------------------------------------------
#Format models (3 way: Age group x Brain Region x Model)-----------------------------------------------------------------------
df_format_models=
  df_collect_subject_cand_r%>%
  filter(model_type=="format_base",
         cand_RDM_names %in% c("g","n","null"))
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
df_mag_abs_models=
  df_collect_subject_cand_r%>%
  filter(model_type=="mag_base",
         cand_RDM_names %in% c("agv2","a"))
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
#Signed magnitude models------------------------------------------)-------------------------------------------------
#Magnitude models - Signed magnitude-----------------------------------------------------------------------
df_mag_sign_models=
  df_collect_subject_cand_r%>%
  filter(model_type=="mag_base",
         cand_RDM_names %in% c("sgv2","s"))
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
#Print out results as csv================================================================)----------------------------------------------------------------------------------
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
  write.csv(file = file.path(PATH_RSA_ANOVA_RESULT,"Format_Way_3_R~agexbrainxmodel_anova_table.csv"))
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
  write.csv(file = file.path(PATH_RSA_ANOVA_RESULT,"Format_Collapse_age_group_Way_2_R~formatxbrainxmodel_anova_table.csv"))
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
    write.csv(file = file.path(PATH_RSA_ANOVA_RESULT,
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
    write.csv(file = file.path(PATH_RSA_ANOVA_RESULT,
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
    write.csv(file = file.path(PATH_RSA_ANOVA_RESULT,
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
  write.csv(file = file.path(PATH_RSA_ANOVA_RESULT,"Abs_mag_Way_3_R~agexbrainxmodel_anova_table.csv"))
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
  write.csv(file = file.path(PATH_RSA_ANOVA_RESULT,"Abs_mag_Collapse_age_group_Way_2_R~formatxbrainxmodel_anova_table.csv"))
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
    write.csv(file = file.path(PATH_RSA_ANOVA_RESULT,
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
    write.csv(file = file.path(PATH_RSA_ANOVA_RESULT,
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
    write.csv(file = file.path(PATH_RSA_ANOVA_RESULT,
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
  write.csv(file = file.path(PATH_RSA_ANOVA_RESULT,"mag_abs_models_anova_table.csv"))
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
  write.csv(file = file.path(PATH_RSA_ANOVA_RESULT,"mag_sign_models_anova_table.csv"))

#Plot out the z_tau ~ model x age group x brain region===================================)===========================
#Constant--------------------------------
LIST_BASE_MODEL_TYPE=c("format_base","mag_base","all_base")
LIST_AGE_GROUP_RAW_NAMES=c("x2_grade","x5_grade","adult","all")
LIST_AGE_GROUP_NEW_NAMES=c("Second Graders","Fifth Graders","Adults","All")
LIST_VOI_NAMES=c("R_IPS","L_IPS","R_V1","L_V1")
# LIST_ORDERED_MAG_CAND_MODELS=c("log-A","linear-A","log-S","linear-S")
# LIST_ORDERED_FORMAT_CAND_MODELS=c("G-format","N-format","Null")
LIST_ORDERED_RDM_RAW_NAMES=c("g","n","null","a","agv2","s","sgv2")
LIST_ORDERED_RDM_NEW_NAMES=c("G-format","N-format","Null",
                             "linear-A","log-A","linear-S","log-S")

#plot--------------------------------------------------
#Formtat models
df_format_models%>%
  mutate(age_group = factor(age_group ,levels=LIST_AGE_GROUP_RAW_NAMES,labels=LIST_AGE_GROUP_NEW_NAMES),
         brain_region = factor(brain_region,levels=LIST_VOI_NAMES),
         cand_RDM_names=factor(cand_RDM_names,
                               levels=LIST_ORDERED_RDM_RAW_NAMES,
                               labels=LIST_ORDERED_RDM_NEW_NAMES))%>%
  group_by(age_group,brain_region,cand_RDM_names)%>%
  summarise(mean_z_tau=mean(z_tau),
            n=n(),
            se_z_tau=sd(z_tau)/sqrt(n))%>%
  ggplot(aes(x=cand_RDM_names,y=mean_z_tau))+
  geom_col(fill="#0033cc")+
  geom_errorbar(aes(ymax=mean_z_tau+se_z_tau,
                    ymin=mean_z_tau-se_z_tau),width=0.5)+
  facet_grid(brain_region~age_group)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5,
                                   face = "bold"))
#Magnitude models
df_collect_subject_cand_r%>%
  filter(model_type=="mag_base",
         cand_RDM_names %in% c("agv2","a","sgv2","s"))%>%
  mutate(age_group = factor(age_group ,levels=LIST_AGE_GROUP_RAW_NAMES,labels=LIST_AGE_GROUP_NEW_NAMES),
         brain_region = factor(brain_region,levels=LIST_VOI_NAMES),
         cand_RDM_names=factor(cand_RDM_names,
                               levels=LIST_ORDERED_RDM_RAW_NAMES,
                               labels=LIST_ORDERED_RDM_NEW_NAMES))%>%
  group_by(age_group,brain_region,cand_RDM_names)%>%
  summarise(mean_z_tau=mean(z_tau),
            n=n(),
            se_z_tau=sd(z_tau)/sqrt(n))%>%
  ggplot(aes(x=cand_RDM_names,y=mean_z_tau))+
  geom_col(aes(fill=cand_RDM_names))+
  scale_fill_manual(values=c("#ffd11a","#ffd11a","#cc2900","#cc2900"))+
  geom_errorbar(aes(ymax=mean_z_tau+se_z_tau,
                    ymin=mean_z_tau-se_z_tau),width=0.5)+
  facet_grid(brain_region~age_group)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5,
                                   face = "bold"))