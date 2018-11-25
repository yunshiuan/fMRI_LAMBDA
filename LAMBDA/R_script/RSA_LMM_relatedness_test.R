#RSA: LMM for relatedness test results[for base models]===========================================
#[Not used: no enough observations for the full model]
library(dplyr)
library(tidyr)
library(R.matlab)
library(stringr)
library(lme4)
library(lmerTest)
library(influence.ME)
#Constants-------------------------------------------------------------------------
#Parameter
NUM_TRIAL=16
NUN_LAYERS=3 #Base model type, age group, brain region
#Path
PATH_RSA_ROOT="D:\\Yun-Shiuan_LAMBDA\\RSA"
PATH_RSA_CURRENT_TRIAL=file.path(PATH_RSA_ROOT,"trial_17")
PATH_RSA_ANOVA_RESULT=file.path(PATH_RSA_CURRENT_TRIAL,
                                "relatedness_test_LMM","base_models")
dir.create(PATH_RSA_ANOVA_RESULT)
#File
FILE_RSA_RESULT_BASE_MODEL_MAT=file.path(PATH_RSA_ROOT,
                                         "trial_14","Part3_2nd_order_analysis","base_models",
                                         "without_ACC","relatedness_test",
                                         "Relatedness_Test_Results_14_30-Apr-2018.mat")
FILE_RELATEDNESS_TEST_HELPER_FUNCTIONS="D:\\GoogleDrive\\Lambda_code\\R_script\\R_functions\\RSA_relatedness_test_visualization_functions.R"
#Helper functions---------------------------------------------------------------
source(FILE_RELATEDNESS_TEST_HELPER_FUNCTIONS)

#Read in relatedness results===========================)================================================
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
#Split into format and magnitude models==============================================)=========================
#Format models-----------------------------------------------------------------------
df_format_models=
  df_collect_subject_cand_r%>%
  filter(model_type=="format_base",
         cand_RDM_names %in% c("g","n","null"))
#LMM
m=lmer(cand_r~cand_RDM_names*age_group*brain_region+(1+cand_RDM_names+brain_region|subject_num),
       data = df_format_models, control=lmerControl(optCtrl = list(maxfun = 1e5)))
a=anova
infl=influence(m, obs = TRUE)
cooks.distance(infl)
