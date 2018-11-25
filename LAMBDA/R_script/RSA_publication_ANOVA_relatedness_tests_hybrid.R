#RSA: ANOVA for relatedness test results[for hybrid models]===========================================
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
#(2018/6/11)
#Constants-------------------------------------------------------------------------
#Parameter
NUN_LAYERS=2 #age group, brain region
BOOLEAN_PRINT_PARAMETRIC_ANOVA_CSV=F
BOOLEAN_PRINT_ROBUST_ANOVA_CSV=T
LIST_NAMES_BRAIN_REGION=c("R_IPS","L_IPS","R_V1","L_V1")
LIST_NAMES_AGE_GROUP=c("x2_grade","x5_grade","adult")
LIST_NAMES_TAU_TYPE=c("cand_r","z_tau")

#Path
PATH_RSA_ROOT="D:\\Yun-Shiuan_LAMBDA\\RSA"
PATH_RSA_CURRENT_TRIAL=file.path(PATH_RSA_ROOT,"trial_22")
PATH_RSA_ANOVA_RESULT_PARAM_RAW_TAU=file.path(PATH_RSA_CURRENT_TRIAL,"Tables",
                                              "relatedness_test_Anova","hybrid_models","parametric","raw_tau")
PATH_RSA_ANOVA_RESULT_ROBUST_RAW_TAU=file.path(PATH_RSA_CURRENT_TRIAL,"Tables",
                                               "relatedness_test_Anova","hybrid_models","robust","raw_tau")
PATH_RSA_ANOVA_RESULT_PARAM_Z_TAU=file.path(PATH_RSA_CURRENT_TRIAL,"Tables",
                                            "relatedness_test_Anova","hybrid_models","parametric","z_tau")
PATH_RSA_ANOVA_RESULT_ROBUST_Z_TAU=file.path(PATH_RSA_CURRENT_TRIAL,"Tables",
                                             "relatedness_test_Anova","hybrid_models","robust","z_tau")
dir.create(PATH_RSA_ANOVA_RESULT_PARAM_RAW_TAU,recursive = T)
dir.create(PATH_RSA_ANOVA_RESULT_PARAM_Z_TAU,recursive = T)

dir.create(PATH_RSA_ANOVA_RESULT_ROBUST_RAW_TAU,recursive = T)
dir.create(PATH_RSA_ANOVA_RESULT_ROBUST_Z_TAU,recursive = T)

#File
FILE_RSA_RESULT_HYBRID_MODEL_MAT=file.path(PATH_RSA_ROOT,
                                           "trial_22","Part3_2nd_order_analysis","hybrid_models",
                                           "without_ACC","relatedness_test",
                                           "Relatedness_Test_Results_22_10-Jun-2018.mat")
FILE_RELATEDNESS_TEST_HELPER_FUNCTIONS="D:\\GoogleDrive\\Lambda_code\\R_script\\R_functions\\RSA_relatedness_test_visualization_functions.R"
#Helper functions---------------------------------------------------------------
source(FILE_RELATEDNESS_TEST_HELPER_FUNCTIONS)

#Read in relatedness results===========================)================================================
#Read in the mat file and flatten it into a list-------------------------------
list_rsa_result=RSA_relatedness_nested_mat_to_single_mats(FILE_RSA_RESULT_HYBRID_MODEL_MAT,NUN_LAYERS)

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
                                     mutate(age_group=str_extract(string=df_name,
                                                                  pattern=".*(?=_R_|_L_)"),
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
#Identify the interested models for ANOVA, 1)best GxA model, 2)best GxS model, 3)S model, 4)A model, and 5)G model---------------
#Get the name of the best GxA/GxS models for each age group and brain region
df_best_hybrid_model=
  df_collect_subject_cand_r%>%
  group_by(age_group,brain_region,cand_RDM_names)%>%
  summarise(mean_cand_r=mean(cand_r))%>%
  arrange(desc(mean_cand_r))%>%
  slice(grep(x=cand_RDM_names,
             pattern="(R_|L_|^g$|^agv2$|sgv2$)",
             invert = T))%>%
  mutate(hybrid_type=ifelse(grepl(x=cand_RDM_names,
                                  pattern="gagv"),
                            yes="GxA",
                            no="GxS"))%>%
  group_by(age_group,brain_region,hybrid_type)%>%
  top_n(n = 1,wt = mean_cand_r)%>%
  #Note that there are some tying models, but it's okay to select a random one.
  #This is because the following ANOVA only concern about whether the hybrid one
  #outperforms the base one.
  slice(1) 

#Filter in the interested models to compare---------------------------------------------------
df_anova=
  df_collect_subject_cand_r%>%
  left_join(df_best_hybrid_model%>%
              select(-mean_cand_r),
            by=c("age_group","brain_region","cand_RDM_names"))%>%
  filter(grepl(x=cand_RDM_names,
               pattern="^g$|^agv2$|sgv2$")|!is.na(hybrid_type))%>%
  mutate(cand_RDM_names=as.character(cand_RDM_names),
         cand_RDM_names=ifelse(grepl(x=cand_RDM_names,
                                     pattern="gagv|gsgv"),
                               yes = hybrid_type,
                               no=cand_RDM_names))
#Z-transform per age group x brain region--------------------------------------------
df_anova=
  df_anova%>%
  group_by(age_group,brain_region)%>%
  mutate(z_tau=scale(cand_r))%>%
  as.data.frame()
#Raw tau & z_tau: Parametric Anova============================================)======================================================
for (tau_type_index in 1:length(LIST_NAMES_TAU_TYPE)){
  #Full model (3-way: Age group x Brain Region x Model: including IPS and V1)-------------------------------
  tau_type=LIST_NAMES_TAU_TYPE[tau_type_index]
  result_anova_full_model=  
    ezANOVA(
      data = df_anova%>%
        rename(dv=tau_type),
      dv = .(dv),
      wid = .(subject_num),
      type = 3,
      within = .(brain_region,cand_RDM_names),
      between = .(age_group),
      detailed=T)
  #When 3-way interaction is not significant: Collpase the age group------------------------------)--------------------------------
  #Collpase the age group (2-way: Brain x Model)-------------------------
  result_anova_2way_collapse_age_group=
    ezANOVA(
      data = df_anova%>%
        rename(dv=tau_type),
      dv = .(dv),
      wid = .(subject_num),
      type = 3,
      within = .(brain_region,cand_RDM_names),
      detailed=T)
  #Collpase the age group: Split by brian region (1-way: Model)-----------------------------------
  list_anova_result_collapse_age_group_split_brain_region=c()
  for (brain_region_index in 1:length(LIST_NAMES_BRAIN_REGION)){
    name_brain_region=LIST_NAMES_BRAIN_REGION[brain_region_index]
    result_anova=
      ezANOVA(
        data = df_anova%>%
          filter(brain_region==name_brain_region)%>%
          rename(dv=tau_type),
        dv = .(dv),
        wid = .(subject_num),
        type = 3,
        within = .(cand_RDM_names),
        detailed=T)
    list_anova_result_collapse_age_group_split_brain_region[[brain_region_index]]=
      result_anova
  }
  names(list_anova_result_collapse_age_group_split_brain_region)=LIST_NAMES_BRAIN_REGION
  #Collpase the age group: Pairwise comparison for model effect--------------------------------------------------------------------------
  #All pairs-----------------------------------------------
  list_pariwise_comparison_collapse_age_group_split_brain_region=c()
  
  for (brain_region_index in 1:length(LIST_NAMES_BRAIN_REGION)){
    name_brain_region=LIST_NAMES_BRAIN_REGION[brain_region_index]
    list_pariwise_comparison_split_brain_region=c()#Collect the age group results per brain region
    
    #Order the model by the cand_r value
    df_model_ordered_level=
      df_anova%>%
      filter(brain_region==name_brain_region)%>%
      group_by(cand_RDM_names)%>%
      rename(dv=tau_type)%>%
      summarise(mean_dv=mean(dv))%>%
      arrange(desc(mean_dv))%>%
      mutate(cand_full_RDM_names=
               case_when(cand_RDM_names=="g" ~ "G-Format",
                         cand_RDM_names=="agv2" ~ "log-A",
                         cand_RDM_names=="sgv2" ~ "log-S",
                         TRUE~cand_RDM_names))
    result_pairwise_format=
      pairwise.t.test(df_anova%>%
                        rename(dv=tau_type)%>%
                        filter(brain_region==name_brain_region)%>%
                        pull(dv),
                      df_anova%>%
                        filter(brain_region==name_brain_region)%>%
                        pull(cand_RDM_names)%>%
                        factor(levels=df_model_ordered_level$cand_RDM_names,
                               labels = df_model_ordered_level$cand_full_RDM_names),
                      paired=T, p.adjust.method="bonferroni")
    
    list_pariwise_comparison_collapse_age_group_split_brain_region[[brain_region_index]]=
      result_pairwise_format
  }
  names(list_pariwise_comparison_collapse_age_group_split_brain_region)=LIST_NAMES_BRAIN_REGION
  
  #Only focus on the best model (compare it to the others)--------------------------------------
  list_best_model_comparison_collapse_age_group_split_brain_region=c()
  
  for (brain_region_index in 1:length(LIST_NAMES_BRAIN_REGION)){
    name_brain_region=LIST_NAMES_BRAIN_REGION[brain_region_index]
    #Order the model by the cand_r value
    df_model_ordered_level=
      df_anova%>%
      filter(brain_region==name_brain_region)%>%
      group_by(cand_RDM_names)%>%
      rename(dv=tau_type)%>%
      summarise(mean_dv=mean(dv))%>%
      arrange(desc(mean_dv))%>%
      mutate(cand_full_RDM_names=
               case_when(cand_RDM_names=="g" ~ "G-Format",
                         cand_RDM_names=="agv2" ~ "log-A",
                         cand_RDM_names=="sgv2" ~ "log-S",
                         TRUE~cand_RDM_names))
    result_pairwise_format=
      pairwise.t.test(df_anova%>%
                        filter(brain_region==name_brain_region)%>%
                        rename(dv=tau_type)%>%
                        pull(dv),
                      df_anova%>%
                        filter(brain_region==name_brain_region)%>%
                        pull(cand_RDM_names)%>%
                        factor(levels=df_model_ordered_level$cand_RDM_names,
                               labels = df_model_ordered_level$cand_full_RDM_names),
                      paired=T, p.adjust.method="none",alternative = "less")#Remember to correct the p value afterwards
    list_best_model_comparison_collapse_age_group_split_brain_region[[brain_region_index]]=
      result_pairwise_format
  }
  names(list_best_model_comparison_collapse_age_group_split_brain_region)=LIST_NAMES_BRAIN_REGION
  
  #When 3-way interaction is significant-----------------------------------)---------------------------
  #Split by brain region (2-way: Age group x Model)------------------------------------------------------------
  list_2way_anova_result_split_brain_region=c()
  
  for (brain_region_index in 1:length(LIST_NAMES_BRAIN_REGION)){
    name_brain_region=LIST_NAMES_BRAIN_REGION[brain_region_index]
    result_anova=
      ezANOVA(
        data = df_anova%>%
          filter(brain_region==name_brain_region)%>%
          rename(dv=tau_type),
        dv = .(dv),
        wid = .(subject_num),
        type = 3,
        within = .(cand_RDM_names),
        between = .(age_group),
        detailed=T)
    list_2way_anova_result_split_brain_region[[brain_region_index]]=result_anova
  }
  names(list_2way_anova_result_split_brain_region)=LIST_NAMES_BRAIN_REGION
  
  #Split by age group (2-way: Brain x Model)------------------------------------------------------------
  list_2way_anova_result_split_age_group=c()
  
  for (age_group_index in 1:length(LIST_NAMES_AGE_GROUP)){
    name_age_group=LIST_NAMES_AGE_GROUP[age_group_index]
    result_anova=
      ezANOVA(
        data = df_anova%>%
          filter(age_group==name_age_group)%>%
          rename(dv=tau_type),
        dv = .(dv),
        wid = .(subject_num),
        type = 3,
        within = .(brain_region,cand_RDM_names),
        detailed=T)
    list_2way_anova_result_split_age_group[[age_group_index]]=result_anova
  }
  names(list_2way_anova_result_split_age_group)=LIST_NAMES_AGE_GROUP
  
  
  #Split by brain region and age group(1-way: Model)------------------------------------------------------------
  list_anova_result_split_brain_region_age_group=c()
  
  for (brain_region_index in 1:length(LIST_NAMES_BRAIN_REGION)){
    list_anova_result_split_age_group=c()#Collect the age group results per brain region
    name_brain_region=LIST_NAMES_BRAIN_REGION[brain_region_index]
    
    for (age_group_index in 1:length(LIST_NAMES_AGE_GROUP)){
      name_age_group=LIST_NAMES_AGE_GROUP[age_group_index]
      result_anova=
        ezANOVA(
          data = df_anova%>%
            filter(brain_region==name_brain_region,age_group==name_age_group)%>%
            rename(dv=tau_type),
          dv = .(dv),
          wid = .(subject_num),
          type = 3,
          within = .(cand_RDM_names),
          detailed=T)
      list_anova_result_split_age_group[[age_group_index]]=result_anova
    }
    names(list_anova_result_split_age_group)=LIST_NAMES_AGE_GROUP
    list_anova_result_split_brain_region_age_group[[brain_region_index]]=list_anova_result_split_age_group
  }
  names(list_anova_result_split_brain_region_age_group)=LIST_NAMES_BRAIN_REGION
  
  #Pairwise comparison for model effect--------------------------------------------------------------------------
  #All pairs-----------------------------------------------
  list_pariwise_comparison_split_brain_region_age_group=c()
  
  for (brain_region_index in 1:length(LIST_NAMES_BRAIN_REGION)){
    name_brain_region=LIST_NAMES_BRAIN_REGION[brain_region_index]
    list_pariwise_comparison_split_brain_region=c()#Collect the age group results per brain region
    
    for (age_group_index in 1:length(LIST_NAMES_AGE_GROUP)){
      name_age_group=LIST_NAMES_AGE_GROUP[age_group_index]
      
      #Order the model by the cand_r value
      df_model_ordered_level=
        df_anova%>%
        filter(brain_region==name_brain_region,
               age_group==name_age_group)%>%
        group_by(cand_RDM_names)%>%
        rename(dv=tau_type)%>%
        summarise(mean_dv=mean(dv))%>%
        arrange(desc(mean_dv))%>%
        mutate(cand_full_RDM_names=
                 case_when(cand_RDM_names=="g" ~ "G-Format",
                           cand_RDM_names=="agv2" ~ "log-A",
                           cand_RDM_names=="sgv2" ~ "log-S",
                           TRUE~cand_RDM_names))
      result_pairwise_format=
        pairwise.t.test(df_anova%>%
                          rename(dv=tau_type)%>%
                          filter(brain_region==name_brain_region,age_group==name_age_group)%>%
                          pull(dv),
                        df_anova%>%
                          filter(brain_region==name_brain_region,age_group==name_age_group)%>%
                          pull(cand_RDM_names)%>%
                          factor(levels=df_model_ordered_level$cand_RDM_names,
                                 labels = df_model_ordered_level$cand_full_RDM_names),
                        paired=T, p.adjust.method="bonferroni")
      list_pariwise_comparison_split_brain_region[[age_group_index]]=result_pairwise_format
    }
    names(list_pariwise_comparison_split_brain_region)=LIST_NAMES_AGE_GROUP
    list_pariwise_comparison_split_brain_region_age_group[[brain_region_index]]=
      list_pariwise_comparison_split_brain_region
  }
  names(list_pariwise_comparison_split_brain_region_age_group)=LIST_NAMES_BRAIN_REGION
  #Only focus on the best model (compare it to the others)--------------------------------------
  LIST_NAMES_BRAIN_REGION=c("R_IPS","L_IPS","R_V1","L_V1")
  LIST_NAMES_AGE_GROUP=c("x2_grade","x5_grade","adult")
  list_best_model_comparison_split_brain_region_age_group=c()
  
  for (brain_region_index in 1:length(LIST_NAMES_BRAIN_REGION)){
    name_brain_region=LIST_NAMES_BRAIN_REGION[brain_region_index]
    list_pariwise_comparison_split_brain_region=c()#Collect the age group results per brain region
    
    for (age_group_index in 1:length(LIST_NAMES_AGE_GROUP)){
      name_age_group=LIST_NAMES_AGE_GROUP[age_group_index]
      
      #Order the model by the cand_r value
      df_model_ordered_level=
        df_anova%>%
        filter(brain_region==name_brain_region,age_group==name_age_group)%>%
        group_by(cand_RDM_names)%>%
        rename(dv=tau_type)%>%
        summarise(mean_dv=mean(dv))%>%
        arrange(desc(mean_dv))%>%
        mutate(cand_full_RDM_names=
                 case_when(cand_RDM_names=="g" ~ "G-Format",
                           cand_RDM_names=="agv2" ~ "log-A",
                           cand_RDM_names=="sgv2" ~ "log-S",
                           TRUE~cand_RDM_names))
      result_pairwise_format=
        pairwise.t.test(df_anova%>%
                          filter(brain_region==name_brain_region,age_group==name_age_group)%>%
                          rename(dv=tau_type)%>%
                          pull(dv),
                        df_anova%>%
                          filter(brain_region==name_brain_region,age_group==name_age_group)%>%
                          pull(cand_RDM_names)%>%
                          factor(levels=df_model_ordered_level$cand_RDM_names,
                                 labels = df_model_ordered_level$cand_full_RDM_names),
                        paired=T, p.adjust.method="none",alternative = "less")#Remember to correct the p value afterwards
      list_pariwise_comparison_split_brain_region[[age_group_index]]=result_pairwise_format
    }
    names(list_pariwise_comparison_split_brain_region)=LIST_NAMES_AGE_GROUP
    list_best_model_comparison_split_brain_region_age_group[[brain_region_index]]=
      list_pariwise_comparison_split_brain_region
  }
  names(list_best_model_comparison_split_brain_region_age_group)=LIST_NAMES_BRAIN_REGION
  
  #Print out results as csv================================)=======================================================
  if(BOOLEAN_PRINT_PARAMETRIC_ANOVA_CSV){
    if(tau_type=="cand_r"){
      PATH_RSA_ANOVA_RESULT_PARAM=PATH_RSA_ANOVA_RESULT_PARAM_RAW_TAU
    }else if(tau_type=="z_tau"){
      PATH_RSA_ANOVA_RESULT_PARAM=PATH_RSA_ANOVA_RESULT_PARAM_Z_TAU
    }
    #Full model (3-way: Age group x Brain Region x Model: including IPS and V1)------------------------------------------------------------
    do.call("smartbind",c(result_anova_full_model,fill=""))%>%
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
      write.csv(file = file.path(PATH_RSA_ANOVA_RESULT_PARAM,"Way_3_R~agexbrainxmodel_anova_table.csv"))
    
    #When 3-way interaction is not significant: Collpase the age group------------------------------)--------------------------------
    #Collpase the age group (2-way: Brain x Model)--------------------------------------------------
    anova_result=result_anova_2way_collapse_age_group
    
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
      write.csv(file = file.path(PATH_RSA_ANOVA_RESULT_PARAM,
                                 paste0("Collapse_age_group_Way_2_R~brainxmodel_anova_table_",
                                        name_age_group,".csv")))
    #Collpase the age group: Split by brian region (1-way: Model)-----------------------------------
    for (brain_region_index in 1:length(LIST_NAMES_BRAIN_REGION)){
      name_brain_region=LIST_NAMES_BRAIN_REGION[brain_region_index]
      anova_result=list_anova_result_collapse_age_group_split_brain_region[[brain_region_index]]
      
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
        write.csv(file = file.path(PATH_RSA_ANOVA_RESULT_PARAM,
                                   paste0("Collapse_age_group_Way_1_R~model_anova_table_",
                                          name_brain_region,".csv")))
      
    }
    #Collpase the age group: Pair wise comparison for model effect---------------------------------
    #All pairs-----------------------------------------------
    for (brain_region_index in 1:length(LIST_NAMES_BRAIN_REGION)){
      name_brain_region=LIST_NAMES_BRAIN_REGION[brain_region_index]
      list_pariwise_comparison_collapse_age_group_split_brain_region[[name_brain_region]]$p.value%>%
        as.data.frame()%>%
        tibble::rownames_to_column(var="Model")%>%
        mutate_if(is.numeric,round,digits=3)%>%
        mutate_if(is.numeric,gsub,pattern="^0$",replacement="<.001")%>%
        write.csv(file = file.path(PATH_RSA_ANOVA_RESULT_PARAM,
                                   paste0("Collapse_age_group_Pairwise_comparison_",
                                          name_brain_region,".csv")))
      
    }
    #Only focus on the best model (compare it to the others)----------------------------------------
    for (brain_region_index in 1:length(LIST_NAMES_BRAIN_REGION)){
      name_brain_region=LIST_NAMES_BRAIN_REGION[brain_region_index]
      list_best_model_comparison_collapse_age_group_split_brain_region[[name_brain_region]]$p.value%>%
        as.data.frame()%>%
        tibble::rownames_to_column(var="Model")%>%
        select(1:2)%>%#The first two columns (model name and the best model)
        mutate_if(is.numeric,p.adjust,"bonferroni")%>%#Correct the raw p-value
        mutate_if(is.numeric,round,digits=3)%>%
        mutate_if(is.numeric,gsub,pattern="^0$",replacement="<.001")%>%
        write.csv(file = file.path(PATH_RSA_ANOVA_RESULT_PARAM,
                                   paste0("Collapse_age_group_Pair_to_best_model_",
                                          name_brain_region,".csv")))
    }        
    #When 3-way interaction is significant-----------------------------------)---------------------------
    #Split by brain region (2-way: Age group x Model)------------------------
    for (brain_region_index in 1:length(LIST_NAMES_BRAIN_REGION)){
      name_brain_region=LIST_NAMES_BRAIN_REGION[brain_region_index]
      anova_result=list_2way_anova_result_split_brain_region[[brain_region_index]]
      
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
        write.csv(file = file.path(PATH_RSA_ANOVA_RESULT_PARAM,
                                   paste0("Way_2_R~agexmodel_anova_table_",
                                          name_brain_region,".csv")))
    }
    #Split by age group (2-way: Brain x Model)------------------------------------------------------------
    for (age_group_index in 1:length(LIST_NAMES_AGE_GROUP)){
      name_age_group=LIST_NAMES_AGE_GROUP[age_group_index]
      anova_result=list_2way_anova_result_split_age_group[[age_group_index]]
      
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
        write.csv(file = file.path(PATH_RSA_ANOVA_RESULT_PARAM,
                                   paste0("Way_2_R~brainxmodel_anova_table_",
                                          name_age_group,".csv")))
    }
    #Split by brain region and age group (1-way anova)-------------------------------------
    for (brain_region_index in 1:length(LIST_NAMES_BRAIN_REGION)){
      name_brain_region=LIST_NAMES_BRAIN_REGION[brain_region_index]
      
      for (age_group_index in 1:length(LIST_NAMES_AGE_GROUP)){
        name_age_group=LIST_NAMES_AGE_GROUP[age_group_index]
        
        anova_result=list_anova_result_split_brain_region_age_group[[brain_region_index]][[age_group_index]]
        
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
          write.csv(file = file.path(PATH_RSA_ANOVA_RESULT_PARAM,
                                     paste0("Way_1_R~model_anova_table_",
                                            name_brain_region,"_",
                                            name_age_group,".csv")))
      }
    }
    #Pair wise comparison for model effect--------------------------------------------------------------------------
    #All pairs-----------------------------------------------
    for (brain_region_index in 1:length(LIST_NAMES_BRAIN_REGION)){
      name_brain_region=LIST_NAMES_BRAIN_REGION[brain_region_index]
      
      for (age_group_index in 1:length(LIST_NAMES_AGE_GROUP)){
        name_age_group=LIST_NAMES_AGE_GROUP[age_group_index]
        
        list_pariwise_comparison_split_brain_region_age_group[[name_brain_region]][[name_age_group]]$p.value%>%
          as.data.frame()%>%
          tibble::rownames_to_column(var="Model")%>%
          mutate_if(is.numeric,round,digits=3)%>%
          mutate_if(is.numeric,gsub,pattern="^0$",replacement="<.001")%>%
          write.csv(file = file.path(PATH_RSA_ANOVA_RESULT_PARAM,
                                     paste0("Pairwise_comparison_",
                                            name_brain_region,"_",
                                            name_age_group,".csv")))
      }
    }
    #Only focus on the best model (compare it to the others)-----------------------------------------------
    for (brain_region_index in 1:length(LIST_NAMES_BRAIN_REGION)){
      name_brain_region=LIST_NAMES_BRAIN_REGION[brain_region_index]
      
      for (age_group_index in 1:length(LIST_NAMES_AGE_GROUP)){
        name_age_group=LIST_NAMES_AGE_GROUP[age_group_index]
        
        list_best_model_comparison_split_brain_region_age_group[[name_brain_region]][[name_age_group]]$p.value%>%
          as.data.frame()%>%
          tibble::rownames_to_column(var="Model")%>%
          select(1:2)%>%#The first two columns (model name and the best model)
          mutate_if(is.numeric,p.adjust,"bonferroni")%>%#Correct the raw p-value
          mutate_if(is.numeric,round,digits=3)%>%
          mutate_if(is.numeric,gsub,pattern="^0$",replacement="<.001")%>%
          write.csv(file = file.path(PATH_RSA_ANOVA_RESULT_PARAM,
                                     paste0("Pair_to_best_model_",
                                            name_brain_region,"_",
                                            name_age_group,".csv")))
      }
    }
  }
}
#Raw Tau & z_tau: Robust Anova==================================================)========================
for (tau_type_index in 1:length(LIST_NAMES_TAU_TYPE)){
  tau_type=LIST_NAMES_TAU_TYPE[tau_type_index]
  #When 3-way interaction is not significant: Collpase the age group------------------------------)--------------------------------
  #Collapse age group: Split by brain region: RDM correlation ~ Model (1w: trimmed)----------------------------------------------
  list_anova_result_collapse_age_group_split_brain_region_trimmed=c()
  
  for (brain_region_index in 1:length(LIST_NAMES_BRAIN_REGION)){
    list_anova_result_split_age_group=c() #Collect the age group results per brain region
    name_brain_region=LIST_NAMES_BRAIN_REGION[brain_region_index]
    df=
      df_anova%>%
      filter(brain_region==name_brain_region)%>%
      rename(dv=tau_type)
    result_anova_trimmed=
      rmanova(y = df$dv,
              groups = df$cand_RDM_names,
              blocks = df$subject_num)
    
    list_anova_result_collapse_age_group_split_brain_region_trimmed[[brain_region_index]]=result_anova_trimmed
  }
  names(list_anova_result_collapse_age_group_split_brain_region_trimmed)=LIST_NAMES_BRAIN_REGION
  #Collapse age group: Post-hoc: RDM correlation ~ Model (post-hoc: trimmed)--------------------------------------------
  #Collapse age group: Pairwise: Only focus on the best model (compare it to the others)--------------------------------------
  list_best_model_comparison_collapse_age_group_split_brain_region_trimmed=c()
  
  #The p criteria for Rom's method (used by the package and recommended by Wilcox, R. (2012))
  list_rom_correction=c(0.05000,0.02500,0.01690,0.01270,0.01020,
                        0.00851,0.00730,0.00639,0.00568,0.00511)
  list_ordered_RDM=c("GxS","G-Format","log-A","log-S")
  for (brain_region_index in 1:length(LIST_NAMES_BRAIN_REGION)){
    list_anova_result_split_age_group=c() #Collect the age group results per brain region
    name_brain_region=LIST_NAMES_BRAIN_REGION[brain_region_index]
    
    df=
      df_anova%>%
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
      mutate(p.crit=list_rom_correction[1:4])%>%
      arrange(Group2)%>%
      mutate(Group1=result_anova$fnames[1],
             Group2=result_anova$fnames[2:length(result_anova$fnames)],
             Significant_Rom=ifelse(p.value<p.crit,yes = "*",no = ""),
             Group1=case_when(Group1=="g" ~ "G-Format",
                              Group1=="agv2" ~ "log-A",
                              Group1=="sgv2" ~ "log-S",
                              TRUE~Group1
             ),
             Group2=case_when(Group2=="g" ~ "G-Format",
                              Group2=="agv2" ~ "log-A",
                              Group2=="sgv2" ~ "log-S",
                              TRUE~Group2))%>%
      rename(p.value.raw=p.value)%>%
      mutate(p.bonferroni=p.value.raw*(length(result_anova$fnames)-1),
             p.bonferroni=ifelse(p.bonferroni>1,yes = 1,no = p.bonferroni))%>%
      mutate_if(is.numeric,round,3)%>%
      mutate_at(vars(matches("^p\\..*")),
                funs(gsub(x=.,pattern="^0$",replacement="<.001")))%>%
      mutate(Group2=factor(Group2,levels = list_ordered_RDM))%>%
      arrange(Group2)
    list_best_model_comparison_collapse_age_group_split_brain_region_trimmed[[brain_region_index]]=result_anova
  }
  names(list_best_model_comparison_collapse_age_group_split_brain_region_trimmed)=LIST_NAMES_BRAIN_REGION
  
  #When 3-way interaction is significant-----------------------------------------------------------)---------------------------
  #Note: This only works for at most: 1b1w, 1w, or 3b.
  #Split by the brain region: RDM correlation ~ Age group x Model (1b1w: trimmed)------------------------------
  list_2way_anova_result_split_brain_region_trimmed=c()
  # list_2way_anova_result_split_brain_region_boot=c()
  for (brain_region_index in 1:length(LIST_NAMES_BRAIN_REGION)){
    name_brain_region=LIST_NAMES_BRAIN_REGION[brain_region_index]
    #Filter in the brain region
    df=
      df_anova%>%
      filter(brain_region==name_brain_region)%>%
      rename(dv=tau_type)
    #Formula
    formula=formula("cand_r~age_group*cand_RDM_names")
    #Trimmed result
    result_anov_trimmed=
      bwtrim(formula=dv~age_group*cand_RDM_names,
             id = subject_num,data = df)
    # #Bootstrap (not working for the sppbb and sppbi)
    # age_group_boot=
    #   sppba(formula=formula,
    #         id = subject_num,data = df)
    # cand_RDM_boot=
    #   sppbb(formula=formula,
    #         id = subject_num,data = df)
    # agexcand_boot=
    #   sppbi(formula=formula,
    #         id = subject_num,data = df)
    list_2way_anova_result_split_brain_region_trimmed[[brain_region_index]]=result_anov_trimmed
  }
  names(list_2way_anova_result_split_brain_region_trimmed)=LIST_NAMES_BRAIN_REGION
  
  #Further split by age group: RDM correlation ~ Model (1w: trimmed)----------------------------------------------
  list_anova_result_split_brain_region_age_group_trimmed=c()
  
  for (brain_region_index in 1:length(LIST_NAMES_BRAIN_REGION)){
    list_anova_result_split_age_group=c() #Collect the age group results per brain region
    name_brain_region=LIST_NAMES_BRAIN_REGION[brain_region_index]
    
    for (age_group_index in 1:length(LIST_NAMES_AGE_GROUP)){
      name_age_group=LIST_NAMES_AGE_GROUP[age_group_index]
      df=
        df_anova%>%
        filter(age_group==name_age_group,
               brain_region==name_brain_region)%>%
        rename(dv=tau_type)
      result_anov_trimmed=
        rmanova(y = df$dv,
                groups = df$cand_RDM_names,
                blocks = df$subject_num)
      list_anova_result_split_age_group[[age_group_index]]=result_anov_trimmed
    }
    names(list_anova_result_split_age_group)=LIST_NAMES_AGE_GROUP
    list_anova_result_split_brain_region_age_group_trimmed[[brain_region_index]]=list_anova_result_split_age_group
  }
  names(list_anova_result_split_brain_region_age_group_trimmed)=LIST_NAMES_BRAIN_REGION
  
  #Post-hoc: RDM correlation ~ Model (post-hoc: trimmed)--------------------------------------------
  #Pairwise: Only focus on the best model (compare it to the others)--------------------------------------
  list_best_model_comparison_split_brain_region_age_group_trimmed=c()
  
  #The p criteria for Rom's method (used by the package and recommended by Wilcox, R. (2012))
  list_rom_correction=c(0.05000,0.02500,0.01690,0.01270,0.01020,
                        0.00851,0.00730,0.00639,0.00568,0.00511)
  list_ordered_RDM=c("GxS","G-Format","log-A","log-S")
  for (brain_region_index in 1:length(LIST_NAMES_BRAIN_REGION)){
    list_anova_result_split_age_group=c() #Collect the age group results per brain region
    name_brain_region=LIST_NAMES_BRAIN_REGION[brain_region_index]
    
    for (age_group_index in 1:length(LIST_NAMES_AGE_GROUP)){
      name_age_group=LIST_NAMES_AGE_GROUP[age_group_index]
      df=
        df_anova%>%
        filter(age_group==name_age_group,
               brain_region==name_brain_region)%>%
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
        mutate(p.crit=list_rom_correction[1:4])%>%
        arrange(Group2)%>%
        mutate(Group1=result_anova$fnames[1],
               Group2=result_anova$fnames[2:length(result_anova$fnames)],
               Significant_Rom=ifelse(p.value<p.crit,yes = "*",no = ""),
               Group1=case_when(Group1=="g" ~ "G-Format",
                                Group1=="agv2" ~ "log-A",
                                Group1=="sgv2" ~ "log-S",
                                TRUE~Group1
               ),
               Group2=case_when(Group2=="g" ~ "G-Format",
                                Group2=="agv2" ~ "log-A",
                                Group2=="sgv2" ~ "log-S",
                                TRUE~Group2))%>%
        rename(p.value.raw=p.value)%>%
        mutate(p.bonferroni=p.value.raw*(length(result_anova$fnames)-1),
               p.bonferroni=ifelse(p.bonferroni>1,yes = 1,no = p.bonferroni))%>%
        mutate_if(is.numeric,round,3)%>%
        mutate_at(vars(matches("^p\\..*")),
                  funs(gsub(x=.,pattern="^0$",replacement="<.001")))%>%
        mutate(Group2=factor(Group2,levels = list_ordered_RDM))%>%
        arrange(Group2)
      
      
      list_anova_result_split_age_group[[age_group_index]]=result_anova
    }
    names(list_anova_result_split_age_group)=LIST_NAMES_AGE_GROUP
    list_best_model_comparison_split_brain_region_age_group_trimmed[[brain_region_index]]=list_anova_result_split_age_group
  }
  names(list_best_model_comparison_split_brain_region_age_group_trimmed)=LIST_NAMES_BRAIN_REGION
  
  
  #Print out results as csv================================)=======================================================
  if(BOOLEAN_PRINT_ROBUST_ANOVA_CSV){
    if(tau_type=="cand_r"){
      PATH_RSA_ANOVA_RESULT_ROBUST=PATH_RSA_ANOVA_RESULT_ROBUST_RAW_TAU
    }else if(tau_type=="z_tau"){
      PATH_RSA_ANOVA_RESULT_ROBUST=PATH_RSA_ANOVA_RESULT_ROBUST_Z_TAU
    }
    #When 3-way interaction is not significant: Collpase the age group------------------------------)--------------------------------
    #Collapse age group: Split by brain region: RDM correlation ~ Model (1w: trimmed)-------------------------------------
    for (brain_region_index in 1:length(LIST_NAMES_BRAIN_REGION)){
      name_brain_region=LIST_NAMES_BRAIN_REGION[brain_region_index]
      anova_result=list_anova_result_collapse_age_group_split_brain_region_trimmed[[brain_region_index]]
      data.frame(Effect="Model",
                 DFn=anova_result$df1,
                 DFd=anova_result$df2,
                 `F`=anova_result$test,
                 p=anova_result$p.value)%>%
        mutate_if(is.numeric,round,3)%>%
        mutate_at(vars(matches("p")),
                  funs(gsub(x=.,pattern="^0$",replacement="<.001")))%>%
        write.csv(file = file.path(PATH_RSA_ANOVA_RESULT_ROBUST,
                                   paste0("Way_1_R~model_anova_table_",
                                          name_brain_region,".csv")))
    }
    #Collapse age group: Pair wise comparison for model effect--------------------------------------------------------------------------
    #Collapse age group: Only focus on the best model (compare it to the others)-----------------------------------------------
    for (brain_region_index in 1:length(LIST_NAMES_BRAIN_REGION)){
      name_brain_region=LIST_NAMES_BRAIN_REGION[brain_region_index]
      anova_result=list_best_model_comparison_collapse_age_group_split_brain_region_trimmed[[name_brain_region]]
      anova_result$comp%>%
        write.csv(file = file.path(PATH_RSA_ANOVA_RESULT_ROBUST,
                                   paste0("Pair_to_best_model_",
                                          name_brain_region,".csv")))
      
    }
    #When 3-way interaction is significant-----------------------------------------------------------)---------------------------
    #Note: This only works for at most: 1b1w, 1w, or 3b.
    #Split by the brain region: RDM correlation ~ Age group x Model (1b1w: trimmed)------------------------
    for (brain_region_index in 1:length(LIST_NAMES_BRAIN_REGION)){
      name_brain_region=LIST_NAMES_BRAIN_REGION[brain_region_index]
      anova_result=list_2way_anova_result_split_brain_region_trimmed[[brain_region_index]]
      
      df=data.frame(Effect=c(anova_result$varnames[2],
                             anova_result$varnames[3],
                             paste0(anova_result$varnames[2]," x ",anova_result$varnames[3])),
                    `F`=c(anova_result$Qa,anova_result$Qb,anova_result$Qab),
                    p=c(anova_result$A.p.value,anova_result$B.p.value,anova_result$AB.p.value))%>%
        mutate(Effect=gsub(x=Effect,pattern="age_group",replacement="Age Group"),
               Effect=gsub(x=Effect,pattern="cand_RDM_names",replacement="Model"))%>%
        mutate_if(is.numeric,round,3)%>%
        mutate_at(vars(matches("p")),
                  funs(gsub(x=.,pattern="^0$",replacement="<.001")))%>%
        # mutate_at(vars(matches("p|p[GG]|p[HF]")),
        #           funs(gsub(x=.,pattern="^0(?=\\.)",replacement="",perl = T)))%>%
        write.csv(file = file.path(PATH_RSA_ANOVA_RESULT_ROBUST,
                                   paste0("Way_2_R~agexmodel_anova_table_",
                                          name_brain_region,".csv")))
    }
    #Split by brain region and age group (1-way anova)-------------------------------------
    for (brain_region_index in 1:length(LIST_NAMES_BRAIN_REGION)){
      name_brain_region=LIST_NAMES_BRAIN_REGION[brain_region_index]
      
      for (age_group_index in 1:length(LIST_NAMES_AGE_GROUP)){
        name_age_group=LIST_NAMES_AGE_GROUP[age_group_index]
        
        anova_result=list_anova_result_split_brain_region_age_group_trimmed[[brain_region_index]][[age_group_index]]
        
        data.frame(Effect="Model",
                   DFn=anova_result$df1,
                   DFd=anova_result$df2,
                   `F`=anova_result$test,
                   p=anova_result$p.value)%>%
          mutate_if(is.numeric,round,3)%>%
          mutate_at(vars(matches("p")),
                    funs(gsub(x=.,pattern="^0$",replacement="<.001")))%>%
          write.csv(file = file.path(PATH_RSA_ANOVA_RESULT_ROBUST,
                                     paste0("Way_1_R~model_anova_table_",
                                            name_brain_region,"_",
                                            name_age_group,".csv")))
      }
    }
    #Pair wise comparison for model effect--------------------------------------------------------------------------
    #Only focus on the best model (compare it to the others)-----------------------------------------------
    for (brain_region_index in 1:length(LIST_NAMES_BRAIN_REGION)){
      name_brain_region=LIST_NAMES_BRAIN_REGION[brain_region_index]
      
      for (age_group_index in 1:length(LIST_NAMES_AGE_GROUP)){
        name_age_group=LIST_NAMES_AGE_GROUP[age_group_index]
        
        anova_result=list_best_model_comparison_split_brain_region_age_group_trimmed[[name_brain_region]][[name_age_group]]
        
        anova_result$comp%>%
          write.csv(file = file.path(PATH_RSA_ANOVA_RESULT_ROBUST,
                                     paste0("Pair_to_best_model_",
                                            name_brain_region,"_",
                                            name_age_group,".csv")))
      }
    }
    
    
  }
}