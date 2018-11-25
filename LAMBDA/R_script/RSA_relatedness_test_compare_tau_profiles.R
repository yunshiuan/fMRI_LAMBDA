#RSA:  Plot the overall tau profile (tau ~ weight) for each eage group x brain region and do the non-parametric test on the distribution
#The goal is to confirm the developmental trend
#y-axis: tau value
#x-axis: weights for the GxA models
library(dplyr)
library(ggplot2)
library(stringr)
library(grid)
library(gridExtra)
####################################################################
#Update the subject list for children: use "inclusive_runs_indexes_new_Jun_10.csv" (2018/6/10)
#Update the subject list for children: use "inclusive_runs_indexes_new_Jun_11.csv" (2018/6/10)

#Constants----------------------------------------------------------------------
#Parameters
NUM_TRIAL=22
NUN_LAYERS=2 #age group,brain region
LIST_ACC_VERSION=c("without_ACC")#,"with_ACC")
# LIST_BASE_MODEL_TYPE=c("format_base","mag_base","all_base")
LIST_AGE_GROUP_RAW_NAMES=c("x2_grade","x5_grade","adult","all")
LIST_AGE_GROUP_NEW_NAMES=c("Second Graders","Fifth Graders","Adults","All")
LIST_VOI_NAMES=c("R_IPS","L_IPS","R_V1","L_V1")
# PATTERN_BASE_FORMAT="^format_base.*"
# PATTERN_BASE_MAG="^mag_base.*"
# PATTERN_BASE_ALL="^all_base.*"

#Path
PATH_RSA_ROOT="D:\\Yun-Shiuan_LAMBDA\\RSA"
PATH_RSA_CURRENT_TRIAL=file.path(PATH_RSA_ROOT,paste0("trial_",NUM_TRIAL))
#(input of this script)
PATH_RELATEDNESS_TEST_RESULTS=file.path(PATH_RSA_ROOT,"trial_22",
                                        "Part3_2nd_order_analysis","hybrid_models",
                                        LIST_ACC_VERSION,"relatedness_test")
#File
FILE_RELATEDNESS_TEST_HELPER_FUNCTIONS="D:\\GoogleDrive\\Lambda_code\\R_script\\R_functions\\RSA_relatedness_test_visualization_functions.R"
file_rsa_relatedness_test_results=file.path(PATH_RELATEDNESS_TEST_RESULTS,
                                            "Relatedness_Test_Results_22_10-Jun-2018.mat")
#Helper functions---------------------------------------------------------------
source(FILE_RELATEDNESS_TEST_HELPER_FUNCTIONS)

#Plot all relatedess test results-----------------------------------------------
acc_version_index=1
file_rsa_result_mat=file_rsa_relatedness_test_results[acc_version_index]

#Read in the mat file and flatten it into a list-------------------------------
list_rsa_result=RSA_relatedness_nested_mat_to_single_mats(file_rsa_result_mat,NUN_LAYERS)

#Convert the flatten list of list to a list of data frame----------------------
list_df_rsa_result=lapply(list_rsa_result,FUN = RSA_relatedness_single_mat_to_df)

# #Merge the "df_relatedness_tidy_plot" with the same type of base model (one dataframe for each base model type)----------
list_merged_df_relatedness_tidy_plot=
  lapply(1,#:length(LIST_BASE_MODEL_TYPE),
         FUN = function(base_model_type_index){
           #Indexing the sub-list of the base model type-----------------------
           list_df_rbind=do.call("rbind",list_df_rsa_result)
           #Extract the dfs for plotting and merge them-------------------------
           list_df_tidy_plot=list_df_rbind[,"df_relatedness_tidy_plot"]
           df_tidy_plot_merged=do.call("rbind.data.frame",list_df_tidy_plot)
           df_tidy_plot_merged=
             df_tidy_plot_merged%>%
             tibble::rownames_to_column(var = "model_name")%>%
             mutate(age_group=unlist(str_extract_all(string=model_name,
                                                     pattern = "^.*(?=_R_|_L_)")),
                    VOI_name=unlist(str_extract_all(string=model_name,
                                                    pattern = "(R_|L_).*(?=\\.\\d)")))%>%
             select(-model_name)
           return(df_tidy_plot_merged)
         })
#Further merge the merged list into a single data frame (one data frame for all base model types)------------
# df_all_df_for_plot=do.call("rbind.data.frame",list_merged_df_relatedness_tidy_plot)
df_all_df_for_plot = list_merged_df_relatedness_tidy_plot[[1]]
df_all_df_for_plot =
  df_all_df_for_plot%>%
  mutate(cand_RDM_names_axis_text=cand_RDM_names)%>%
  # Relabel the model names so that they look neater as axis texts
  mutate(cand_RDM_names_axis_text=gsub(x=cand_RDM_names_axis_text,pattern="^g$",replacement="G-format"),
         cand_RDM_names_axis_text=gsub(x=cand_RDM_names_axis_text,pattern="^gagv2(?=\\d+$)",replacement="GxA-",perl = T),
         cand_RDM_names_axis_text=gsub(x=cand_RDM_names_axis_text,pattern="^gsgv2(?=\\d+$)",replacement="GxS-",perl = T),
         cand_RDM_names_axis_text=gsub(x=cand_RDM_names_axis_text,pattern="^agv2$",replacement="log-A"),
         cand_RDM_names_axis_text=gsub(x=cand_RDM_names_axis_text,pattern="^sgv2$",replacement="log-S"))%>% 
  mutate(age_group=factor(age_group,levels = LIST_AGE_GROUP_RAW_NAMES),
         age_group=factor(age_group,labels = LIST_AGE_GROUP_NEW_NAMES),
         VOI_name=factor(VOI_name,levels = LIST_VOI_NAMES))
#Add in the weight variable-----------------------------------------------------------------------------------
df_all_df_for_plot=
  df_all_df_for_plot%>%
  mutate(weight=str_extract(string=cand_RDM_names_axis_text,
                            pattern="(?<=(GxA|GxS)-)\\d+"),
         weight=case_when(cand_RDM_names_axis_text=="log-A" ~ 100,
                          cand_RDM_names_axis_text=="G-format" ~ 0,
                          TRUE ~ as.numeric(weight)))
#Transform tau to Z score for each age group x brain region------------------------------------------------------------------------------------
df_all_df_for_plot=
  df_all_df_for_plot%>%
    mutate(tau=cand_relatedness_r_avg)%>%
    group_by(age_group,VOI_name)%>%
    mutate(z_tau=scale(tau))
  
#Plot the tau profile (gridExtra approach)---------------------------------------------------------------------------------------
#3 age groups on the same plot
list_dv=c("tau","z_tau")
list_VOI_scope=list(all=c("R_IPS","L_IPS","R_V1","L_V1"),
                    ips=c("R_IPS","L_IPS"))
plot_index=1
gg_list=c()
for(dv_index in 1:length(list_dv)){
  for(voi_scope_index in 1:length(list_VOI_scope)){
    dv_name=list_dv[dv_index]
    VOI_scope=list_VOI_scope[[voi_scope_index]]
    gg_list[[plot_index]]=
      df_all_df_for_plot%>%
        filter(weight%%10==0,
               age_group!="All",
               VOI_name%in% VOI_scope)%>%
        slice(grep(x = cand_RDM_names_axis_text,
                   pattern="GxA|G-format|log-A"))%>%
        as.data.frame()%>%
        rename(dv=!!(dv_name))%>%
        ggplot(aes(x=weight,
                   y=dv,
                   color=age_group))+
        geom_line(size=1.5)+
        scale_x_continuous(breaks = seq(from=0,to=100,by=10),
                           labels = 0:10)+
        labs(y=dv_name)+
        theme_bw()+
        facet_grid(VOI_name~.)
    plot_index=plot_index+1
  }
}
grid.arrange(grobs=gg_list[c(1,3)],ncol=2)
grid.arrange(grobs=gg_list[c(2,4)],ncol=2)


  