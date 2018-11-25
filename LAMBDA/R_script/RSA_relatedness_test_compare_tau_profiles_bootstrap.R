#RSA:  By bootstapping, compare the overall tau profile (tau ~ weight) for each eage group x brain region and do the non-parametric test on the distribution
#The goal is to confirm the developmental trend
#y-axis: tau value
#x-axis: weights for the GxA models
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(grid)
library(gridExtra)


#Constants----------------------------------------------------------------------
#Parameters
NUM_TRIAL=22
NUN_LAYERS=3 #age group,brain region, bootstrap index
NUM_BOOTSTRAP=1000
LIST_ACC_VERSION=c("without_ACC")#,"with_ACC")
# LIST_BASE_MODEL_TYPE=c("format_base","mag_base","all_base")
LIST_AGE_GROUP_RAW_NAMES=rev(c("x2_grade","x5_grade","adult"))
# LIST_AGE_GROUP_RAW_NAMES=c("x2_grade","x5_grade","adult","all")

LIST_AGE_GROUP_NEW_NAMES=rev(c("Second Graders","Fifth Graders","Adults"))
# LIST_AGE_GROUP_NEW_NAMES=c("Second Graders","Fifth Graders","Adults","All")

LIST_VOI_RAW_NAMES=c("R_IPS","L_IPS","R_V1","L_V1")
LIST_VOI_NEW_NAMES=c("Right IPS","Left IPS","Right V1","Left V1")

# PATTERN_BASE_FORMAT="^format_base.*"
# PATTERN_BASE_MAG="^mag_base.*"
# PATTERN_BASE_ALL="^all_base.*"

#Path
PATH_RSA_ROOT="D:\\Yun-Shiuan_LAMBDA\\RSA"
PATH_RSA_CURRENT_TRIAL=file.path(PATH_RSA_ROOT,paste0("trial_",NUM_TRIAL))
#(input of this script)
PATH_BOOTSTRAP_RELATEDNESS_TEST_RESULTS=file.path(PATH_RSA_ROOT,"trial_22",
                                        "Part3_2nd_order_analysis","hybrid_models",
                                        "bootstrap",
                                        LIST_ACC_VERSION,"relatedness_test")
PATH_BOORSTRAP_FIGURE_OUTPUT=file.path(PATH_RSA_CURRENT_TRIAL,"R_figures",
                                       "bootstrap") 
dir.create(PATH_BOORSTRAP_FIGURE_OUTPUT,recursive = T)
#File
FILE_RELATEDNESS_TEST_HELPER_FUNCTIONS="D:\\GoogleDrive\\Lambda_code\\R_script\\R_functions\\RSA_relatedness_test_visualization_functions.R"
file_rsa_relatedness_test_results=file.path(PATH_BOOTSTRAP_RELATEDNESS_TEST_RESULTS,
                                            "Relatedness_Test_Results_22_12-Jun-2018.mat")
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
                                                    pattern = "(R_|L_).*(?=_b\\d+)")),
                    boot_strap_index=unlist(str_extract_all(string=model_name,
                                                            pattern = "b\\d+(?=\\.\\d)")))%>%
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
  mutate(age_group=factor(age_group,levels = LIST_AGE_GROUP_RAW_NAMES,labels = LIST_AGE_GROUP_NEW_NAMES),
         VOI_name=factor(VOI_name,levels = LIST_VOI_RAW_NAMES,labels=LIST_VOI_NEW_NAMES))
#Add in the weight variable-----------------------------------------------------------------------------------
df_all_df_for_plot=
  df_all_df_for_plot%>%
  mutate(weight=str_extract(string=cand_RDM_names_axis_text,
                            pattern="(?<=(GxA|GxS)-)\\d+"),
         weight=case_when(cand_RDM_names_axis_text=="log-A" ~ 100,
                          cand_RDM_names_axis_text=="G-format" ~ 0,
                          TRUE ~ as.numeric(weight)))
#Transform tau to Z score for each age group x brain region x bootstrap trial------------------------------------------------------------------------------------
df_all_df_for_plot=
  df_all_df_for_plot%>%
  mutate(tau=cand_relatedness_r_avg)%>%
  group_by(age_group,VOI_name,boot_strap_index)%>%
  mutate(z_tau=scale(tau)) #Z-transoremd for each bootstrapping iteration (since priority is defined this way)
#Derive 95% CI and mean of the tau profiles------------------------------------------------------------------------------------
df_all_df_for_plot_summarise=
  df_all_df_for_plot%>%
    filter(weight%%10==0)%>%
    arrange(age_group,VOI_name,desc(cand_relatedness_r_avg))%>%
    group_by(age_group,VOI_name,weight)%>%
    summarise(mean_tau=mean(tau),
              mean_z_tau=mean(z_tau),
              # n=n(),
              se_tau=sd(tau),
              se_z_tau=sd(z_tau))
              #This is not the correct way to derive bootstrapped CI
              # ci_tau=se_tau*1.96,
              # ci_z_tau=se_z_tau*1.96)
#Correct way to derive bootstrapped CI
df_ci_tau=
  df_all_df_for_plot%>%
    filter(weight%%10==0)%>%
    arrange(age_group,VOI_name,desc(cand_relatedness_r_avg))%>%
    group_by(age_group,VOI_name,weight)%>%
    mutate(rank_tau=rank(tau,ties.method = "first"))%>%
    filter(rank_tau%in%c(25,975))%>%
    mutate(ci_tau_type=case_when(rank_tau==50 ~ "ci_upper_tau",
                                 rank_tau==950 ~ "ci_lower_tau",
                                 TRUE~"none"))%>%
    filter(ci_tau_type=="ci_upper_tau"|ci_tau_type=="ci_lower_tau")%>%
    group_by(age_group,VOI_name,weight,ci_tau_type)%>%
    summarise(ci_tau=unique(tau))%>%
    spread(ci_tau_type,ci_tau)
df_ci_z_tau=
  df_all_df_for_plot%>%
  filter(weight%%10==0)%>%
  arrange(age_group,VOI_name,desc(cand_relatedness_r_avg))%>%
  group_by(age_group,VOI_name,weight)%>%
  mutate(rank_z_tau=rank(z_tau,ties.method = "first"))%>%
  filter(rank_z_tau%in%c(25,975))%>%
  mutate(ci_z_tau_type=case_when(rank_z_tau==25 ~ "ci_upper_z_tau",
                                 rank_z_tau==975 ~ "ci_lower_z_tau",
                                 TRUE~"none"))%>%
  filter(ci_z_tau_type=="ci_upper_z_tau"|ci_z_tau_type=="ci_lower_z_tau")%>%
  group_by(age_group,VOI_name,weight,ci_z_tau_type)%>%
  summarise(ci_z_tau=unique(z_tau))%>%
  spread(ci_z_tau_type,ci_z_tau)
df_all_df_for_plot_summarise=
  df_all_df_for_plot_summarise%>%
    left_join(df_ci_tau,by=c("age_group","VOI_name","weight"))%>%
    left_join(df_ci_z_tau,by=c("age_group","VOI_name","weight"))
#Plot the tau profile (gridExtra approach)---------------------------------------------------------------------------------------
#3 age groups on the same plot
#Unused: plot with 95% CI-------------------------------------------------------------
list_mean_dv=c("mean_tau","mean_z_tau")
list_ci_upper_dv=c("ci_upper_tau","ci_upper_z_tau")
list_ci_lower_dv=c("ci_lower_tau","ci_lower_z_tau")

list_VOI_scope=list(all=c("Right IPS","Left IPS","Right V1","Left V1"),
                    ips=c("Right IPS","Left IPS"))
# dodge= position_dodge(0.5) #"bar" for dodging
plot_index=1
gg_list=c()
for(dv_index in 1:length(list_mean_dv)){
  for(voi_scope_index in 1:length(list_VOI_scope)){
    dv_mean_name=list_mean_dv[dv_index]#mean of the tau type
    dv_ci_upper_name=list_ci_upper_dv[dv_index]#ci of the tau type
    dv_ci_lower_name=list_ci_lower_dv[dv_index]#ci of the tau type
    
    VOI_scope=list_VOI_scope[[voi_scope_index]]
    gg_list[[plot_index]]=
      df_all_df_for_plot_summarise%>%
      rename(dv_mean=!!(dv_mean_name),
             dv_ci_upper=!!(dv_ci_upper_name),
             dv_ci_lower=!!(dv_ci_lower_name))%>%#Choose the dv according to the tau type
      filter(weight%%10==0,
             age_group!="All",
             VOI_name%in% VOI_scope)%>%
      as.data.frame()%>%
      # mutate(bar=paste0(weight,age_group,VOI_name))%>% #"bar" for dodging
      ggplot(aes(x=weight,
                 y=dv_mean))+
                 #color=age_group))+
      # geom_ribbon(aes(ymin=dv_mean-dv_ci,
      #                 ymax=dv_mean+dv_ci),
      #             fill = "grey70")+
      geom_errorbar(aes(ymin=dv_ci_upper,
                        ymax=dv_ci_lower),
                        # group = bar),#"bar" for dodging
                    width=2)+
                    # position=dodge)+
      geom_point(aes(shape=age_group),size=3,position=dodge)+
      geom_line(aes(linetype=age_group),size=1)+
      scale_x_continuous(breaks = seq(from=0,to=100,by=10),
                         labels = 0:10)+
      scale_linetype_manual(name="Age Groups",
                            values=c("solid","longdash","dotted"))+
      scale_shape_manual(name="Age Groups",
                         values=c(15,17,16))+
      labs(y=dv_mean_name)+
      theme_bw()+
      facet_grid(VOI_name~.)+
      theme(legend.key.width = unit(3,"line"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
    plot_index=plot_index+1
  }
}
grid.arrange(grobs=gg_list[c(1,3)],ncol=2)
grid.arrange(grobs=gg_list[c(2,4)],ncol=2)


#plot with 95% SE-------------------------------------------------------------
list_mean_dv=c("mean_tau","mean_z_tau")
list_se_dv=c("se_tau","se_z_tau")

list_VOI_scope=list(all=c("Right IPS","Left IPS","Right V1","Left V1"),
                    ips=c("Right IPS","Left IPS"))
dodge= position_dodge(0.5) #"bar" for dodging
plot_index=1
gg_list=c()
for(dv_index in 1:length(list_mean_dv)){
  for(voi_scope_index in 1:length(list_VOI_scope)){
    dv_mean_name=list_mean_dv[dv_index]#mean of the tau type
    dv_se_name=list_se_dv[dv_index]#se of the tau type

    VOI_scope=list_VOI_scope[[voi_scope_index]]
    gg_list[[plot_index]]=
      df_all_df_for_plot_summarise%>%
      rename(dv_mean=!!(dv_mean_name),
             dv_se=!!(dv_se_name))%>%#Choose the dv according to the tau type
      filter(weight%%10==0,
             age_group!="All",
             VOI_name%in% VOI_scope)%>%
      as.data.frame()%>%
      mutate(bar=paste0(weight,age_group,VOI_name))%>% #"bar" for dodging
      ggplot(aes(x=weight,
                 y=dv_mean,
             color=age_group))+
      geom_errorbar(aes(ymin=dv_mean+dv_se,
                        ymax=dv_mean-dv_se,
                    group = bar),#"bar" for dodging
                    width=3,
                    position=dodge)+
      geom_point(aes(shape=age_group),size=3,position=dodge)+
      geom_line(size=1)+
      scale_x_continuous(breaks = seq(from=0,to=100,by=10),
                         labels = 0:10)+
      # scale_linetype_manual(name="Age Groups",
      #                       values=c("solid","longdash","dotted"))+
      scale_color_discrete(name="Age Groups")+
      scale_shape_manual(name="Age Groups",
                         values=c(15,17,16))+
      labs(y=dv_mean_name)+
      theme_bw()+
      facet_grid(VOI_name~.)+
      theme(legend.key.width = unit(3,"line"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
    plot_index=plot_index+1
  }
}
#Visualize to the console
# grid.arrange(grobs=gg_list[c(1,3)],ncol=2)
# grid.arrange(grobs=gg_list[c(2,4)],ncol=2)
for (file_type in c(".png",".pdf")){
  g=arrangeGrob(grobs=gg_list[c(1,3)],ncol=2)
  ggsave(g,filename = file.path(PATH_BOORSTRAP_FIGURE_OUTPUT,
                                paste0("Bootstrap_tau_profile_IPS",file_type)),
         width = 12,height = 8)
  g=arrangeGrob(grobs=gg_list[c(2,4)],ncol=2)
  ggsave(g,filename = file.path(PATH_BOORSTRAP_FIGURE_OUTPUT,
                                paste0("Bootstrap_tau_profile_IPS_and_V1",file_type)),
         width = 10,height = 8)
}

