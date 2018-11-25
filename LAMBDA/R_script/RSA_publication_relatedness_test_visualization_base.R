#RSA figure for publicstion: Visualize the relatedness test result [for base models]------------------------------------------
library(stringr)
library(grid)
library(gridExtra)
#Update the subject list for children: use "inclusive_runs_indexes_new_Jun_10.csv" (2018/6/10)
#Update the subject list for children: use "inclusive_runs_indexes_new_Jun_11.csv" (2018/6/10)
#Collect all base models in the same plot.
#Plot with both the raw_tau and z_tau. (2018/6/11)
#Constants----------------------------------------------------------------------
#Parameters
NUM_TRIAL=22
NUN_LAYERS=3 #Base model type,age group,brain region
LIST_ACC_VERSION=c("without_ACC")
LIST_TAU_TYPE=c("z_tau","raw_tau")
LIST_SE_TYPE=c("scaled_cand_SEs","raw_cand_SEs")
LIST_BASE_MODEL_TYPE=c("format_base","mag_base","all_base")
LIST_AGE_GROUP_RAW_NAMES=c("x2_grade","x5_grade","adult","all")
LIST_AGE_GROUP_NEW_NAMES=c("Second Graders","Fifth Graders","Adults","All")
LIST_VOI_NAMES=c("R_IPS","L_IPS","R_V1","L_V1")
LIST_ORDERED_MAG_CAND_MODELS=c("log-A","linear-A","log-S","linear-S")
LIST_ORDERED_FORMAT_CAND_MODELS=c("G-format","N-format","Null")
LIST_ORDERED_ALL_BASE_MODELS=c("G-format","N-format","Null",
                               "log-A","linear-A","log-S","linear-S")
#Could be adjusted manually
acc_version_index = 1
tau_index=1
DF_RDM_COLOR=data.frame(RDM_name=c("R_V1", "L_V1", "R_IPS", "L_IPS",
                                   #Format models
                                   "g","n","null",
                                   #Magnitude models
                                   "a","s",
                                   "agv2","sgv2",
                                   #ACC
                                   "ACC"),
                        RDM_color=c("#000000","#000000","#000000","#000000",
                                    
                                    #Format models
                                    rep("#0033cc",3),#Dark blue
                                    
                                    #Magnitude models
                                    #a 
                                    "#ffd11a",#Yellow
                                    # "#2d862d",#Dark green
                                    #s
                                    "#cc2900",#Dark red
                                    
                                    #For agv2 
                                    "#ffd11a",#Yellow
                                    # "#2d862d",#Dark green
                                    #For sgv2 
                                    "#cc2900",#Dark red
                                    
                                    #ACC
                                    "#cc0000"),
                        stringsAsFactors = F)

PATTERN_BASE_FORMAT="^format_base.*"
PATTERN_BASE_MAG="^mag_base.*"
PATTERN_BASE_ALL="^all_base.*"

#Path
PATH_RSA_ROOT="D:\\Yun-Shiuan_LAMBDA\\RSA"
PATH_RSA_CURRENT_TRIAL=file.path(PATH_RSA_ROOT,paste0("trial_",NUM_TRIAL))
#input
PATH_RELATEDNESS_TEST_RESULTS=file.path(PATH_RSA_CURRENT_TRIAL,
                                        "Part3_2nd_order_analysis","base_models",
                                        LIST_ACC_VERSION,"relatedness_test")
#output
PATH_RELATEDNESS_BAR_PLOT=file.path(PATH_RSA_CURRENT_TRIAL,
                                    "R_figures","For_publication","Supplementary",
                                    "Figure_relatedness_test_base_models_FigS1A",
                                    "base_models",LIST_ACC_VERSION,LIST_TAU_TYPE)
#File
FILE_RELATEDNESS_TEST_HELPER_FUNCTIONS="D:\\GoogleDrive\\Lambda_code\\R_script\\R_functions\\RSA_relatedness_test_visualization_functions.R"
file_rsa_relatedness_test_results=
  unlist(lapply(PATH_RELATEDNESS_TEST_RESULTS,
                FUN = function(path){
                  file=list.files(path,
                                  pattern = "Relatedness_Test_Results",
                                  full.names = T)
                  return(file)
                }))
#Helper functions---------------------------------------------------------------
source(FILE_RELATEDNESS_TEST_HELPER_FUNCTIONS)

#Plot all relatedess test results-----------------------------------------------

# for(acc_version_index in 1:length(LIST_ACC_VERSION)){
tau_type=LIST_TAU_TYPE[tau_index]
se_type=LIST_SE_TYPE[tau_index]
file_rsa_result_mat=file_rsa_relatedness_test_results[acc_version_index]

#Read in the mat file and flatten it into a list-------------------------------
list_rsa_result=RSA_relatedness_nested_mat_to_single_mats(file_rsa_result_mat,NUN_LAYERS)

#Convert the flatten list of list to a list of data frame----------------------
list_df_rsa_result=lapply(list_rsa_result,FUN = RSA_relatedness_single_mat_to_df)

#Merge the "df_relatedness_tidy_plot" with the same type of base model (one dataframe for each base model type)----------
index_format_model=grep(x = names(list_df_rsa_result),pattern = PATTERN_BASE_FORMAT)
index_mag_model=grep(x = names(list_df_rsa_result),pattern = PATTERN_BASE_MAG)
index_all_model=grep(x = names(list_df_rsa_result),pattern = PATTERN_BASE_ALL)
list_index_model_type=list(index_format_model,index_mag_model,index_all_model)

list_merged_df_relatedness_tidy_plot=
  lapply(1:length(LIST_BASE_MODEL_TYPE),
         FUN = function(base_model_type_index){
           #Base model type name
           base_model_type_name=LIST_BASE_MODEL_TYPE[base_model_type_index]
           base_model_type_indexes=list_index_model_type[[base_model_type_index]]
           
           #Indexing the sub-list of the base model type-----------------------
           list_df_this_base_model=list_df_rsa_result[base_model_type_indexes]
           list_df_this_base_model_rbind=do.call("rbind",list_df_this_base_model)
           
           #Extract the dfs for plotting and merge them-------------------------
           list_df_this_base_model_tidy_plot=
             list_df_this_base_model_rbind[,"df_relatedness_tidy_plot"]
           df_this_base_model_tidy_plot_merged=
             do.call("rbind.data.frame",list_df_this_base_model_tidy_plot)
           df_this_base_model_tidy_plot_merged=
             df_this_base_model_tidy_plot_merged%>%
             tibble::rownames_to_column(var = "model_name")%>%
             mutate(base_model_type=unlist(str_extract_all(string=model_name,
                                                           pattern = "^.*_base")),
                    age_group=unlist(str_extract_all(string=model_name,
                                                     pattern = "(?<=base_).*(?=_R_|_L_)")),
                    VOI_name=unlist(str_extract_all(string=model_name,
                                                    pattern = "(R_|L_).*(?=\\.\\d)")))%>%
             select(-model_name)
           return(df_this_base_model_tidy_plot_merged)
         })

#Further merge the merged list into a single data frame (one data frame for all base model types)------------
df_all_df_for_plot=do.call("rbind.data.frame",list_merged_df_relatedness_tidy_plot)
df_all_df_for_plot=
  df_all_df_for_plot%>%
  filter(VOI_name %in% LIST_VOI_NAMES,
         cand_RDM_names !="R_DLPFC")%>%
  mutate(cand_RDM_names_axis_text=cand_RDM_names)%>%
  # Relabel the model names so that they look neater as axis texts
  mutate(cand_RDM_names_axis_text=gsub(x=cand_RDM_names_axis_text,pattern="^g$",replacement="G-format"),
         cand_RDM_names_axis_text=gsub(x=cand_RDM_names_axis_text,pattern="^n$",replacement="N-format"),
         cand_RDM_names_axis_text=gsub(x=cand_RDM_names_axis_text,pattern="^null$",replacement="Null"),
         cand_RDM_names_axis_text=gsub(x=cand_RDM_names_axis_text,pattern="^a$",replacement="linear-A"),
         cand_RDM_names_axis_text=gsub(x=cand_RDM_names_axis_text,pattern="^s$",replacement="linear-S"),
         cand_RDM_names_axis_text=gsub(x=cand_RDM_names_axis_text,pattern="^agv2$",replacement="log-A"),
         cand_RDM_names_axis_text=gsub(x=cand_RDM_names_axis_text,pattern="^sgv2$",replacement="log-S"))%>% 
  mutate(age_group=factor(age_group,levels = LIST_AGE_GROUP_RAW_NAMES),
         age_group=factor(age_group,labels = LIST_AGE_GROUP_NEW_NAMES),
         VOI_name=factor(VOI_name,levels = LIST_VOI_NAMES))%>%
  #Z-transform the tau
  group_by(age_group,VOI_name)%>%
  mutate(scale=attr(scale(cand_relatedness_r_avg),"scaled:scale"),#Save for scaling the SE
         z_tau=scale(cand_relatedness_r_avg)[,],
         raw_tau=cand_relatedness_r_avg,
         #scale the SE accordingly
         scaled_cand_SEs=cand_SEs/scale,
         raw_cand_SEs=cand_SEs
         )%>%
  as.data.frame()
  
#Plot the bar plot===================================================)-------------------------------
#Plot parameters--------------------------------------------------------
##(Not used)Facet approach
##-Could not enable the ceiling band to be adjusted automatically for each panel,
##-because facet_wrap() will not work for geom_polygon(), which is based on another data frame.
# pdf_size=c(15,12)#Width/Height
# bar_plot_mag_base=RSA_relatedness_bar_plot(df_all_df_for_plot%>%
#                                            filter(base_model_type=="mag_base"),
#                                           as_panel = T,axis_text_size = 10,
#                                           df_rdm_color = DF_RDM_COLOR,
#                                           bar_order = LIST_ORDERED_MAG_CAND_MODELS)
# bar_plot_mag_base+
#   facet_wrap(VOI_name~age_group,
#              ncol = 4,
#              scales = "free_x")+
#   theme(axis.text.x = element_text(size = 10))
# 
# bar_plot_mag_base+
#   facet_grid(VOI_name~age_group,
#              switch = "y")+
#   theme(axis.text.x = element_text(size = 10))

#Grid Extra Approach------------------------------------------------------
#Get the common limit_y for all panels
common_limit_y=c(min(min(df_all_df_for_plot%>%pull(tau_type)-df_all_df_for_plot%>%pull(se_type)),0),
                 max(df_all_df_for_plot%>%pull(tau_type)+df_all_df_for_plot%>%pull(se_type)))
#Plot and collect the figures into a list
collect_bar_plot=c()
plot_index=1
for(model_type_index in 1:length(LIST_BASE_MODEL_TYPE)){
  this_model_type_name=LIST_BASE_MODEL_TYPE[model_type_index]
  
  for(age_index in 1:length(LIST_AGE_GROUP_NEW_NAMES)){
    this_age_group_name=LIST_AGE_GROUP_NEW_NAMES[age_index]
    
    for(VOI_index in 1:length(LIST_VOI_NAMES)){
      this_VOI_name=LIST_VOI_NAMES[VOI_index]
      tryCatch({
        
        #Decide the order of bars
        if(this_model_type_name=="all_base"){
          # bar_order=NULL #Decide by the y value (Kendall's tau)
          bar_order=LIST_ORDERED_ALL_BASE_MODELS
        }else if(this_model_type_name=="format_base"){
          bar_order=LIST_ORDERED_FORMAT_CAND_MODELS
        }else if(this_model_type_name=="mag_base"){
          bar_order=LIST_ORDERED_MAG_CAND_MODELS
        }
        #Plot
        collect_bar_plot[[plot_index]]=
          df_all_df_for_plot%>%
          filter(base_model_type==this_model_type_name,
                 age_group==this_age_group_name,
                 VOI_name==this_VOI_name)%>%
          select(-cand_relatedness_r_avg,-cand_SEs)%>%# Use the tau type instead of the tau in the mat
          rename(cand_relatedness_r_avg=tau_type,
                 cand_SEs=se_type)%>%
          RSA_relatedness_bar_plot(.,as_panel = T,axis_text_size = 10,
                                   df_rdm_color = DF_RDM_COLOR,
                                   noise_ceiling = F,
                                   limit_y = common_limit_y,
                                   remove_panel_title = T,
                                   bar_order = bar_order)
      },warning=function(w){
        stop("converted from warning: ", conditionMessage(w),
             "\n model type:",this_model_type_name,
             "\n age group:",this_age_group_name,
             "\n VOI:", this_VOI_name)
      })
      plot_index=plot_index+1
    }
  }
}

#Name the plot list
list_plot_names=
  expand.grid(LIST_VOI_NAMES,LIST_AGE_GROUP_NEW_NAMES,LIST_BASE_MODEL_TYPE)%>%
  unite(plot_name,everything())%>%
  mutate(plot_name=gsub(x=plot_name,pattern=" ",replacement="_"))%>%
  pull(plot_name)

names(collect_bar_plot)=list_plot_names

#Put panels into canvases---------------------------------------------------------------
#Set up parameters
# top_title=textGrob("How related is each of the candidate RDMs to the reference RDM?",
#                    gp=gpar(fontsize=15,fontface="bold"))
# left_title=textGrob(paste0("RDM correlation\n","Kendall-taua, averaged across"," subjects"),
#                     rot = 90,
#                     gp=gpar(fontsize=15))


#Separate the full list into sub-list by different base model type
list_plot_all_base=collect_bar_plot[grep(x = names(collect_bar_plot),
                                         pattern = "all_base")]
list_plot_format_base=collect_bar_plot[grep(x = names(collect_bar_plot),
                                            pattern = "format_base")]
list_plot_mag_base=collect_bar_plot[grep(x = names(collect_bar_plot),
                                         pattern = "mag_base")]
# #Plot to the window--------------------------------------------------
# grid.arrange(grobs=list_plot_all_base,
#              ncol=length(LIST_AGE_GROUP_NEW_NAMES),as.table=F,
#              top=top_title,
#              left=left_title)
# grid.arrange(grobs=list_plot_format_base,
#              ncol=length(LIST_AGE_GROUP_NEW_NAMES),as.table=F,
#              top=top_title,
#              left=left_title)
# grid.arrange(grobs=list_plot_mag_base,
#              ncol=length(LIST_AGE_GROUP_NEW_NAMES),as.table=F,
#              top=top_title,
#              left=left_title)
#Save as files------------------------------------------------
pdf_size=c(14,9)#Width/Height
for (file_type in c(".pdf",".png")){
  path_output=PATH_RELATEDNESS_BAR_PLOT[tau_index]
  dir.create(path_output,recursive = T)
  file_suffix=file_type
  g=arrangeGrob(grobs=list_plot_all_base,
                ncol=length(LIST_AGE_GROUP_NEW_NAMES),as.table=F)
                # top=top_title,
                # left=left_title)
  ggsave(file=file.path(path_output,paste0("Relatedness_bar_plot_all_bases_agexVOI",file_suffix)),
         plot = g,width =  pdf_size[1],height =  pdf_size[2])
  
  g=arrangeGrob(grobs=list_plot_format_base,
                ncol=length(LIST_AGE_GROUP_NEW_NAMES),as.table=F)
                # top=top_title,
                # left=left_title)
  ggsave(file=file.path(path_output,paste0("Relatedness_bar_plot_format_bases_agexVOI.pdf",file_suffix)),
         plot = g,width =  pdf_size[1],height =  pdf_size[2])
  
  g=arrangeGrob(grobs=list_plot_mag_base,
                ncol=length(LIST_AGE_GROUP_NEW_NAMES),as.table=F)
                # top=top_title,
                # left=left_title)
  ggsave(file=file.path(path_output,paste0("Relatedness_bar_plot_mag_bases_agexVOI.pdf",file_suffix)),
         plot = g,width =  pdf_size[1],height =  pdf_size[2])
}
# }