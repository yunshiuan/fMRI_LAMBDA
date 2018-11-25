#Visualize the relatedness test result [for hybrid models]------------------------------------------
library(stringr)
library(grid)
library(gridExtra)
#Constants----------------------------------------------------------------------
#Parameters
NUM_TRIAL=16
NUN_LAYERS=2 #age group,brain region
LIST_ACC_VERSION=c("without_ACC")#,"with_ACC")
# LIST_BASE_MODEL_TYPE=c("format_base","mag_base","all_base")
LIST_AGE_GROUP_RAW_NAMES=c("x2_grade","x5_grade","adult","all")
LIST_AGE_GROUP_NEW_NAMES=c("Second Graders","Fifth Graders","Adults","All")
LIST_VOI_NAMES=c("R_V1","L_V1","R_IPS","L_IPS")
# LIST_VOI_NAMES=c("R_V1","L_V1","R_IPS","L_IPS","R_DLPFC")
DF_RDM_COLOR=data.frame(RDM_name=c("R_V1", "L_V1", "R_IPS", "L_IPS", "R_DLPFC",
                                   "g",
                                   "agv2","sgv2",
                                   paste0("gagv2",1:9),
                                   paste0("gsgv2",1:9),
                                   "ACC"),
                        RDM_color=c("#000000","#000000","#000000","#000000","#000000",
                                    # For G-Format model
                                    "#0033cc",#Version3
                                    # "#009900",#Version 1 & 2
                                    
                                    #For agv2 and sgv2
                                    "#2d862d","#cc2900",#Version3
                                    # "#002db3","#e68a00",#Version2
                                    #"#995c00","#995c00", #Version1
                                    
                                    #For hybrid models
                                    rep("#009999",9), #Version3
                                    rep("#cc0066",9), #Version3
                                    # rep("#00ffcc",9), #Version2
                                    # rep("#ccff33",9), #Version2
                                    # rep("#cccc00",9), #Version1
                                    # rep("#cccc00",9), #Version1
                                    "#cc0000"),
                        stringsAsFactors = F)
# PATTERN_BASE_FORMAT="^format_base.*"
# PATTERN_BASE_MAG="^mag_base.*"
# PATTERN_BASE_ALL="^all_base.*"

#Path
PATH_RSA_ROOT="D:\\Yun-Shiuan_LAMBDA\\RSA"
PATH_RSA_CURRENT_TRIAL=file.path(PATH_RSA_ROOT,paste0("trial_",NUM_TRIAL))
PATH_RELATEDNESS_TEST_RESULTS=file.path(PATH_RSA_ROOT,"trial_16",
                                        "Part3_2nd_order_analysis","hybrid_models",
                                        LIST_ACC_VERSION,"relatedness_test")
PATH_RELATEDNESS_BAR_PLOT=file.path(PATH_RSA_CURRENT_TRIAL,
                                    "R_figures","relatedness_test","hybrid_models",
                                    LIST_ACC_VERSION)
#File
FILE_RELATEDNESS_TEST_HELPER_FUNCTIONS="D:\\GoogleDrive\\Lambda_code\\R_script\\R_functions\\RSA_relatedness_test_visualization_functions.R"
file_rsa_relatedness_test_results=file.path(PATH_RELATEDNESS_TEST_RESULTS,
                                            "Relatedness_Test_Results_16_28-May-2018.mat")
# file_rsa_relatedness_test_results=
#   unlist(lapply(PATH_RELATEDNESS_TEST_RESULTS,
#                 FUN = function(path){
#                   file=list.files(path,
#                                   pattern = "Relatedness_Test_Results",
#                                   full.names = T)
#                   return(file)
#                 }))

#Helper functions---------------------------------------------------------------
source(FILE_RELATEDNESS_TEST_HELPER_FUNCTIONS)

#Plot all relatedess test results-----------------------------------------------

for(acc_version_index in 1:length(LIST_ACC_VERSION)){
  
  file_rsa_result_mat=file_rsa_relatedness_test_results[acc_version_index]
  
  #Read in the mat file and flatten it into a list-------------------------------
  list_rsa_result=RSA_relatedness_nested_mat_to_single_mats(file_rsa_result_mat,NUN_LAYERS)
  
  #Convert the flatten list of list to a list of data frame----------------------
  list_df_rsa_result=lapply(list_rsa_result,FUN = RSA_relatedness_single_mat_to_df)
  
  # #Merge the "df_relatedness_tidy_plot" with the same type of base model (one dataframe for each base model type)----------
  # index_format_model=grep(x = names(list_df_rsa_result),pattern = PATTERN_BASE_FORMAT)
  # index_mag_model=grep(x = names(list_df_rsa_result),pattern = PATTERN_BASE_MAG)
  # index_all_model=grep(x = names(list_df_rsa_result),pattern = PATTERN_BASE_ALL)
  # list_index_model_type=list(index_format_model,index_mag_model,index_all_model)
  
  list_merged_df_relatedness_tidy_plot=
    lapply(1,#:length(LIST_BASE_MODEL_TYPE),
           FUN = function(base_model_type_index){
             # #Base model type name
             # base_model_type_name=LIST_BASE_MODEL_TYPE[base_model_type_index]
             # base_model_type_indexes=list_index_model_type[[base_model_type_index]]
             
             #Indexing the sub-list of the base model type-----------------------
             # list_df_this_base_model=list_df_rsa_result[base_model_type_indexes]
             # list_df_this_base_model_rbind=do.call("rbind",list_df_this_base_model)
             
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
  # #Further merge the merged list into a single data frame (one data frame for all base model types)------------
  # df_all_df_for_plot=do.call("rbind.data.frame",list_merged_df_relatedness_tidy_plot)
  df_all_df_for_plot = list_merged_df_relatedness_tidy_plot[[1]]
  df_all_df_for_plot =
    df_all_df_for_plot%>%
    mutate(cand_RDM_names_axis_text=cand_RDM_names)%>%
    # Relabel the model names so that they look neater as axis texts
    mutate(cand_RDM_names_axis_text=gsub(x=cand_RDM_names_axis_text,pattern="^g$",replacement="G-format"),
           cand_RDM_names_axis_text=gsub(x=cand_RDM_names_axis_text,pattern="^gagv2(?=\\d$)",replacement="GxA-",perl = T),
           cand_RDM_names_axis_text=gsub(x=cand_RDM_names_axis_text,pattern="^gsgv2(?=\\d$)",replacement="GxS-",perl = T),
           cand_RDM_names_axis_text=gsub(x=cand_RDM_names_axis_text,pattern="^agv2$",replacement="log-A"),
           cand_RDM_names_axis_text=gsub(x=cand_RDM_names_axis_text,pattern="^sgv2$",replacement="log-S"))%>% 
    mutate(age_group=factor(age_group,levels = LIST_AGE_GROUP_RAW_NAMES),
           age_group=factor(age_group,labels = LIST_AGE_GROUP_NEW_NAMES),
           VOI_name=factor(VOI_name,levels = LIST_VOI_NAMES))
  #Plot the bar plot----------------------------------------------------------------------------------------------
  #Plot parameters--------------------------------------------------------
  
  ##(Not used)Facet approach
  ##-Could not enable the ceiling band to be adjusted automatically for each panel,
  ##-because facet_wrap() will not work for geom_polygon(), which is based on another data frame.
  # pdf_size=c(15,12)#Width/Height
  # bar_plot_mag_base=RSA_relatedness_bar_plot(df_all_df_for_plot%>%
  #                                           filter(base_model_type=="mag_base"))
  # bar_plot_mag_base+
  #   facet_wrap(VOI_name~age_group,
  #              ncol = 3,
  #              scales = "free_x")+
  #   theme(axis.text.x = element_text(size = 10))
  
  #Grid Extra Approach------------------------------------------------------
  
  collect_bar_plot=c()
  plot_index=1

  for(age_index in 1:length(LIST_AGE_GROUP_NEW_NAMES)){
    this_age_group_name=LIST_AGE_GROUP_NEW_NAMES[age_index]
    
    for(VOI_index in 1:length(LIST_VOI_NAMES)){
      this_VOI_name=LIST_VOI_NAMES[VOI_index]
      tryCatch({
        #Plot
        collect_bar_plot[[plot_index]]=
          df_all_df_for_plot%>%
          filter(age_group==this_age_group_name,
                 VOI_name==this_VOI_name)%>%
          RSA_relatedness_bar_plot(.,as_panel = T,
                                   axis_text_size = 8,
                                   df_rdm_color = DF_RDM_COLOR)
      },warning=function(w){
        stop("converted from warning: ", conditionMessage(w),
             "\n age group:",this_age_group_name,
             "\n VOI:", this_VOI_name)
      })
      plot_index=plot_index+1
    }
  }

  #Name the plot list
  list_plot_names=
    expand.grid(LIST_VOI_NAMES,LIST_AGE_GROUP_NEW_NAMES)%>%
    unite(plot_name,everything())%>%
    mutate(plot_name=gsub(x=plot_name,pattern=" ",replacement="_"))%>%
    pull(plot_name)
  
  names(collect_bar_plot)=list_plot_names
  
  #Put panels into canvases---------------------------------------------------------------
  #Set up parameters
  top_title=textGrob("How closely is the reference RDM \nrelated to each of the candidate RDMs?",
                     gp=gpar(fontsize=15,fontface="bold"))
  left_title=textGrob(paste0("RDM correlation\n","Kendall-taua, averaged across"," subjects"),
                      rot = 90,
                      gp=gpar(fontsize=15))
  
  
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
  path_output=PATH_RELATEDNESS_BAR_PLOT[acc_version_index]
  dir.create(path_output,recursive = T)
  file_suffix=".png"
  g=arrangeGrob(grobs=collect_bar_plot,
                ncol=length(LIST_AGE_GROUP_NEW_NAMES),as.table=F,
                top=top_title,
                left=left_title)
  ggsave(file=file.path(path_output,paste0("Relatedness_bar_plot_all_bases_agexVOI",file_suffix)),
         plot = g,width =  pdf_size[1],height =  pdf_size[2])
}