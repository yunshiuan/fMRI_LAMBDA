#RSA Figure for publication [RDM - Show how mixing works]-------------------------------------------------
#Showing 
#Will be Figure 2 - B (2018/5/2)
#Removing the titles and the color bar (rely on external program instead)
library(R.matlab)
library(dplyr)
library(tidyr)
library(gridExtra)
library(grid)
############################################################################
# Re-define the scale: set w from 0.5 to 1.0 (Trial 19, 2018/6/7)

#Constants---------------------------------------------------------------------
#Parameters
LIST_INTERESTED_CONCEPTUAL_RDMS=list(w_2_5_8=c("g",
                                               "gagv220","gagv250","gagv280",
                                               "agv2"),
                                     w_4_6_8=c("g",
                                               "gagv240","gagv260","gagv280",
                                               "agv2"),
                                     w_2_4_6_8=c("g",
                                                 "gagv220","gagv240","gagv260","gagv280",
                                                 "agv2"),
                                     w_1to9=c("g",
                                              paste0("gagv2",1:9,"0"),
                                              "agv2"))
BOOLEAN_SUB_PLOT_TITLE=F #Decide whether to include RDM plot title
#Path
PATH_ROOT_RSA="D:\\Yun-Shiuan_LAMBDA\\RSA"
PATH_PLOT_OUTPUT=file.path(PATH_ROOT_RSA,"trial_19","R_figures","For_publication")
#File
FILE_RDMS=file.path(PATH_ROOT_RSA,
                    "trial_19\\Part1_neural_and_conceptual_RDMs\\RDMs\\RDMs_and_options_neural_and_conceptual_trial_19_08-Jun-2018.mat")
# FILE_RDMS=file.path(PATH_ROOT_RSA,
#                     "trial_13\\Part1_neural_and_conceptual_RDMs\\RDMs\\RDMs_and_options_neural_and_conceptual_trial_13_18-Apr-2018.mat")
FILE_RSA_FUNCTIONS="D:\\GoogleDrive\\Lambda_code\\R_script\\R_functions\\RSA_conceptual_models_functions.R"

#Helper function-----------------------------
source(FILE_RSA_FUNCTIONS)
generate_ggplot_list_model_mixing=function(list_interested_conceptual_RDMs,include_sub_plot_title){
  #Read in RDMs
  list_collect_conceptual_RDMs=list()
  list_collect_conceptual_names=c()
  RDM_index=1
  
  for (conceptual_type_index in 1:length(list_collect_RDMs_conceptual)){
    
    #Name the dim names of the array by corresponding VOI names to facilitate indexing in the inner loop
    conceptual_names_in_mat=unlist(list_collect_RDMs_conceptual[[conceptual_type_index]]["name",,])
    dimnames(list_collect_RDMs_conceptual[[conceptual_type_index]])[[3]]=conceptual_names_in_mat
    conceptual_type_name=names(list_collect_RDMs_conceptual)[conceptual_type_index]
    
    for(conceptual_model_index in 1:length(list_interested_conceptual_RDMs)){
      conceptual_model_name_interested=list_interested_conceptual_RDMs[conceptual_model_index]
      RDM=list_collect_RDMs_conceptual[[conceptual_type_index]]["RDM",1,conceptual_model_name_interested]
      list_collect_conceptual_RDMs[[RDM_index]]=RDM
      list_collect_conceptual_names[RDM_index]=paste0(conceptual_type_name,
                                                      "_",
                                                      conceptual_model_name_interested)
      RDM_index=RDM_index+1
    }
  }
  names(list_collect_conceptual_RDMs)=list_collect_conceptual_names
  
  #Plot conceptual RDMs and save them into a list
  list_plot_title=rep(list_interested_conceptual_RDMs,times=length(list_collect_RDMs_conceptual))
  #Rename the plot titles
  if(include_sub_plot_title){
    list_plot_title=
      list_plot_title%>%
      gsub(.,pattern = "^g$",replacement="G-format")%>%
      # gsub(.,pattern = "^n$",replacement="N-format")%>%
      # gsub(.,pattern = "null",replacement="Null")%>%
      # gsub(.,pattern = "^a$",replacement="linear-A")%>%
      # gsub(.,pattern = "^s$",replacement="linear-S")%>%
      gsub(.,pattern = "^agv2$",replacement="log-A")%>%
      # gsub(.,pattern = "^sgv2$",replacement="log-S")%>%
      gsub(.,pattern = "^gagv2",replacement="GxA-")
  }else{
    list_plot_title=rep("",times=length(list_plot_title))
  }

  list_collect_ggplot_conceptual_RDMs=
    lapply(X = 1:length(list_collect_conceptual_RDMs),
           FUN = function(conceptual_RDM_index){
             RDM_matrix=list_collect_conceptual_RDMs[conceptual_RDM_index]
             plot_title=list_plot_title[conceptual_RDM_index]
             
             RDM=RDM_plot(RDM_matrix = RDM_matrix,
                          plot_title = plot_title,
                          plot_rank_transform_and_scale = T,
                          show_legend = F,
                          show_axis_text = F)
             return(RDM$plot)
           })
  names(list_collect_ggplot_conceptual_RDMs)=list_collect_conceptual_names
  return(list_collect_ggplot_conceptual_RDMs)
}

#Read in RDMs------------------------------------------------------------------
RDMs=readMat(FILE_RDMS)
RDMs_conceptual=RDMs[["conceptual.RDMs"]]
list_collect_RDMs_conceptual=list(mixed=RDMs_conceptual)

#Collect plots: Generate plots of all interested conceptual RDMs---------------
list_collect_ggplot_conceptual_RDMs=list()
list_collect_grob_conceptual_mixing=list()
for(w_version_index in 1:length(LIST_INTERESTED_CONCEPTUAL_RDMS)){

  #Collect plots: Generate plots of all interested conceptual RDMs---------------
  list_collect_ggplot_conceptual_RDMs[[w_version_index]]=
    generate_ggplot_list_model_mixing(list_interested_conceptual_RDMs = LIST_INTERESTED_CONCEPTUAL_RDMS[[w_version_index]],
                                      include_sub_plot_title=BOOLEAN_SUB_PLOT_TITLE)
  
  #Arrange all plots onto the same canvas------------------------------------------
  margin = theme(plot.margin = unit(c(0.4,0.4,0.4,0.4), "cm"))

  list_collect_grob_conceptual_mixing[[w_version_index]]=
    arrangeGrob(grobs=lapply(list_collect_ggplot_conceptual_RDMs[[w_version_index]],"+", margin),
                nrow=1)
}
names(list_collect_ggplot_conceptual_RDMs)=names(LIST_INTERESTED_CONCEPTUAL_RDMS)
names(list_collect_grob_conceptual_mixing)=names(LIST_INTERESTED_CONCEPTUAL_RDMS)


#Save all the plots--------------------------------------------------------------
dir.create(file.path(PATH_PLOT_OUTPUT,"Figure_RDMs_Fig2B"),recursive = T)
scale=1.7
for (w_version_index in 1:length(LIST_INTERESTED_CONCEPTUAL_RDMS)){
  for (filetype in c(".pdf",".png")){
    
    w_version_name=names(LIST_INTERESTED_CONCEPTUAL_RDMS)[w_version_index]
    g_conceptual_mixing=list_collect_grob_conceptual_mixing[[w_version_index]]
    
    plot_name=file.path(PATH_PLOT_OUTPUT,"Figure_RDMs_Fig2B",paste0("Mixing_models_",w_version_name,filetype))
    ggsave(plot = g_conceptual_mixing,filename = plot_name,
           width = 7*scale,height = 3)
  }
}