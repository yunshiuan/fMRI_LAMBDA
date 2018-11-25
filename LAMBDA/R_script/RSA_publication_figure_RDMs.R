#RSA Figure for publication [RDMs]-------------------------------------------------
#Including neural RDMs and conceptual RDMs in the same canvas
#Will be Figure 2 - A (2018/5/1)
#Update the subject list for children: use "inclusive_runs_indexes_new_Jun_10.csv" (2018/6/10)
#Update the subject list for children: use "inclusive_runs_indexes_new_Jun_11.csv" (2018/6/10)

library(R.matlab)
library(dplyr)
library(tidyr)
library(gridExtra)
library(grid)
#Constants---------------------------------------------------------------------
#Parameters
LIST_INTERESTED_NEURAL_RDMS=c("R-IPS","L-IPS","R-V1","L-V1")
LIST_INTERESTED_CONCEPTUAL_RDMS=c("g","n","null",
                                  "a","s",
                                  "agv2","sgv2")
BOOLEAN_PLOT_TOP_TITLE=F #Decide whether to include plot titles on the top
BOOLEAN_PLOT_SUB_TITLE=F #Decide whether to include plot titles of each RDM

#Path
PATH_ROOT_RSA="D:\\Yun-Shiuan_LAMBDA\\RSA"
PATH_PLOT_OUTPUT=file.path(PATH_ROOT_RSA,"trial_22","R_figures","For_publication")
#File
FILE_RDMS=file.path(PATH_ROOT_RSA,
                    "trial_22","Part1_neural_and_conceptual_RDMs","RDMs",
                    "RDMs_and_options_neural_and_conceptual_trial_22_10-Jun-2018.mat")
FILE_RSA_FUNCTIONS="D:\\GoogleDrive\\Lambda_code\\R_script\\R_functions\\RSA_conceptual_models_functions.R"

#Helper function
source(FILE_RSA_FUNCTIONS)
#Read in RDMs------------------------------------------------------------------
RDMs=readMat(FILE_RDMS)
RDMs_neural_2_grade=RDMs[["neural.RDMs.avg.2.grade"]]
RDMs_neural_5_grade=RDMs[["neural.RDMs.avg.5.grade"]]
RDMs_neural_adult=RDMs[["neural.RDMs.avg.adult"]]
list_collect_RDMs_neural=list(grade_2=RDMs_neural_2_grade,
                              grade_5=RDMs_neural_5_grade,
                              adult=RDMs_neural_adult)
RDMs_conceptual=RDMs[["conceptual.RDMs"]]
list_collect_RDMs_conceptual=list(base=RDMs_conceptual)

#Collect plots: Generate plots of all neural RDMs------------------------------
#Read in average neural RDMs
list_collect_neural_RDMs=
  RSA_read_avg_nerual_RDMS(list_collect_RDMs_neural = list_collect_RDMs_neural,
                           list_interested_neural_RDMs = LIST_INTERESTED_NEURAL_RDMS)
#Plot neural RDMs and save them into a list
list_plot_title=rep(LIST_INTERESTED_NEURAL_RDMS,times=length(list_collect_RDMs_neural))
list_collect_ggplot_neural_RDMs=
  lapply(X = 1:length(list_collect_neural_RDMs),
         FUN = function(neural_RDM_index){
           RDM_matrix=list_collect_neural_RDMs[neural_RDM_index]
           #Name of each neural RDM
           if(BOOLEAN_PLOT_SUB_TITLE){
              plot_title=list_plot_title[neural_RDM_index]
           }else{
              plot_title=""
            }
           RDM=RDM_plot(RDM_matrix = RDM_matrix,
                        plot_title = plot_title,
                        plot_rank_transform_and_scale = T,
                        show_legend = F,
                        show_axis_text = F)
           return(RDM$plot)
         })
names(list_collect_ggplot_neural_RDMs)=names(list_collect_neural_RDMs)

#Collect plots: Generate plots of all interested conceptual RDMs---------------
#Read in RDMs
list_collect_conceptual_RDMs=list()
list_collect_conceptual_names=c()
RDM_index=1

for (conceptual_type_index in 1:length(list_collect_RDMs_conceptual)){
  
  #Name the dim names of the array by corresponding VOI names to facilitate indexing in the inner loop
  conceptual_names_in_mat=unlist(list_collect_RDMs_conceptual[[conceptual_type_index]]["name",,])
  dimnames(list_collect_RDMs_conceptual[[conceptual_type_index]])[[3]]=conceptual_names_in_mat
  conceptual_type_name=names(list_collect_RDMs_conceptual)[conceptual_type_index]
  
  for(conceptual_model_index in 1:length(LIST_INTERESTED_CONCEPTUAL_RDMS)){
    conceptual_model_name_interested=LIST_INTERESTED_CONCEPTUAL_RDMS[conceptual_model_index]
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
list_plot_title=rep(LIST_INTERESTED_CONCEPTUAL_RDMS,times=length(list_collect_RDMs_conceptual))
#Specify the plot title for each RDM
if(BOOLEAN_PLOT_SUB_TITLE){
  list_plot_title=
    list_plot_title%>%
      gsub(.,pattern = "^g$",replacement="G-format")%>%
      gsub(.,pattern = "^n$",replacement="N-format")%>%
      gsub(.,pattern = "null",replacement="Null")%>%
      gsub(.,pattern = "^a$",replacement="linear-A")%>%
      gsub(.,pattern = "^s$",replacement="linear-S")%>%
      gsub(.,pattern = "^agv2$",replacement="log-A")%>%
      gsub(.,pattern = "^sgv2$",replacement="log-S")
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
#Name the title of each plot-----------------------------------------------------
if(BOOLEAN_PLOT_TOP_TITLE){
  title_all_in_one=textGrob("Conceptual RDMs and Neural RDMs",
                                gp=gpar(fontsize=15,fontface="bold"))
  title_grade_2=textGrob("Second Grader",gp=gpar(fontsize=15,fontface="bold"),vjust = 4)
  title_grade_5=textGrob("Fifth Grader",gp=gpar(fontsize=15,fontface="bold"),vjust = 4)
  title_adult=textGrob("Adult",gp=gpar(fontsize=15,fontface="bold"),vjust = 4)
  title_conceptual=title_all_in_one
}else{
  title_all_in_one=textGrob("")
  title_grade_2=textGrob("")
  title_grade_5=textGrob("")
  title_adult=textGrob("")
  title_conceptual=textGrob("")
}
#Arrange all plots onto the same canvas------------------------------------------
if(BOOLEAN_PLOT_TOP_TITLE){
  plot_title_all="G-format"
}else{
  plot_title_all=""
}
  
RDM_first=RDM_plot(RDM_matrix = list_collect_conceptual_RDMs$base_g,
                   plot_title = plot_title_all,
                   plot_rank_transform_and_scale = T,
                   show_legend = T,
                   show_axis_text = T)
layout=rbind(c(1,1,1,2),
             c(1,1,1,3),
             c(4,5,6,7),
             seq(from=8,length.out = 4),
             seq(from=12,length.out = 4),
             seq(from=16,length.out = 4))

g_all_in_one=arrangeGrob(grobs=c(list(RDM_first$plot),
                     list_collect_ggplot_conceptual_RDMs[c("base_n","base_null")],
                     list_collect_ggplot_conceptual_RDMs[c("base_a","base_s","base_agv2","base_sgv2")],
                     list_collect_ggplot_neural_RDMs),
            ncol=4,as.table=T,layout_matrix=layout,
            top=title_all_in_one)
#Visualize it
# grid.arrange(g_all_in_one)

#Add in subtitles and split the plots---------------------------
layout_conceptual=rbind(c(1,1,1,2),
                        c(1,1,1,3),
                        c(4,5,6,7))
#By age group
g_grade_2_neural=arrangeGrob(grobs = list_collect_ggplot_neural_RDMs[c("grade_2_R-IPS","grade_2_L-IPS",
                                                                       "grade_2_R-V1","grade_2_L-V1")],
                             ncol = 4,
                             top = title_grade_2)
g_grade_5_neural=arrangeGrob(grobs = list_collect_ggplot_neural_RDMs[c("grade_5_R-IPS","grade_5_L-IPS",
                                                                       "grade_5_R-V1","grade_5_L-V1")],
                             ncol = 4,
                             top = title_grade_5)
g_adult=arrangeGrob(grobs = list_collect_ggplot_neural_RDMs[c("adult_R-IPS","adult_L-IPS",
                                                              "adult_R-V1","adult_L-V1")],
                    ncol = 4,
                    top = title_adult)

#Conceptual models
g_conceptual=arrangeGrob(grobs=c(list(RDM_first$plot),
                               list_collect_ggplot_conceptual_RDMs[c("base_n","base_null")],
                               list_collect_ggplot_conceptual_RDMs[c("base_a","base_s","base_agv2","base_sgv2")]),
                               ncol=4,as.table=T,layout_matrix=layout_conceptual,
                               top= title_conceptual)
layout_overall=rbind(1,
                     1,
                     1,
                     2,
                     3,
                     4)
g_age_label=arrangeGrob(g_conceptual,
              g_adult,g_grade_5_neural,g_grade_2_neural,
              ncol=1,layout_matrix=layout_overall)
#Save all the plots--------------------------------------------------------------
path_plot_output_fig2a=file.path(PATH_PLOT_OUTPUT,"Figure_RDMs_Fig2A")
dir.create(path_plot_output_fig2a,recursive = T)
for (filetype in c(".pdf",".png")){
  #All without label
  plot_name=file.path(path_plot_output_fig2a,paste0("All_in_one",filetype))
  ggsave(plot = g_all_in_one,filename = plot_name,
         width = 6,height = 9)
  #All with label
  plot_name=file.path(path_plot_output_fig2a,paste0("All_in_one_with_age_label",filetype))
  ggsave(plot = g_age_label,filename = plot_name,
         width = 6,height = 9)
  #Conceptual
  plot_name=file.path(path_plot_output_fig2a,paste0("Conceptual_RDMs",filetype))
  ggsave(plot = g_conceptual,filename = plot_name,bg = "transparent",
         width = 6.5,height = 5.3)
  #Adult
  plot_name=file.path(path_plot_output_fig2a,paste0("Adult_neural_RDMs",filetype))
  ggsave(plot = g_adult,filename = plot_name,bg = "transparent",
         width = 6,height = 3)
  #5 grader
  plot_name=file.path(path_plot_output_fig2a,paste0("Fifth_grader_neural_RDMs",filetype))
  ggsave(plot = g_grade_5_neural,filename = plot_name,bg = "transparent",
         width = 6,height = 3)
  #2 grader
  plot_name=file.path(path_plot_output_fig2a,paste0("Second_grader_neural_RDMs",filetype))
  ggsave(plot = g_grade_2_neural,filename = plot_name,bg = "transparent",
         width = 6,height = 3)
}