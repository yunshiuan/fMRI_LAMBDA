# Representational Similarity Analysis(RSA)[MDS visualization]----------------------------------
#Note:
#Including Shepard plots for 1st order MDS Results
#The input of MDS results is from the rsatoolbox in the MATLAB environment

library("dplyr")
library("tidyr")
library("stringr")
library("R.matlab")
library("ggplot2")

#Constants--------------------------------
#Path
PATH_ROOT="D:\\Yun-Shiuan_LAMBDA"
PATH_FIGURES_OUTPUT=file.path(PATH_ROOT,"RSA","trial_12","R_figures","MDS_1st_order_analysis")
#File
FILE_MDS_RESULT=file.path(PATH_ROOT,"RSA","trial_12","Part2_1st_oder_analysis","MDS_Results_12_13-Apr-2018.mat")
MDS_result=readMat(FILE_MDS_RESULT)$list.MDS.results
FILE_HELPER_FUNCTIONS="D:\\GoogleDrive\\Lambda_code\\R_script\\R_functions\\RSA_MDS_visualization_functions.R"
source(FILE_HELPER_FUNCTIONS)
#Parameter
LIST_STRESS=gsub(x = dimnames(MDS_result)[[1]],pattern = "\\.",replacement = "_")
LIST_AGE_GROUP=gsub(x = dimnames(MDS_result[[1]])[[1]],pattern = "\\.",replacement = "_")
LIST_VOI=gsub(x = dimnames(MDS_result[[1]][[1]])[[1]],pattern = "\\.",replacement = "_")
LIST_RESULT_FIELDS=gsub(x = dimnames(MDS_result[[1]][[1]][[1]])[[1]],pattern = "\\.",replacement = "_")

# INDEX_COORDINATE=1
INDEX_DISSIMILARITIES=3
INDEX_DISTANCES=4
INDEX_DISPARITIES=5
INDEX_STRESS=6

NUM_CONDITION=18
PAR_MEASURE_COLOR=c("#ff0000","#0000ff")
PAR_TYPE_SHAPE=c(20,1)
PAR_TYPE_LINETYPE=c("solid","blank")
#1st order: Shepard plot===============================================================================================
#Read in MDS results and preprocessing-----------------------------------------------------------
#convert to a data frame (df_shepard_collect)
df_shepard_collect=
  shepard_mat_to_df_1_order(MDS_result = MDS_result,
                            list_stress = LIST_STRESS,
                            list_age_group = LIST_AGE_GROUP,
                            list_VOI = LIST_VOI)

#Convert to ordered facter to faciliate visualization
df_shepard_collect=
  df_shepard_collect%>%
  mutate(age_group=factor(age_group,levels=c("all_sub","x2_grade","x5_grade","adult","child")),
         VOI_name=factor(VOI_name,levels=c("R_V1","L_V1","R_IPS","L_IPS","R_DLPFC")))%>%
  mutate(age_group=factor(age_group,labels=c("All subjects","Second Graders","Fifth Graders","Adults","Children")),
         VOI_name=factor(VOI_name,labels=c("Right V1","Left V1","Right IPS","Left IPS","Right DLPFC")))
#1st order: Visualization of Shepard plot of MDS ---------------------------------------------------------------------
#Detailed preprocessing and setting up for visualization (df_df_mds_visulization)
list_interested_age_groups=c("All subjects","Second Graders","Fifth Graders","Adults")
list_interested_VOI_names=list(all_VOI=c("Right V1","Left V1","Right IPS","Left IPS","Right DLPFC"),
                               without_DLPFC=c("Right V1","Left V1","Right IPS","Left IPS"))
list_interested_stress_type=c("stress","metricstress")
list_pdf_size=list(c(15,25),c(12,17)) #Height/ Width; Change across different VOI list (the more VOIs contained, the bigger the pdf)
par_text_font_size=c(6,4)
# list_interested_VOI_names=list_interested_VOI_scope
# list_pdf_size=list_pdf_size_for_VOI_scope
# par_text_font_size=text_font_size_for_VOI_scope
# name_plot_base=file_type_suffix
g=shepard_plot_1_order(df_shepard_collect = df_shepard_collect,
                     list_interested_VOI_scope = list_interested_VOI_names,
                     list_pdf_size_for_VOI_scope = list_pdf_size,
                     text_font_size_for_VOI_scope = par_text_font_size,
                     list_interested_age_groups = list_interested_age_groups,
                     list_interested_stress_type = list_interested_stress_type,
                     path_figures_output = PATH_FIGURES_OUTPUT,
                     file_type_suffix = ".pdf",
                     peek_first_plot = T,
                     save_output = F,
                     return_plot = F
                     )