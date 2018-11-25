# Representational Similarity Analysis(RSA)[MDS visualization - Hybrid Models]----------------------------------
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
PATH_FIGURES_OUTPUT=file.path(PATH_ROOT,"RSA","trial_14","R_figures","MDS_2nd_order_analysis","hybrid_models")
#File
NAMES_ACC_VERSION=c("without_ACC")
FILE_MDS_RESULT=file.path(PATH_ROOT,"RSA","trial_14","Part3_2nd_order_analysis","hybrid_models",
                          NAMES_ACC_VERSION,"MDS","MDS_Results_14_27-Apr-2018.mat")
FILE_HELPER_FUNCTIONS="D:\\GoogleDrive\\Lambda_code\\R_script\\R_functions\\RSA_MDS_visualization_functions.R"
source(FILE_HELPER_FUNCTIONS)
#Parameter
INDEX_DISSIMILARITIES=3
INDEX_DISTANCES=4
INDEX_DISPARITIES=5
INDEX_STRESS=6

PAR_MEASURE_COLOR=c("#ff0000","#0000ff")
PAR_TYPE_SHAPE=c(20,1)
PAR_TYPE_LINETYPE=c("solid","blank")
NAMES_AGE_GROUP_OLD=c("x2_grade","x5_grade","adult","all")
NAMES_AGE_GROUP_NEW=c("Second Graders","Fifth Graders","Adults","All")
# Quasi-constants
list_mds_result=lapply(X = setNames(FILE_MDS_RESULT,
                                    NAMES_ACC_VERSION),
                       FUN = function(file_mds_result){
                         mds_result=readMat(file_mds_result)$list.MDS.results
                         return(mds_result)})
mds_any=list_mds_result[[1]]
list_stress=gsub(x = dimnames(mds_any)[[1]],pattern = "\\.",replacement = "_")
list_age_group=gsub(x = dimnames(mds_any[[1]])[[1]],pattern = "\\.",replacement = "_")
list_result_field=gsub(x = dimnames(mds_any[[1]][[1]])[[1]],pattern = "\\.",replacement = "_")

#2nd order: Shepard plot===============================================================================================
#Read in MDS results and preprocessing-----------------------------------------------------------
#convert to a data frame (df_shepard_collect)
df_shepard_collect=
  shepard_mat_to_df_2_order(list_MDS_result =list_mds_result,
                            list_stress = list_stress,
                            list_age_group = list_age_group,
                            list_ACC_version = NAMES_ACC_VERSION)

#Convert to ordered facter to faciliate visualization
df_shepard_collect=
  df_shepard_collect%>%
  mutate(age_group=factor(age_group,levels=NAMES_AGE_GROUP_OLD))%>%
  mutate(age_group=factor(age_group,labels=NAMES_AGE_GROUP_NEW))
#2nd order: Visualization of Shepard plot of MDS ---------------------------------------------------------------------
#Detailed preprocessing and setting up for visualization (df_df_mds_visulization)
list_interested_age_groups=c("Second Graders","Fifth Graders","Adults","All")
list_interested_VOI_names=list(all_VOI=c("Right V1","Left V1","Right IPS","Left IPS","Right DLPFC"))
#Note that the version without DLPFC should based on the MDS analysis without DLPFC,
#rather than basing on the version with DLPFC but just hide the dof of DLPFC.
# list_interested_VOI_names=list(all_VOI=c("Right V1","Left V1","Right IPS","Left IPS","Right DLPFC"),
#                                without_DLPFC=c("Right V1","Left V1","Right IPS","Left IPS"))

list_interested_stress_type=c("stress","metricstress")
list_pdf_size=list(c(12,6))
par_text_font_size=c(4)
# list_interested_VOI_names=list_interested_VOI_scope
# list_pdf_size=list_pdf_size_for_VOI_scope
# par_text_font_size=text_font_size_for_VOI_scope
# name_plot_base=file_type_suffix
g=shepard_plot_2_order(df_shepard_collect = df_shepard_collect,
                       list_interested_VOI_scope = list_interested_VOI_names,
                       list_ACC_version = NAMES_ACC_VERSION,
                       list_pdf_size_for_VOI_scope = list_pdf_size,
                       text_font_size_for_VOI_scope = par_text_font_size,
                       list_interested_age_groups = list_interested_age_groups,
                       list_interested_stress_type = list_interested_stress_type,
                       path_figures_output = PATH_FIGURES_OUTPUT,
                       file_type_suffix = ".pdf",
                       peek_first_plot = F,
                       save_output = T,
                       return_plot = F)