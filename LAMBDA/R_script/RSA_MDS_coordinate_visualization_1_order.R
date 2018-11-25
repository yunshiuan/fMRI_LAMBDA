# Representational Similarity Analysis(RSA)[1st order MDS visualization]----------------------------------
#Note:
#Including 1st order MDS visualization
#The input of MDS results is from the rsatoolbox in the MATLAB environment
#Note: I separate 1st and 2nd order MDS into different scripts, because they are essentaully distinct.
# -1st MDS is plotted in a 3D fashion (in addition to the spatial 2D),
# and the features are at condition level. (color:format,size:magnitude,shape:sign)
# -2nd MDS is plotted in a 1D fashion (in addition to the spatial 2D),
# and the features are at RDM level. (color: RDM type)
library("dplyr")
library("tidyr")
library("stringr")
library("R.matlab")
library("ggplot2")

#Constants--------------------------------
#Path
PATH_ROOT="D:\\Yun-Shiuan_LAMBDA"
PATH_FIGURES_OUTPUT=file.path(PATH_ROOT,"RSA","trial_12","Part2_1st_order_analysis")
#File
FILE_MDS_RESULT=file.path(PATH_ROOT,"RSA","trial_12","Part2_1st_order_analysis","MDS_Results_12_11-Apr-2018.mat")
MDS_result=readMat(FILE_MDS_RESULT)$list.MDS.results
FILE_MDS_HELPER="D:\\GoogleDrive\\Lambda_code\\R_script\\R_functions\\RSA_MDS_visualization_functions.R"

#Parameter
list_stress=gsub(x = dimnames(MDS_result)[[1]],pattern = "\\.",replacement = "_")
list_age_group=gsub(x = dimnames(MDS_result[[1]])[[1]],pattern = "\\.",replacement = "_")
list_VOI=gsub(x = dimnames(MDS_result[[1]][[1]])[[1]],pattern = "\\.",replacement = "_")
list_results_fields=gsub(x = dimnames(MDS_result[[1]][[1]][[1]])[[1]],pattern = "\\.",replacement = "_")

INDEX_COORDINATE=1
# INDEX_DISSIMILARITIES=2
# INDEX_DISTANCES=3
# INDEX_DISPARITIES=4
# INDEX_STRESS=5

NUM_CONDITION=18
NAMES_RSA_FORMAT=c("FF","CN","LL")
NAMES_RSA_MAGNITUDE=c("N","M","F")
NAMES_RDS_SIGNS=c("+","-")#+:Left is larger
PAR_FORMAT_COLOR=c("#1a1aff","#ff4dff","#ff0000")
PAR_SIGN_SHAPE=c(1,16) # codes for shapes
PAR_MAG_SIZE=c(1,3,6)
#Helper function-------------------------------------------------------------------
source(FILE_MDS_HELPER)

#Read in MDS results and preprocessing-----------------------------------------------------------
df_mds_coordinates=
  mds_mat_to_df_1_order(MDS_result,
                        list_stress,
                        list_age_group,
                        list_VOI)
#Convert to ordered facter to faciliate visualization
df_mds_coordinates=
  df_mds_coordinates%>%
    mutate(age_group=factor(age_group,levels=c("all_sub","adult","x5_grade","x2_grade","child")),
           VOI_name=factor(VOI_name,levels=c("R_IPS","L_IPS","R_V1","L_V1","R_DLPFC")),
           RSA_magnitude=factor(RSA_magnitude,levels=c("N","M","F")),
           RSA_format=factor(RSA_format,levels=c("FF","CN","LL")))%>%
    mutate(age_group=factor(age_group,labels=c("All subjects","Adults","Fifth Graders","Second Graders","Children")),
           VOI_name=factor(VOI_name,labels=c("Right IPS","Left IPS","Right V1","Left V1","Right DLPFC")))

#MDS Visualization==========================================================================================
#Detailed preprocessing and setting up for visualization (df_df_mds_visualization)
list_interested_age_groups=c("All subjects","Second Graders","Fifth Graders","Adults")
list_interested_VOI_names=list(without_DLPFC=c("Right V1","Left V1","Right IPS","Left IPS"))
list_interested_stress_type=c("stress","metricstress")
list_pdf_size=list(c(8,8)) #Change acrossdifferent VOI list (the more VOIs contained, the bigger the pdf)
name_plot_title_base="Conditions MDS for each neural RDM (18 conditions)"
name_file_name_base="First_order_MDS_agexVOI_"
file_type=".pdf"
par_scale=1.2
mds_plot_1_order(df_mds_coordinates,
                 list_interested_age_groups,
                 list_interested_stress_type,
                 list_interested_VOI_names,
                 list_pdf_size,
                 file_type_suffix=file_type,
                 par_mag_size = c(1.2/(par_scale),3/(par_scale),6/(par_scale)),
                 path_figures_output=PATH_FIGURES_OUTPUT,
                 name_plot_title_base=name_plot_title_base,
                 name_file_name_base=name_file_name_base,
                 save_output=TRUE,
                 peek_first_plot=FALSE,
                 return_plots=FALSE)