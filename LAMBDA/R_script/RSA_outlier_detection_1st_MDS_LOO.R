#RSA: Detect potential outlier(s) in the RSA analysis based on 1st MDS results (with LOO)-------------------------------------------------
#Compute the tau value with LOO. Detect outliers based on the tau value changes while LOO.
library("dplyr")
library("tidyr")
library("stringr")
library("R.matlab")
library("ggplot2")
library("psych")
####################################################################
#Update the subject list for children: use "inclusive_runs_indexes_new_Jun_10.csv" (2018/6/10)
#Update the subject list for children: use "inclusive_runs_indexes_new_Jun_11.csv" (2018/6/10)

#Constants--------------------------------
#Path
PATH_ROOT="D:\\Yun-Shiuan_LAMBDA"
PATH_FIGURES_OUTPUT=file.path(PATH_ROOT,"RSA","trial_22","Part2_1st_order_analysis","Leave_one_out")
#File
FILE_MDS_LOO_RESULT=file.path(PATH_ROOT,"RSA","trial_22","Part2_1st_order_analysis",
                              "Leave_one_out","LOO_MDS_Results_22_10-Jun-2018.mat")
FILE_MDS_WHOLE_SAMPLE_RESULT=file.path(PATH_ROOT,"RSA","trial_22",
                                       "Part2_1st_order_analysis","MDS_Results_22_10-Jun-2018.mat")
MDS_result_LOO=readMat(FILE_MDS_LOO_RESULT)$list.MDS.results
MDS_result_whole_sample=readMat(FILE_MDS_WHOLE_SAMPLE_RESULT)$list.MDS.results
FILE_MDS_HELPER="D:\\GoogleDrive\\Lambda_code\\R_script\\R_functions\\RSA_MDS_visualization_functions.R"

#Parameter
list_stress=gsub(x = dimnames(MDS_result_LOO)[[1]],pattern = "\\.",replacement = "_")
list_age_group_LOO=gsub(x = dimnames(MDS_result_LOO[[1]])[[1]],pattern = "\\.",replacement = "_")
list_VOI_LOO=gsub(x = dimnames(MDS_result_LOO[[1]][[1]][[1]])[[1]],pattern = "\\.",replacement = "_")
list_age_group_whole_sample=gsub(x = dimnames(MDS_result_whole_sample[[1]])[[1]],pattern = "\\.",replacement = "_")
list_VOI_whole_sample=gsub(x = dimnames(MDS_result_whole_sample[[1]][[1]])[[1]],pattern = "\\.",replacement = "_")
list_results_fields=gsub(x = dimnames(MDS_result_LOO[[1]][[1]][[1]][[1]])[[1]],pattern = "\\.",replacement = "_")

INDEX_COORDINATE=1
# INDEX_DISSIMILARITIES=2
# INDEX_DISTANCES=3
# INDEX_DISPARITIES=4
# INDEX_STRESS=5

NUM_CONDITION=18
NAMES_RSA_FORMAT=c("FF","CN","LL")
NAMES_RSA_MAGNITUDE=c("N","M","F")
NAMES_RDS_SIGNS=c("+","-")#+:Left is larger
NAMES_ORDERED_AGE_GROUP_RAW_LOO=rev(c("adult","grade_5","grade_2"))
NAMES_ORDERED_AGE_GROUP_PLOT_LOO=rev(c("Adults","Fifth Graders","Second Graders"))

NAMES_ORDERED_AGE_GROUP_RAW_WHOLE_SAMPLE=rev(c("all_sub","adult","x5_grade","x2_grade","child"))
NAMES_ORDERED_AGE_GROUP_PLOT_WHOLE_SAMPLE=rev(c("All subjects","Adults","Fifth Graders","Second Graders","Children"))
NAMES_ORDERED_VOI_NAME_RAW=c("R_IPS","L_IPS","R_V1","L_V1","R_DLPFC")
NAMES_ORDERED_VOI_NAME_PLOT=c("Right IPS","Left IPS","Right V1","Left V1","Right DLPFC")
NAMES_INTERESTED_VOI=c("Right IPS","Left IPS","Right V1","Left V1")

PAR_FORMAT_COLOR=c("#1a1aff","#ff4dff","#ff0000")
PAR_SIGN_SHAPE=c(16,1) # codes for shapes
PAR_MAG_SIZE=c(1,3,6)
#Helper function-----------------------------------------------------------------------------------
source(FILE_MDS_HELPER)

#Read in LOO MDS results and preprocessing==============================================)---------------------------
df_mds_coordinates_LOO=
  mds_mat_to_df_1_order_LOO(MDS_result_LOO,
                           list_stress,
                           list_age_group_LOO,
                           list_VOI_LOO)
#Convert to ordered facter to faciliate visualization
df_mds_coordinates_LOO=
  df_mds_coordinates_LOO%>%
  mutate(age_group=factor(age_group,levels=NAMES_ORDERED_AGE_GROUP_RAW_LOO),
         VOI_name=factor(VOI_name,levels=NAMES_ORDERED_VOI_NAME_RAW),
         RSA_magnitude=factor(RSA_magnitude,levels=NAMES_RSA_MAGNITUDE),
         RSA_format=factor(RSA_format,levels=NAMES_RSA_FORMAT))%>%
  mutate(age_group=factor(age_group,labels=NAMES_ORDERED_AGE_GROUP_PLOT_LOO),
         VOI_name=factor(VOI_name,labels=NAMES_ORDERED_VOI_NAME_PLOT))
#Compute tau---------------------------------------------------------------------------------------
df_mds_rotation_LOO=
  df_mds_coordinates_LOO%>%
  mutate(magnitude_numeric=case_when(RSA_magnitude=="N"~1,
                                     RSA_magnitude=="M"~2,
                                     RSA_magnitude=="F"~3))%>%
  group_by(MDS_stress_type,age_group,VOI_name,sub_name_LOO)%>%
  do(radian={
    rotation_parameter=mds_find_mag_axis(cord_x=.$x,
                                         cord_y=.$y,
                                         magnitude=.$magnitude_numeric)
    radian=rotation_parameter$axis_radian},
    kendall_tau={
      rotation_parameter=mds_find_mag_axis(cord_x=.$x,
                                           cord_y=.$y,
                                           magnitude=.$magnitude_numeric)
      kendall_tau=rotation_parameter$kendall_tau}
  )%>%
  as.data.frame()%>%
  mutate_if(is.list,unlist)
#Order the subject by the correspondent tau
df_mds_rotation_LOO=
  df_mds_rotation_LOO%>%
    group_by(MDS_stress_type,age_group,VOI_name)%>%
    mutate(tau_rank=rank(desc(kendall_tau),ties.method ="random"),
           sub_name_LOO_abrv=gsub(x=sub_name_LOO,
                                  pattern="(XFC)|(df10)",
                                  replacement=""))
#Read in WHOLE SAMPLE MDS results and preprocessing=====================================)---------------------
df_mds_coordinates_whole_sample=
  mds_mat_to_df_1_order(MDS_result_whole_sample,
                        list_stress,
                        list_age_group_whole_sample,
                        list_VOI_whole_sample)
#Convert to ordered facter to faciliate visualization
df_mds_coordinates_whole_sample=
  df_mds_coordinates_whole_sample%>%
  mutate(age_group=factor(age_group,levels=NAMES_ORDERED_AGE_GROUP_RAW_WHOLE_SAMPLE),
         VOI_name=factor(VOI_name,levels=NAMES_ORDERED_VOI_NAME_RAW),
         RSA_magnitude=factor(RSA_magnitude,levels=NAMES_RSA_MAGNITUDE),
         RSA_format=factor(RSA_format,levels=NAMES_RSA_FORMAT))%>%
  mutate(age_group=factor(age_group,labels=NAMES_ORDERED_AGE_GROUP_PLOT_WHOLE_SAMPLE),
         VOI_name=factor(VOI_name,labels=NAMES_ORDERED_VOI_NAME_PLOT))
#Rotate the MDS plot----------------------------------------------------------------------------------
#Get the rotation parameters
df_mds_rotation_whole_sample=
  df_mds_coordinates_whole_sample%>%
  mutate(magnitude_numeric=case_when(RSA_magnitude=="N"~1,
                                     RSA_magnitude=="M"~2,
                                     RSA_magnitude=="F"~3))%>%
  group_by(MDS_stress_type,age_group,VOI_name)%>%
  do(radian={
    rotation_parameter=mds_find_mag_axis(cord_x=.$x,
                                         cord_y=.$y,
                                         magnitude=.$magnitude_numeric)
    radian=rotation_parameter$axis_radian},
    kendall_tau={
      rotation_parameter=mds_find_mag_axis(cord_x=.$x,
                                           cord_y=.$y,
                                           magnitude=.$magnitude_numeric)
      kendall_tau=rotation_parameter$kendall_tau}
  )%>%
  as.data.frame()%>%
  mutate_if(is.list,unlist)
#Merge the whole sample MDS result to the LOO MDS result (so that the full-sample-line could be plotted)-------------
df_mds_rotation_LOO=
  df_mds_rotation_LOO%>%
    left_join(df_mds_rotation_whole_sample%>%
                rename(kendall_tau_whole_sample=kendall_tau,
                       radian_whole_sample=radian),
              by=c("MDS_stress_type","age_group","VOI_name"))%>%
  as.data.frame()%>%
  mutate(age_group=factor(age_group,levels=NAMES_ORDERED_AGE_GROUP_PLOT_LOO))
#Visualize the LOO results (tau ~ subject being LOO, facet by VOI and age group)-------------------
df_mds_rotation_LOO%>%
  # filter(age_group=="Second Graders",VOI_name=="Right IPS")%>%
  filter(VOI_name%in%NAMES_INTERESTED_VOI)%>%
  ggplot(aes(x=tau_rank,y=kendall_tau))+
  geom_point()+
  ggrepel::geom_text_repel(aes(label=sub_name_LOO_abrv))+
  facet_grid(VOI_name~age_group)+
  geom_hline(aes(yintercept=kendall_tau_whole_sample),color="red")+
  theme_bw()

