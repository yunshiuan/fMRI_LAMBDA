#RSA Figure for publication [1st order MDS visualization]-------------------------------------------------
#Including neural RDMs and conceptual RDMs in the same canvas
########################################################################
#Will be Figure 2 - A (2018/5/2)

#Include the Kendall's tau test (2018/5/29)
#Remove the title, increase the space between the three dots in the legend to accomodate icons, and switch the R>L as filled
#Note:
#Including 1st order MDS visualization
#The input of MDS results is from the rsatoolbox in the MATLAB environment
#Note: I separate 1st and 2nd order MDS into different scripts, because they are essentaully distinct.
# -1st MDS is plotted in a 3D fashion (in addition to the spatial 2D),
# and the features are at condition level. (color:format,size:magnitude,shape:sign)
# -2nd MDS is plotted in a 1D fashion (in addition to the spatial 2D),
# and the features are at RDM level. (color: RDM type)

#Update the subject list for children: use "inclusive_runs_indexes_new_Jun_10.csv" (2018/6/10)
#Update the subject list for children: use "inclusive_runs_indexes_new_Jun_11.csv" (2018/6/10)

library("dplyr")
library("tidyr")
library("stringr")
library("R.matlab")
library("ggplot2")
library("psych")

#Constants--------------------------------
#Path
PATH_ROOT="D:\\Yun-Shiuan_LAMBDA"
PATH_FIGURES_OUTPUT=file.path(PATH_ROOT,"RSA","trial_22","R_figures","For_publication","Figure_MDS_1st_order_analysis_Fig3")
#File
FILE_MDS_RESULT=file.path(PATH_ROOT,"RSA","trial_22","Part2_1st_order_analysis","MDS_Results_22_10-Jun-2018.mat")
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
NAMES_ORDERED_AGE_GROUP_RAW=rev(c("all_sub","adult","x5_grade","x2_grade","child"))
NAMES_ORDERED_AGE_GROUP_PLOT=rev(c("All subjects","Adults","Fifth Graders","Second Graders","Children"))
NAMES_ORDERED_VOI_NAME_RAW=c("R_IPS","L_IPS","R_V1","L_V1","R_DLPFC")
NAMES_ORDERED_VOI_NAME_PLOT=c("Right IPS","Left IPS","Right V1","Left V1","Right DLPFC")

PAR_FORMAT_COLOR=c("#1a1aff","#ff4dff","#ff0000")
PAR_SIGN_SHAPE=c(16,1) # codes for shapes
PAR_MAG_SIZE=c(1,3,6)
#Helper function-----------------------------------------------------------------------------------
source(FILE_MDS_HELPER)

#Read in MDS results and preprocessing-------------------------------------------------------------
df_mds_coordinates=
  mds_mat_to_df_1_order(MDS_result,
                        list_stress,
                        list_age_group,
                        list_VOI)
#Convert to ordered facter to faciliate visualization
df_mds_coordinates=
  df_mds_coordinates%>%
  mutate(age_group=factor(age_group,levels=NAMES_ORDERED_AGE_GROUP_RAW),
         VOI_name=factor(VOI_name,levels=NAMES_ORDERED_VOI_NAME_RAW),
         RSA_magnitude=factor(RSA_magnitude,levels=NAMES_RSA_MAGNITUDE),
         RSA_format=factor(RSA_format,levels=NAMES_RSA_FORMAT))%>%
  mutate(age_group=factor(age_group,labels=NAMES_ORDERED_AGE_GROUP_PLOT),
         VOI_name=factor(VOI_name,labels=NAMES_ORDERED_VOI_NAME_PLOT))
#Rotate the MDS plot----------------------------------------------------------------------------------
#Get the rotation parameters
df_mds_rotation=
  df_mds_coordinates%>%
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
#Join the rotation parameters back to the df
df_mds_coordinates=
  df_mds_coordinates%>%
    left_join(df_mds_rotation,
              by = c("MDS_stress_type","age_group","VOI_name"))
#Rotate the MDS plots
df_mds_coordinates_rotated=
  df_mds_coordinates%>%
    group_by(MDS_stress_type,age_group,VOI_name,RSA_id)%>%
    do(rotated_cord_x={
        result=mds_rotate(cord_x=.$x,
                   cord_y=.$y,
                   radian=.$radian)
        rotated_cord_x=result$cord_x_rotated},
       rotated_cord_y={
         result=mds_rotate(cord_x=.$x,
                           cord_y=.$y,
                           radian=.$radian)
         rotated_cord_y=result$cord_y_rotated})%>%
  as.data.frame()%>%
  mutate_if(is.list,unlist)
#Add back the rotated coordinates
df_mds_coordinates=
  df_mds_coordinates%>%
    left_join(df_mds_coordinates_rotated,
              by = c("MDS_stress_type","age_group","VOI_name","RSA_id"))
#Reflect the plot to fix the format order----------------------------------------------
#Decide whether reflection is needed
df_mds_need_reflection=
  df_mds_coordinates%>%
    group_by(MDS_stress_type,age_group,VOI_name)%>%
    do(need_reflection={
        result=mds_format_need_reflection_or_not(cord_x=.$rotated_cord_x,
                                                 cord_y=.$rotated_cord_y,
                                                 RSA_format=.$RSA_format,
                                                 RSA_order_desired = c("FF","LL"))})%>%
    mutate_if(is.list,unlist)
df_mds_coordinates=
  df_mds_coordinates%>%
    left_join(df_mds_need_reflection,
              by = c("MDS_stress_type","age_group","VOI_name"))
#Reflection
df_mds_coordinates=
  df_mds_coordinates%>%
    group_by(MDS_stress_type,age_group,VOI_name)%>%
    do(reflected_cord_y={
         result=mds_reflect(cord_x=.$rotated_cord_x,
                            cord_y=.$rotated_cord_y,
                            mirror_vertical =unique(.$need_reflection),
                            mirror_horizontal = F)
         reflected_cord_y=result$cord_y_reflected},
       RSA_id=.$RSA_id)%>%
    unnest()%>%
    as.data.frame()%>%
    mutate_if(is.list,unlist)%>%
    right_join(df_mds_coordinates,
               by = c("MDS_stress_type","age_group","VOI_name","RSA_id"))%>%
    mutate(reflected_cord_x=rotated_cord_x)

#MDS Visualization===================================)==========================================================================================
#Detailed preprocessing and setting up for visualization (df_df_mds_visualization)
list_interested_age_groups=c("All subjects","Second Graders","Fifth Graders","Adults")
list_interested_VOI_names=list(without_DLPFC=c("Right V1","Left V1","Right IPS","Left IPS"))
list_interested_stress_type=c("stress")
list_pdf_size=list(c(8,8)) #Change across different VOI list (the more VOIs contained, the bigger the pdf)
name_plot_title_base="Conditions MDS for each neural RDM (18 conditions)"
name_file_name_base="First_order_MDS_agexVOI_"
par_scale=1.2
# #MDS Unrotated version----------------------------------------------------------
# mds_plot_1_order(df_mds_coordinates,
#                  list_interested_age_groups,
#                  list_interested_stress_type,
#                  list_interested_VOI_names,
#                  list_pdf_size,
#                  file_type_suffix=file_type,
#                  par_mag_size = c(1.2/(par_scale),3/(par_scale),6/(par_scale)),
#                  path_figures_output=PATH_FIGURES_OUTPUT,
#                  name_plot_title_base=name_plot_title_base,
#                  name_file_name_base=name_file_name_base,
#                  save_output=F,
#                  peek_first_plot=T,
#                  return_plots=FALSE)
# #MDS Rotated version----------------------------------------------------------
# df_mds_coordinates%>%
#   select(-x,-y)%>%
#   rename(x=rotated_cord_x,
#          y=rotated_cord_y)%>%
#   mds_plot_1_order(df_mds_coordinates = .,
#                    list_interested_age_groups,
#                    list_interested_stress_type,
#                    list_interested_VOI_names,
#                    list_pdf_size,
#                    file_type_suffix=".png",
#                    par_mag_size = c(1.2/(par_scale),3/(par_scale),6/(par_scale)),
#                    path_figures_output=PATH_FIGURES_OUTPUT,
#                    name_plot_title_base=name_plot_title_base,
#                    name_file_name_base=name_file_name_base,
#                    with_kendall_tau = T,
#                    save_output=T,
#                    peek_first_plot=T,
#                    return_plots=FALSE)
#Visualize: MDS Rotated and Reflected version----------------------------------------------------------
# loadfonts(device = "win")
# fonts()
# fonttable()
for (file_type in c(".pdf",".png")){
  df_mds_coordinates%>%
    select(-x,-y)%>%
    rename(x=reflected_cord_x,
           y=reflected_cord_y)%>%
    mds_plot_1_order(df_mds_coordinates = .,
                     list_interested_age_groups,
                     list_interested_stress_type,
                     list_interested_VOI_names,
                     list_pdf_size,
                     file_type_suffix=file_type,
                     par_mag_size = c(1.2/(par_scale),3/(par_scale),6/(par_scale)),
                     path_figures_output=PATH_FIGURES_OUTPUT,
                     name_plot_title_base=name_plot_title_base,
                     name_file_name_base=name_file_name_base,
                     remove_plot_title=T,
                     with_kendall_tau = T,
                     save_output=T,
                     peek_first_plot=T,
                     return_plots=FALSE)
}

#Test the Kendall's tau===============================)-----------------------------
#Constant---------------------------------------------------------------------------
#Parameter
N_SECOND_GRADE=18
N_FIFTH_GRADE=21
N_ADULT=24
N_ALL=sum(N_SECOND_GRADE,N_FIFTH_GRADE,N_ADULT)
list_interested_VOI_names=c("Right V1","Left V1","Right IPS","Left IPS")
#Path
PATH_ROOT="D:\\Yun-Shiuan_LAMBDA"
PATH_CSV_OUTPUT=file.path(PATH_ROOT,"RSA","trial_22","Tables","Table_MDS_1st_order_analysis")
dir.create(PATH_CSV_OUTPUT,recursive = T)
#Convert Kendall's tau to Pearson's R-----------------------------------------------
df_mds_coordinates=
  df_mds_coordinates%>%
    mutate(pearson_r=sin(0.5*pi*kendall_tau))
#Compare the correlation coefficients across age groups-----------------------------------------------
#n subjects of 2nd grader, 5th grader, and adult: 18/21/24
df_collect_compare_correlation_result=data.frame()
for(brain_region_index in 1:length(list_interested_VOI_names)){
  brain_region=list_interested_VOI_names[brain_region_index]
  #Extract the Pearson's R
  r_2_grade=
    df_mds_coordinates%>%
    filter(MDS_stress_type=="stress",
           VOI_name==brain_region,
           age_group=="Second Graders")%>%
    pull(pearson_r)%>%
    unique()
  r_5_grade=
    df_mds_coordinates%>%
    filter(MDS_stress_type=="stress",
           VOI_name==brain_region,
           age_group=="Fifth Graders")%>%
    pull(pearson_r)%>%
    unique()
  r_adult=
    df_mds_coordinates%>%
    filter(MDS_stress_type=="stress",
           VOI_name==brain_region,
           age_group=="Adults")%>%
    pull(pearson_r)%>%
    unique()
  
  #Test of difference between two independent correlations
  diff_r_2_5=
    r.test(n = N_SECOND_GRADE,r12 = r_2_grade,
           r34 = r_5_grade,n2 = N_FIFTH_GRADE)
  diff_r_5_adult=
    r.test(n = N_FIFTH_GRADE,r12 = r_5_grade,
           r34 = r_adult,n2 = N_ADULT)
  diff_r_2_adult=
    r.test(n = N_SECOND_GRADE,r12 = r_2_grade,
           r34 = r_adult,n2 = N_ADULT)
  
  #Collect the result
  df_compare_correlation_result=data.frame(comparison=c("diff_r_2_5","diff_r_5_adult","diff_r_2_adult"),
                                           z=c(diff_r_2_5$z,diff_r_5_adult$z,diff_r_2_adult$z),
                                           p=c(diff_r_2_5$p,diff_r_5_adult$p,diff_r_2_adult$p),
                                           brain_region=brain_region,
                                           stringsAsFactors = F)
  df_collect_compare_correlation_result=
    df_collect_compare_correlation_result%>%
    rbind.data.frame(df_compare_correlation_result)
}
df_collect_compare_correlation_result=
  df_collect_compare_correlation_result%>%
    mutate(z=round(z,3),
           p=round(p,3),
           p=ifelse(p==0,yes = "<.001",no = p))
write.csv(x = df_collect_compare_correlation_result,
          file = file.path(PATH_CSV_OUTPUT,"Compare_kendall_tau_across_age_group.csv"))

#Compare the correlation coefficients across brain region-----------------------------------------------
#n subjects of 2nd grader, 5th grader, and adult: 18/21/24
names_variable_brain_region=c("r_all_r_ips","r_all_l_ips",
                              "r_all_r_v1","r_all_l_v1")
for(brain_region_index in 1:length(list_interested_VOI_names)){
  brain_region=list_interested_VOI_names[brain_region_index]
  name_variable=names_variable_brain_region[brain_region_index]
  assign(x = name_variable,
         value = 
          df_mds_coordinates%>%
          filter(MDS_stress_type=="stress",
                 VOI_name==brain_region,
                 age_group=="All subjects")%>%
          pull(pearson_r)%>%
          unique())
}

#Test of difference between two independent correlations
diff_r_right_ips_v1=
  r.test(n = N_ALL,r12 = r_all_r_ips,
         r34 = r_all_r_v1)
diff_l_left_ips_v1=
  r.test(n = N_ALL,r12 = r_all_l_ips,
         r34 = r_all_l_v1)

#Collect the result
df_compare_correlation_result=data.frame(comparison=c("diff_r_right_ips_v1","diff_l_left_ips_v1"),
                                         z=c(diff_r_right_ips_v1$z,diff_l_left_ips_v1$z),
                                         p=c(diff_r_right_ips_v1$p,diff_l_left_ips_v1$p),
                                         stringsAsFactors = F)
df_compare_correlation_result=
  df_compare_correlation_result%>%
  mutate(z=round(z,3),
         p=round(p,3),
         p=ifelse(p==0,yes = "<.001",no = p))
write.csv(x = df_compare_correlation_result,
          file = file.path(PATH_CSV_OUTPUT,"Compare_kendall_tau_across_brain_region.csv"))
