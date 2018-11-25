#RSA Figure for publication [2nd order MDS visualization - Hybrid Models]-------------------------------------------------
#Including neural RDMs and conceptual RDMs in the same canvas
#Will be Figure 2 - A (2018/5/2)
#Remove title (2018/6/7)
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

#Constants--------------------------------
#Path
PATH_ROOT="D:\\Yun-Shiuan_LAMBDA"
PATH_FIGURES_OUTPUT=file.path(PATH_ROOT,"RSA","trial_22","R_figures",
                              "For_publication","Figure_MDS_2nd_order_analysis_Fig4")
#File
NAMES_ACC_VERSION=c("without_ACC")
FILE_MDS_RESULT=file.path(PATH_ROOT,"RSA","trial_22","Part3_2nd_order_analysis","hybrid_models",
                          NAMES_ACC_VERSION,"MDS","MDS_Results_22_10-Jun-2018.mat")
#Parameter
INDEX_COORDINATE=1
INDEX_DOT_NAME=2
list_w=c(10,40,70)
LIST_HYBRID_MODELS=c("g",
                     "agv2","sgv2",
                     paste0("gagv2",list_w),
                     paste0("gsgv2",list_w))
NAMES_ACC_FOR_AGE_GROUP=c("acS","acF","acA")
NAMES_AGE_GROUP_OLD=c("x2_grade","x5_grade","adult","all")
NAMES_AGE_GROUP_NEW=c("Second Graders","Fifth Graders","Adults","All Subjects")
DF_RDM_COLOR=data.frame(RDM_name=c("Right V1", "Left V1", "Right IPS", "Left IPS", #"Right DLPFC",
                                   "g",
                                   "agv2","sgv2",
                                   paste0("gagv2",list_w),
                                   paste0("gsgv2",list_w),
                                   "ACC"),
                        RDM_color=c("#000000","#000000","#000000","#000000",#"#000000",
                                    # For G-Format model
                                    "#0033cc",#Dark blue
                                    
                                    #For agv2 
                                    "#ffd11a", #Yellow
                                    # "#2d862d",#Dark green

                                    #For sgv2
                                    "#cc2900",#Dark red
                                    
                                    #For hybrid models
                                    #For GxA
                                    rep("#33cc33",3), #Green
                                    # rep("#009999",3), #Bluish-green
                                    
                                    #For GxS
                                    rep("#b300b3",3), #Purple
                                    # rep("#cc0066",3), #Pink
                                    
                                    #ACC
                                    "#cc0000"),
                        stringsAsFactors = F)

#Helper functions
FILE_MDS_HELPER="D:\\GoogleDrive\\Lambda_code\\R_script\\R_functions\\RSA_MDS_visualization_functions.R"
source(FILE_MDS_HELPER)
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
#Read in MDS results and preprocessing-----------------------------------------------------------
#And convert to a data frame (df_mds_coordinates)
df_mds_coordinates=mds_mat_to_df_2_order(list_MDS_result = list_mds_result,
                                         list_stress = list_stress,
                                         list_age_group = list_age_group,
                                         list_ACC_version = NAMES_ACC_VERSION)
#Further process to faciliate visualization
df_mds_coordinates=
  df_mds_coordinates%>%
  mutate(RDM_name_as_label=RDM_name,
         RDM_name_as_label=
           gsub(pattern="gagv2",
                x=RDM_name_as_label,replacement="GxA-",perl = T),
         RDM_name_as_label=
           gsub(pattern="gsgv2",
                x=RDM_name_as_label,replacement="GxS-",perl = T),
         #Get rid of the following zeros
         RDM_name_as_label=
           gsub(pattern="(?<=(GxA|GxS)-\\d)0",
                x=RDM_name_as_label,replacement="",perl = T),
         RDM_name_as_label=
           case_when(RDM_name=="g" ~ "G-format",
                     
                     RDM_name=="agv2" ~ "log-A",
                     RDM_name=="sgv2" ~ "log-S",
                     TRUE ~ as.character(RDM_name_as_label)),
         RDM_name=
           case_when(RDM_name=="R-V1" ~ "Right V1",
                     RDM_name=="L-V1" ~ "Left V1",
                     RDM_name=="R-IPS" ~ "Right IPS",
                     RDM_name=="L-IPS" ~ "Left IPS",
                     RDM_name=="R-DLPFC" ~ "Right DLPFC",
                     TRUE ~ as.character(RDM_name)),
         RDM_name=
           ifelse(RDM_name%in%NAMES_ACC_FOR_AGE_GROUP,yes = "ACC",no = RDM_name))%>%
  # ACC_var_name_age_group=
  #    factor(case_when(age_group==NAMES_AGE_GROUP_NEW[1]~NAMES_ACC_FOR_AGE_GROUP[1],
  #                     age_group==NAMES_AGE_GROUP_NEW[2]~NAMES_ACC_FOR_AGE_GROUP[2],
  #                     age_group==NAMES_AGE_GROUP_NEW[3]~NAMES_ACC_FOR_AGE_GROUP[3])))    
  #Convert to factor to faciliate visualization
  mutate_at(.vars = vars(c("age_group","RDM_name","acc_version","MDS_stress_type")),
            .funs = funs(factor(.,levels=unique(.))))%>%
  #Rename levels
  mutate(age_group=factor(age_group,levels=NAMES_AGE_GROUP_OLD))%>%
  mutate(age_group=factor(age_group,labels=NAMES_AGE_GROUP_NEW))%>%
  #Reorder levels (to correspond to the color order in the DF_RDM_COLOR)
  mutate(RDM_name=factor(RDM_name,levels = DF_RDM_COLOR$RDM_name))

#MDS axis rotation---------------------------------------------------------------------------
#Compute the radian
df_mds_rotation=
  df_mds_coordinates%>%
    group_by(age_group,MDS_stress_type)%>%
    do(radian={
      result=mds_find_2nd_order_axis(cord_x=.$x,cord_y=.$y,
                              RDM_name_as_label=.$RDM_name_as_label)
      })%>%
    as.data.frame()%>%
    mutate_if(is.list,unlist)
df_mds_coordinates=
  df_mds_coordinates%>%
    left_join(df_mds_rotation,
              by=c("age_group","MDS_stress_type"))
#Rotate the MDS result
df_mds_coordinates_rotated=
  df_mds_coordinates%>%
    group_by(age_group,MDS_stress_type,RDM_name_as_label)%>%
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
    unnest()%>%
    as.data.frame()
#Add back the rotated coordinates
df_mds_coordinates=
  df_mds_coordinates%>%
  left_join(df_mds_coordinates_rotated,
            by = c("MDS_stress_type","age_group","RDM_name_as_label"))
#MDS plot reflection---------------------------------------------------------------------------
#Decide whether reflection is needed
df_mds_need_reflection=
  df_mds_coordinates%>%
  group_by(age_group,MDS_stress_type)%>%
  do(need_reflection={
    result=mds_2nd_need_reflection_or_not(cord_x=.$rotated_cord_x,
                                          cord_y=.$rotated_cord_y,
                                          RDM_names=.$RDM_name_as_label,
                                          RSA_order_desired = c("log-A","log-S"))})%>%
  mutate_if(is.list,unlist)
df_mds_coordinates=
  df_mds_coordinates%>%
  left_join(df_mds_need_reflection,
            by = c("MDS_stress_type","age_group"))

#Reflection
df_mds_coordinates=
  df_mds_coordinates%>%
  group_by(MDS_stress_type,age_group)%>%
  do(reflected_cord_y={
    result=mds_reflect(cord_x=.$rotated_cord_x,
                       cord_y=.$rotated_cord_y,
                       mirror_vertical =unique(.$need_reflection),
                       mirror_horizontal = F)
    reflected_cord_y=result$cord_y_reflected},
    RDM_name_as_label=.$RDM_name_as_label)%>%
  unnest()%>%
  as.data.frame()%>%
  right_join(df_mds_coordinates,
             by = c("MDS_stress_type","age_group","RDM_name_as_label"))%>%
  mutate(reflected_cord_x=rotated_cord_x)
#MDS Visualization==========================================================================================
#Detailed preprocessing and setting up for visualization (df_df_mds_visualization)
list_interested_age_groups=c("Second Graders","Fifth Graders","Adults","All Subjects")
list_interested_VOI_names=list(all_VOI=c("Right V1","Left V1","Right IPS","Left IPS","Right DLPFC"))
#Note that the version without DLPFC should based on the MDS analysis without DLPFC,
#rather than basing on the version with DLPFC but just hiding the dot of DLPFC.
# list_interested_VOI_names=list(all_VOI=c("Right V1","Left V1","Right IPS","Left IPS","Right DLPFC"),
#                                without_DLPFC=c("Right V1","Left V1","Right IPS","Left IPS"))
par_scale=1.5
list_interested_stress_type=c("stress")
list_pdf_size=list(c(17/par_scale,6/par_scale))
name_plot_title_base="Second order MDS for both the conceptual and neural RDMs"
name_file_name_base="Second_order_MDS_"

for(file_type in c(".png",".pdf")){
  df_mds_coordinates%>%
    select(-x,-y)%>%
    rename(x=reflected_cord_x,
           y=reflected_cord_y)%>%
    # rename(x=rotated_cord_x,
    #        y=rotated_cord_y)%>%
    mds_plot_2_order(df_mds_coordinates = .,
                     list_interested_age_groups = list_interested_age_groups,
                     list_interested_stress_type = list_interested_stress_type,
                     list_interested_VOI_names = list_interested_VOI_names,
                     list_interested_conceptual_models = LIST_HYBRID_MODELS,
                     list_ACC_version = NAMES_ACC_VERSION,
                     list_pdf_size = list_pdf_size,
                     # list_pdf_size_for_VOI_scope,
                     # text_font_size_for_VOI_scope,
                     file_type_suffix=file_type,
                     path_figures_output = PATH_FIGURES_OUTPUT,
                     name_plot_title_base=name_plot_title_base,
                     name_file_name_base=name_file_name_base,
                     df_rdm_color=DF_RDM_COLOR,
                     remove_plot_title=T,
                     save_output=T,
                     peek_first_plot=T,
                     return_plots=FALSE)
}
