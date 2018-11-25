# Representational Similarity Analysis(RSA)[2nd order MDS visualization - Base Models]----------------------------------
#Note:
#Including 2nd order MDS visualization with base model results
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
PATH_FIGURES_OUTPUT=file.path(PATH_ROOT,"RSA","trial_13","R_figures","MDS_2nd_order_analysis")
#File
NAMES_ACC_VERSION=c("with_ACC","without_ACC")
FILE_MDS_RESULT=file.path(PATH_ROOT,"RSA","trial_13","Part3_2nd_order_analysis","base_models",
                          NAMES_ACC_VERSION,"MDS","MDS_Results_13_18-Apr-2018.mat")
#Parameter
INDEX_COORDINATE=1
INDEX_DOT_NAME=2
LIST_BASE_MODELS=c("g","n","null",
                   "a","s",
                   "agv1","agv2","agv3","sgv1","sgv2")
NAMES_ACC_FOR_AGE_GROUP=c("acS","acF","acA")
NAMES_AGE_GROUP_OLD=c("x2_grade","x5_grade","adult","all")
NAMES_AGE_GROUP_NEW=c("Second Graders","Fifth Graders","Adults","All")
DF_RDM_COLOR=data.frame(RDM_name=c("Right V1", "Left V1", "Right IPS", "Left IPS", "Right DLPFC",
                                   "g","n","null",
                                   "a","s",
                                   "agv1", "agv2", "agv3","sgv1","sgv2",
                                   "ACC"),
                        RDM_color=c("#000000","#000000","#000000","#000000","#000000",
                                    "#009900","#0000b3","#ffccff",
                                    "#ffad33","#ffad33",
                                    "#995c00","#995c00","#995c00","#995c00","#995c00",
                                    "#cc0000"),
                        stringsAsFactors = F)
# INDEX_DISSIMILARITIES=3
# INDEX_DISTANCES=4
# INDEX_DISPARITIES=5
# INDEX_STRESS=6

# NUM_CONDITION=18
# NAMES_RSA_FORMAT=c("FF","CN","LL")
# NAMES_RSA_MAGNITUDE=c("N","M","F")
# NAMES_RDS_SIGNS=c("+","-")#+:Left is larger
# PAR_SIGN_SHAPE=c(1,16) # codes for shapes
# PAR_MAG_SIZE=c(1,3,6)

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
  mutate(RDM_name=
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
             

#MDS Visualization==========================================================================================
#Detailed preprocessing and setting up for visualization (df_df_mds_visualization)
list_interested_age_groups=c("Second Graders","Fifth Graders","Adults","All")
list_interested_VOI_names=list(all_VOI=c("Right V1","Left V1","Right IPS","Left IPS","Right DLPFC"))
#Note that the version without DLPFC should based on the MDS analysis without DLPFC,
#rather than basing on the version with DLPFC but just hide the dof of DLPFC.
# list_interested_VOI_names=list(all_VOI=c("Right V1","Left V1","Right IPS","Left IPS","Right DLPFC"),
#                                without_DLPFC=c("Right V1","Left V1","Right IPS","Left IPS"))

list_interested_stress_type=c("stress","metricstress")
list_pdf_size=list(c(15,6))
name_plot_title_base="Second order MDS for RDMs \nCriterion = "
name_file_name_base="Second_order_MDS_"
g=mds_plot_2_order(df_mds_coordinates = df_mds_coordinates,
                 list_interested_age_groups = list_interested_age_groups,
                 list_interested_stress_type = list_interested_stress_type,
                 list_interested_VOI_names = list_interested_VOI_names,
                 list_interested_conceptual_models = LIST_BASE_MODELS,
                 list_ACC_version=NAMES_ACC_VERSION,
                 # list_pdf_size_for_VOI_scope,
                 # text_font_size_for_VOI_scope,
                 file_type_suffix=".pdf",
                 path_figures_output = PATH_FIGURES_OUTPUT,
                 name_plot_title_base=name_plot_title_base,
                 name_file_name_base=name_file_name_base,
                 df_rdm_color=DF_RDM_COLOR,
                 save_output=FALSE,
                 peek_first_plot=TRUE,
                 return_plots=FALSE)

# Already wraped the following code in the function mds_plot_2_order()
# for(acc_version_index in 1:length(NAMES_ACC_VERSION)){
#   
#   name_acc_version=NAMES_ACC_VERSION[acc_version_index]
#   
#   df_mds_visualization_ACC=
#     df_mds_coordinates%>%
#       filter(acc_version%in%name_acc_version)
#   
#   for(list_VOI_index in 1:length(list_interested_VOI_names)){
#     
#     VOI_list_type=list_interested_VOI_names[[list_VOI_index]]
#     VOI_list_name=names(list_interested_VOI_names)[list_VOI_index]
#     
#     #The interested RDMs for this iteration BASE MODELS + (with/without ACC x with/without DLPFC)
#     list_RDM_interested=c(LIST_BASE_MODELS,VOI_list_type)
#     
#     #Include ACC if intrested
#     if (name_acc_version ==  "with_ACC"){
#       list_RDM_interested=c(list_RDM_interested,"ACC")
#       #No need to include ACC if not intrested
#     }
#   
#     pdf_width=list_pdf_size[[list_VOI_index]][1]
#     pdf_height=list_pdf_size[[list_VOI_index]][2]
#     
#       for(stress_type_index in 1:length(list_interested_stress_type)){
#         
#         stress_type=list_interested_stress_type[[stress_type_index]]
#         path_output=file.path(PATH_FIGURES_OUTPUT,name_acc_version,stress_type,"MDS_plot")
#         dir.create(path_output,recursive = T)
#         file_output=file.path(path_output,paste0(name_file_name_base,VOI_list_name,file_type))
#         #Filter in  the interested scope of RDMs
#         df_mds_visualization=
#           df_mds_visualization_ACC%>%
#           filter(MDS_stress_type %in% stress_type,
#                  age_group%in%list_interested_age_groups,
#                  RDM_name%in%list_RDM_interested)
#         #To make sure both axes share the same limits in the plot
#         limit_x_range=range(df_mds_visualization$x)
#         limit_y_range=range(df_mds_visualization$y)
#         limit_common_range=c(min(limit_x_range[1],limit_y_range[1]),
#                              max(limit_x_range[2],limit_y_range[2]))
#         #plot title
#         plot_title=paste0(name_plot_title_base,stress_type)
#         #RDM colors
#         plot_RDM_color=
#           DF_RDM_COLOR%>%
#             filter(RDM_name%in%list_RDM_interested)%>%
#             pull(RDM_color)
#         
#         #plot
#         df_mds_visualization%>%
#           ggplot(aes(x=x,y=y))+
#           geom_point(aes(color=RDM_name),size=6)+ #Plot the format and mag.
#           facet_wrap(~age_group,ncol = length(list_interested_age_groups))+
#           coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE)+ #Fix the x/y axis as having same units
#           # geom_point(aes(shape=RSA_sign),color="black")+ #Overlay the sign
#           scale_x_continuous(limits = limit_common_range)+
#           scale_y_continuous(limits = limit_common_range)+
#           scale_color_manual(name="RDM Name",
#                              values =  plot_RDM_color)+
#           labs(title=plot_title)+
#           theme_bw()+
#           theme(axis.title = element_blank(),
#                 axis.ticks = element_blank(),
#                 axis.text = element_blank(),
#                 plot.title = element_text(hjust=0.5,size = 14,face = "bold"),
#                 # strip.background=element_rect(fill=NA, color=NA),
#                 strip.text = element_text(size=12),
#                 #Remove legend
#                 legend.position = "none",
#                 #Remove grid
#                 panel.grid.major =  element_blank(),
#                 panel.grid.minor =  element_blank())+
#           # ggrepel::geom_text_repel(aes(label=c(45,43)),size=4)+ #Mark the sign
#           ggrepel::geom_text_repel(aes(label=RDM_name),size=6)+ #Mark the Rdm NAME
#           ggsave(filename = file_output,
#                  width = pdf_width,height = pdf_height)
#       }
#   }
# }