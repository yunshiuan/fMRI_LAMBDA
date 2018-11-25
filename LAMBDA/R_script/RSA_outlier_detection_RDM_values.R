#RSA: Detect potential outlier(s) in the RSA analysis based on RDM values-------------------------------------------------
library("dplyr")
library("tidyr")
library("ggplot2")
library("R.matlab")

#Constants----------------------------------------------------------------------------------------
#Parameter
LIST_INTERESTED_NEURAL_RDMS=c("R-IPS","L-IPS","R-V1","L-V1")
#Path
PATH_RSA_ROOT="D:\\Yun-Shiuan_LAMBDA\\RSA"
PATH_OUTPUT=file.path(PATH_RSA_ROOT,"trial_17","Outlier detection")
# dir.create(PATH_OUTPUT,recursive = T)
#File
FILE_ROOT="D:\\Yun-Shiuan_LAMBDA"
FILE_NEURAL_RDM=file.path(FILE_ROOT,"RSA","trial_13",
                          "Part1_neural_and_conceptual_RDMs","RDMs",
                          "RDMs_and_options_neural_and_conceptual_trial_13_18-Apr-2018.mat")
FILE_RSA_FUNCTIONS="D:\\GoogleDrive\\Lambda_code\\R_script\\R_functions\\RSA_conceptual_models_functions.R"
#subject list
FILE_SUBJECT_LIST_AND_DEMOGRAPHIC=file.path(FILE_ROOT,
                                            "Run_inclusion_info",
                                            "inclusive_runs_indexes_with_demographic_April_5.csv")
#Helper function
source(FILE_RSA_FUNCTIONS)
#Read in RDMs---------------------------------------------------------------------------------
RDMs=readMat(FILE_NEURAL_RDM)
RDMs_neural_2_grade=RDMs[["neural.RDMs.avg.2.grade"]]
RDMs_neural_5_grade=RDMs[["neural.RDMs.avg.5.grade"]]
RDMs_neural_adult=RDMs[["neural.RDMs.avg.adult"]]

list_collect_avg_neural_RDMs=list(grade_2=RDMs_neural_2_grade,
                                  grade_5=RDMs_neural_5_grade,
                                  adult=RDMs_neural_adult)
list_collect_individual_neural_RDMs=list(RDMs[["neural.RDMs"]])
#Average neural RDMs
list_collect_avg_neural_RDMs=
  RSA_read_avg_nerual_RDMS(list_collect_RDMs_neural = list_collect_avg_neural_RDMs,
                       list_interested_neural_RDMs = LIST_INTERESTED_NEURAL_RDMS)
#Individual neural RDMs
list_collect_individual_neural_RDMs=
  RSA_read_ind_nerual_RDMS(list_collect_RDMs_neural = list_collect_individual_neural_RDMs,
                       list_interested_neural_RDMs = LIST_INTERESTED_NEURAL_RDMS)
#Read in subject list and the demographic information----------------------------------------------------------------
df_subject_list=read.csv(FILE_SUBJECT_LIST_AND_DEMOGRAPHIC,
                         header = T,
                         stringsAsFactors = F)
df_subject_list=
  df_subject_list%>%
    mutate(sub_id_abrv=gsub(x = sub_id,pattern="df10|XFC",replacement=""))%>%
  select(-X,-run_num)%>%
  distinct()
#Outliers detection==========================================================)====================
#Detect outliers based on RDM values=========================================)====================
#Ecludiean distance from a RDM to the average RDM-------------------------------------------------
#Compare each RDM with its age group x brain region average RDM
list_age_group=c("grade_2","grade_5","adult")
list_brain_region=c("R-IPS","L-IPS","R-V1","L-V1")
df_collect_distance=data.frame()
for(brain_region_index in 1:length(list_brain_region)){
  brain_region=list_brain_region[brain_region_index]
  
  for(age_group_index in 1:length(list_age_group)){
    age_group=list_age_group[age_group_index]
    age_group_brain_region=paste0(age_group,"_",brain_region)
    
    # The average RDM per age group x brain region
    RDM_age_group_brain_region_avg=
      list_collect_avg_neural_RDMs[[age_group_brain_region]][[1]]
    RDM_age_group_brain_region_avg=
      vectorize_symmetric_matrix(matrix = RDM_age_group_brain_region_avg)
    
    # The subject list of this age group (for the inner subject loop to use)
    age_group_subject_list_abrv=
      df_subject_list%>%
      filter(grade_character==age_group)%>%
      pull(sub_id_abrv)
    age_group_subject_list_full=
      df_subject_list%>%
      filter(grade_character==age_group)%>%
      pull(sub_id)
    # cat(age_group_subject_list_full)
    
    for (sub_id_index in 1:length(age_group_subject_list_abrv)){
      subject_id_abrv=age_group_subject_list_abrv[sub_id_index]
      subject_id_full=age_group_subject_list_full[sub_id_index]
      
      brain_region_subject_id=paste0(brain_region,"_",subject_id_abrv)
      
      # The individual RDM of this suject per age group x brain region    
      RDM_age_group_brain_region_ind=
        list_collect_individual_neural_RDMs[[brain_region_subject_id]][[1]]
      RDM_age_group_brain_region_ind=
        vectorize_symmetric_matrix(matrix = RDM_age_group_brain_region_ind)
      
      # Distance between the subject to the average
      distance_sub_to_average=
        matrix(c(RDM_age_group_brain_region_avg,
                 RDM_age_group_brain_region_ind),
               byrow = T,
               nrow = 2)
      distance_sub_to_average=
        as.numeric(dist(distance_sub_to_average))
      # Collect the result
      df_distance=data.frame(brain_region=brain_region,
                             age_group=age_group,
                             distance=distance_sub_to_average,
                             sub_id=subject_id_full)
      df_collect_distance=
        df_collect_distance%>%
        rbind.data.frame(df_distance)
    }
  }
}
#Visualization------------------------------------------------------------------------------------
df_id=
  df_collect_distance%>%
  group_by(brain_region,age_group)%>%
  summarise(n_end=n())
sub_id_by_group=
  lapply(1:nrow(df_id),
         FUN = function(group_index){
           seq=1:df_id$n_end[group_index]
         })%>%
  unlist()
df_collect_distance$id_numeric=
  sub_id_by_group
g=df_collect_distance%>%
  mutate(sub_id=gsub(x=sub_id,pattern="df10|XFC",replacement=""))%>%
  group_by(brain_region,age_group)%>%
  ggplot(aes(x=id_numeric,y=distance))+
  geom_point(size=3)+
  ggrepel::geom_text_repel(aes(label=sub_id))+
  labs(title="The distance between each subject's RDM and the average RDM of the age group x brain region")+
  facet_grid(brain_region~age_group,switch = "y")+
  theme_bw()+
  theme(strip.background=element_rect(fill=NA, color=NA),
        #Remove grid
        panel.grid.major =  element_blank(),
        panel.grid.minor =  element_blank(),
        #Remove tics/titles/texts
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x= element_blank(),
        plot.title = element_text(hjust=0.5,size = 14,face = "bold"))

for(file_type in c(".pdf",".png")){
  file_outlier_neural_RDM=
    file.path(PATH_OUTPUT,paste0("Neural_RDMs",file_type))
  ggsave(plot = g,
         filename = file_outlier_neural_RDM,
         width = 10,height = 8)
}

