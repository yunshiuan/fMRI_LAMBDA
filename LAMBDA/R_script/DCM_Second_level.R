#===========================================================================================
# LMER approach for DCM oarameters( See Dustance effect and age effect on distance effect)
# Distance Effect: t_Distance=Parameter_Far-Pamarater_Near
# Age Effect:  t_distance~Age+Gender
#===========================================================================================
library("dplyr")
library("tidyr")
library("lme4")
library("lmerTest")
library("pbapply")
#Constant----------------------------------------------------
#To be adjusted according to need
path_output="D:\\Yun-Shiuan_LAMBDA\\DCM\\Part3_second_level_LMER_no_SMA\\full_duration\\Houde_resliced_peak_4mm_C_matrix_only_V1\\all_possible_paths\\resliced_epi"
file_DCM_paths="D:\\Yun-Shiuan_LAMBDA\\DCM\\Part3_second_level_LMER_no_SMA\\full_duration\\Houde_resliced_peak_4mm_C_matrix_only_V1\\collect_all_path_coefficients.csv"
file_demographic="D:\\Yun-Shiuan_LAMBDA\\demographic\\tidy_demographic.csv"
# save.image(file = "D:\\GoogleDrive\\Lambda_code\\R_script\\RData\\DCM_Second_level.RData")
# (Part1) Load in demographic information---------------------------------
df_demo=read.csv(file = file_demographic,stringsAsFactors = F)


# (Part2) LMER: Intrinsic connections, Age effect, and Distance Effect------------------
# Read in the data frame (with all path coefficients)
df_path_collect_all_id=read.csv(file = file_DCM_paths,
                                stringsAsFactors = F)
# Some Neccesary preprocessing
# (1)Add in demographic info (only run if haven't overwrite the csv file)
df_path_collect_all_id=
  df_path_collect_all_id%>%
  left_join(df_demo%>%select(scan_age_precise,gender,sub_id),
            by = "sub_id")
# write.csv(x = df_path_collect_all_id,
#           file = "D:\\Yun-Shiuan_LAMBDA\\DCM\\Part3_second_level_LMER\\collect_all_path_coefficients.csv")
# (2)Centered age and gender (so that intercept could be interprated)
# (Not Used) Centered by subjects' mean instead of runs' mean
# (inappropriate since the mean age  of all observations is no longer zero given unequal amount of runs)
  # mean_age=
  #   df_path_collect_all_id%>%
  #   group_by(sub_id)%>%
  #   summarise(precise_age=head(scan_age_precise,1))%>%
  #   pull(precise_age)%>%mean()
  # mean_gender=
  #   df_path_collect_all_id%>%
  #   mutate(genger_dummy=ifelse(gender=="Male",yes = 1,no = -1))%>%
  #   group_by(sub_id)%>%
  #   summarise(genger_dummy=head(genger_dummy,1))%>%
  #   pull(genger_dummy)%>%mean()
  # df_path_collect_all_id=
  #   df_path_collect_all_id%>%
  #   mutate(genger_dummy=ifelse(gender=="Male",yes = 1,no = -1))%>%
  #   mutate(centered_scan_age_precise=scan_age_precise-mean_age,
  #          centered_gender_dummy=genger_dummy-mean_gender)
# (Used) Centered by runs' mean instead of subjects' mean
df_path_collect_all_id=
  df_path_collect_all_id%>%
  mutate(genger_dummy=ifelse(gender=="Male",yes = 1,no = -1))%>%
  mutate(centered_scan_age_precise=scale(scan_age_precise,scale = F)[,],
         centered_gender_dummy=scale(genger_dummy,scale = F)[,])
# Make sure to set the correct interested paths
# # (Trail1) Assume bi-hemispheres bottom-up pathways (e.g. R_V1 -> L_IPS) but no bi-hemispheres connections between same regions (e.g.,R_IPS -> L_IPS)
# A & B matrix
# interested_path=list(c("L_V1","L_IPS"),c("L_V1","R_IPS"),
#                      c("R_V1","L_IPS"),c("R_V1","R_IPS"),
#                      c("L_IPS","L_SMA"),c("L_IPS","R_SMA"),
#                      c("R_IPS","L_SMA"),c("R_IPS","R_SMA"),
#                      c("L_SMA","L_M1"),c("R_SMA","L_M1"),
#                      c("R_DLPFC","L_V1"),c("R_DLPFC","R_V1"),
#                      c("R_DLPFC","L_IPS"),c("R_DLPFC","R_IPS"),
#                      c("R_DLPFC","L_SMA"),c("R_DLPFC","R_SMA"),
#                      c("R_DLPFC","L_M1"))
# (Trail2) Assume no bi-hemispheres bottom-up pathways but bi-hemispheres connections between same regions
# interested_path=list(
#                      # uni-hemispheres bottom-up pathways
#                      c("L_V1","L_IPS"),c("R_V1","R_IPS"),
#                      c("L_IPS","L_SMA"),c("R_IPS","R_SMA"),
#                      c("L_SMA","L_M1"),
#                      # bi-hemispheres connections between same regions
#                      c("R_V1","L_V1"),c("L_V1","R_V1"),
#                      c("R_IPS","L_IPS"),c("L_IPS","R_IPS"),
#                      c("R_SMA","L_SMA"),c("L_SMA","R_SMA"),
#                      # top-down monitoring from DLPFC
#                      c("R_DLPFC","L_V1"),c("R_DLPFC","R_V1"),
#                      c("R_DLPFC","L_IPS"),c("R_DLPFC","R_IPS"),
#                      c("R_DLPFC","L_SMA"),c("R_DLPFC","R_SMA"),
#                      c("R_DLPFC","L_M1"))
# (Trial3) EDA - Test all connections
interested_path=
  df_path_collect_all_id%>%
  filter(matrix_type=="A")%>%
  select(region_from,region_to)%>%
  t()%>%
  as.data.frame(stringsAsFactors = F)%>%
  c()%>%
  unique()

#C matrix
#(Trail1) All possible driving inputs
# interested_region=c("L_V1","R_V1","L_IPS","R_IPS",
#                     "L_SMA","R_SMA","L_M1","R_DLPFC")
#(Trial2)driving inputs only on V1
interested_region=c("L_V1","R_V1")
# (Part2-A) Test the interested intrinsic connecitons (-and age-related effect) (A matrix)------------------
A_models=pblapply(interested_path,
           FUN = function(paths){
             m=lmer(formula = value~centered_gender_dummy+centered_scan_age_precise+(1|sub_id),
                    data = df_path_collect_all_id%>%
                      filter(matrix_type=="A"& region_from==paths[1]& region_to==paths[2]))
             return(m)
         })
names(A_models)=unlist(lapply(interested_path,
                              FUN = function(paths){paste0(paths[1],"-",paths[2])}))
A_models_summarized=pblapply(seq_along(A_models),
                           FUN = function(model_index){
                             model=A_models[[model_index]]
                             model_name=names(A_models)[model_index]

                             s=summary(model)
                             coeff=data.frame(s$coefficients)
                             coeff$parameter=rownames(s$coefficients)
                             names(coeff)=c("Estimate","Std_Error","df","t_value","p_value","parameter")
                             coeff$interested_path=model_name
                             return(coeff)
                           })
df_A_models_summarized=do.call("rbind",A_models_summarized)
df_A_models_summarized=
  df_A_models_summarized%>%
    mutate(significant=ifelse(p_value<0.01,yes = "**",
                              no = ifelse(p_value<0.05,yes = "*",
                                          no = ifelse(p_value<0.1,yes = "+",no = ""))))%>%
    mutate(direction=ifelse(significant!="",
                            yes = ifelse(t_value>0,yes = "pos",no = "neg"),
                            no = ""))
write.csv(x = df_A_models_summarized,
          file = paste0(path_output,"\\A_matrix_result.csv"))
# (Part2-B) Test the modulating effect - Distance effect and Age effect (B matrix)------------------
# Derive the distance effect (Near - Far)
df_B_distance_path_all_id=
  df_path_collect_all_id%>%
    select(-X)%>%
    filter(matrix_type%in%c("B_Far","B_Medium","B_Near"))%>%
    spread(key=matrix_type,value = value)%>%
    mutate(B_Distance=B_Near-B_Far)

B_models=pblapply(interested_path,
           FUN = function(paths){
             m=lmer(formula = B_Distance~centered_gender_dummy+centered_scan_age_precise+(1|sub_id),
                    data = df_B_distance_path_all_id%>%
                      filter(region_from==paths[1]& region_to==paths[2]))
             return(m)
          })
names(B_models)=unlist(lapply(interested_path,
                              FUN = function(paths){paste0(paths[1],"-",paths[2])}))
B_models_summarized=pblapply(seq_along(B_models),
                           FUN = function(model_index){
                             model=B_models[[model_index]]
                             model_name=names(B_models)[model_index]

                             s=summary(model)
                             coeff=data.frame(s$coefficients)
                             coeff$parameter=rownames(s$coefficients)
                             names(coeff)=c("Estimate","Std_Error","df","t_value","p_value","parameter")
                             coeff$interested_path=model_name
                             return(coeff)
                           })
df_B_models_summarized=do.call("rbind",B_models_summarized)
df_B_models_summarized=
  df_B_models_summarized%>%
  mutate(significant=ifelse(p_value<0.01,yes = "**",
                            no = ifelse(p_value<0.05,yes = "*",
                                        no = ifelse(p_value<0.1,yes = "+",no = ""))))%>%
  mutate(direction=ifelse(significant!="",
                          yes = ifelse(t_value>0,yes = "pos",no = "neg"),
                          no = ""))
write.csv(x = df_B_models_summarized,
          file = paste0(path_output,"\\B_Dist_effect_matrix_result.csv"))
# (Part2-C) Test the Driving effect - Distance effect and Age effect (C matrix)------------------

# Derive the distance effect (Near - Far)
df_C_distance_path_all_id=
  df_path_collect_all_id%>%
  select(-X,-region_from)%>%
  filter(matrix_type%in%c("C_Far","C_Medium","C_Near"))%>%
  spread(key=matrix_type,value = value)%>%
  mutate(C_Distance=C_Near-C_Far)

C_models=pblapply(interested_region,
                  FUN = function(region){
                    m=lmer(formula = C_Distance~centered_gender_dummy+centered_scan_age_precise+(1|sub_id),
                           data = df_C_distance_path_all_id%>%
                             filter(region_to==region))
                    return(m)
                  })
names(C_models)=interested_region

C_models_summarized=pblapply(seq_along(C_models),
                           FUN = function(model_index){
                             model=C_models[[model_index]]
                             model_name=names(C_models)[model_index]

                             s=summary(model)
                             coeff=data.frame(s$coefficients)
                             coeff$parameter=rownames(s$coefficients)
                             names(coeff)=c("Estimate","Std_Error","df","t_value","p_value","parameter")
                             coeff$interested_path=model_name
                             return(coeff)
                           })
df_C_models_summarized=do.call("rbind",C_models_summarized)
df_C_models_summarized=
  df_C_models_summarized%>%
  mutate(significant=ifelse(p_value<0.01,yes = "**",
                            no = ifelse(p_value<0.05,yes = "*",
                                        no = ifelse(p_value<0.1,yes = "+",no = ""))))%>%
  mutate(direction=ifelse(significant!="",
                          yes = ifelse(t_value>0,yes = "pos",no = "neg"),
                          no = ""))
write.csv(x = df_C_models_summarized,
          file = paste0(path_output,"\\C_Dist_effect_matrix_result.csv"))


# (Part3) Split LMER results by age group------------------------------------------------------
# (Part3-A) By age group: Test the interested intrinsic connecitons(A matrix)------------------
A_models_5_grader=pblapply(interested_path,
                  FUN = function(paths){
                    m=lmer(formula = value~centered_gender_dummy+centered_scan_age_precise+(1|sub_id),
                           data = df_path_collect_all_id%>%
                             filter(matrix_type=="A"& region_from==paths[1]& region_to==paths[2] & scan_age_precise>10)%>%
                             mutate(centered_gender_dummy=scale(centered_gender_dummy,scale = F)[,],
                                    centered_scan_age_precise=scale(centered_scan_age_precise)[,]))
                    return(m)
                  })

A_models_2_grader=pblapply(interested_path,
                           FUN = function(paths){
                             m=lmer(formula = value~centered_gender_dummy+centered_scan_age_precise+(1|sub_id),
                                    data = df_path_collect_all_id%>%
                                      filter(matrix_type=="A"& region_from==paths[1]& region_to==paths[2] & scan_age_precise<=10)%>%
                                      mutate(centered_gender_dummy=scale(centered_gender_dummy,scale = F)[,],
                                      centered_scan_age_precise=scale(centered_scan_age_precise)[,]))
                             return(m)
                           })
names(A_models_5_grader)=unlist(lapply(interested_path,
                              FUN = function(paths){paste0(paths[1],"-",paths[2])}))
names(A_models_2_grader)=unlist(lapply(interested_path,
                              FUN = function(paths){paste0(paths[1],"-",paths[2])}))

A_models_summarized_5_grader=pblapply(seq_along(A_models_5_grader),
                               FUN = function(model_index){
                                 model=A_models_5_grader[[model_index]]
                                 model_name=names(A_models_5_grader)[model_index]

                                 s=summary(model)
                                 coeff=data.frame(s$coefficients)
                                 coeff$parameter=rownames(s$coefficients)
                                 names(coeff)=c("Estimate","Std_Error","df","t_value","p_value","parameter")
                                 coeff$interested_path=model_name
                                 return(coeff)
                               })
A_models_summarized_2_grader=pblapply(seq_along(A_models_2_grader),
                                    FUN = function(model_index){
                                      model=A_models_2_grader[[model_index]]
                                      model_name=names(A_models_2_grader)[model_index]

                                      s=summary(model)
                                      coeff=data.frame(s$coefficients)
                                      coeff$parameter=rownames(s$coefficients)
                                      names(coeff)=c("Estimate","Std_Error","df","t_value","p_value","parameter")
                                      coeff$interested_path=model_name
                                      return(coeff)
                                    })
df_A_models_summarized_5_grader=do.call("rbind",A_models_summarized_5_grader)
df_A_models_summarized_2_grader=do.call("rbind",A_models_summarized_2_grader)

df_A_models_summarized_5_grader=
    df_A_models_summarized_5_grader%>%
    mutate(significant=ifelse(p_value<0.01,yes = "**",
                              no = ifelse(p_value<0.05,yes = "*",
                                          no = ifelse(p_value<0.1,yes = "+",no = ""))))%>%
    mutate(direction=ifelse(significant!="",
                            yes = ifelse(t_value>0,yes = "pos",no = "neg"),
                            no = ""))
df_A_models_summarized_2_grader=
    df_A_models_summarized_2_grader%>%
    mutate(significant=ifelse(p_value<0.01,yes = "**",
                              no = ifelse(p_value<0.05,yes = "*",
                                          no = ifelse(p_value<0.1,yes = "+",no = ""))))%>%
    mutate(direction=ifelse(significant!="",
                            yes = ifelse(t_value>0,yes = "pos",no = "neg"),
                            no = ""))

write.csv(x = df_A_models_summarized_5_grader,
          file = paste0(path_output,"\\A_matrix_result_5_grader.csv"))

write.csv(x = df_A_models_summarized_2_grader,
          file = paste0(path_output,"\\A_matrix_result_2_grader.csv"))

# (Part3-B) By age group: Test the interested distance-modulated connecitons(B matrix)------------------

B_models_5_grader=pblapply(interested_path,
                  FUN = function(paths){
                    m=lmer(formula = B_Distance~centered_gender_dummy+centered_scan_age_precise+(1|sub_id),
                           data = df_B_distance_path_all_id%>%
                             filter(region_from==paths[1]& region_to==paths[2] & scan_age_precise>10)%>%
                      mutate(centered_gender_dummy=scale(centered_gender_dummy,scale = F)[,],
                             centered_scan_age_precise=scale(centered_scan_age_precise)[,]))
                    return(m)
                  })

B_models_2_grader=pblapply(interested_path,
                           FUN = function(paths){
                             m=lmer(formula = B_Distance~centered_gender_dummy+centered_scan_age_precise+(1|sub_id),
                                    data = df_B_distance_path_all_id%>%
                                      filter(region_from==paths[1]& region_to==paths[2] & scan_age_precise<10)%>%
                               mutate(centered_gender_dummy=scale(centered_gender_dummy,scale = F)[,],
                                      centered_scan_age_precise=scale(centered_scan_age_precise)[,]))
                             return(m)
                           })
names(B_models_5_grader)=unlist(lapply(interested_path,
                              FUN = function(paths){paste0(paths[1],"-",paths[2])}))
names(B_models_2_grader)=unlist(lapply(interested_path,
                              FUN = function(paths){paste0(paths[1],"-",paths[2])}))
B_models_summarized_5_grader=pblapply(seq_along(B_models_5_grader),
                                 FUN = function(model_index){
                                   model=B_models_5_grader[[model_index]]
                                   model_name=names(B_models_5_grader)[model_index]

                                   s=summary(model)
                                   coeff=data.frame(s$coefficients)
                                   coeff$parameter=rownames(s$coefficients)
                                   names(coeff)=c("Estimate","Std_Error","df","t_value","p_value","parameter")
                                   coeff$interested_path=model_name
                                   return(coeff)
                              })
B_models_summarized_2_grader=pblapply(seq_along(B_models_2_grader),
                                    FUN = function(model_index){
                                      model=B_models_2_grader[[model_index]]
                                      model_name=names(B_models_2_grader)[model_index]

                                      s=summary(model)
                                      coeff=data.frame(s$coefficients)
                                      coeff$parameter=rownames(s$coefficients)
                                      names(coeff)=c("Estimate","Std_Error","df","t_value","p_value","parameter")
                                      coeff$interested_path=model_name
                                      return(coeff)
                                    })
df_B_models_summarized_5_grader=do.call("rbind",B_models_summarized_5_grader)
df_B_models_summarized_2_grader=do.call("rbind",B_models_summarized_2_grader)

df_B_models_summarized_5_grader=
  df_B_models_summarized_5_grader%>%
  mutate(significant=ifelse(p_value<0.01,yes = "**",
                            no = ifelse(p_value<0.05,yes = "*",
                                        no = ifelse(p_value<0.1,yes = "+",no = ""))))%>%
  mutate(direction=ifelse(significant!="",
                          yes = ifelse(t_value>0,yes = "pos",no = "neg"),
                          no = ""))

df_B_models_summarized_2_grader=
  df_B_models_summarized_2_grader%>%
  mutate(significant=ifelse(p_value<0.01,yes = "**",
                            no = ifelse(p_value<0.05,yes = "*",
                                        no = ifelse(p_value<0.1,yes = "+",no = ""))))%>%
  mutate(direction=ifelse(significant!="",
                          yes = ifelse(t_value>0,yes = "pos",no = "neg"),
                          no = ""))

write.csv(x = df_B_models_summarized_5_grader,
          file = paste0(path_output,"\\B_Dist_effect_matrix_result_5_grader.csv"))
write.csv(x = df_B_models_summarized_2_grader,
          file = paste0(path_output,"\\B_Dist_effect_matrix_result_2_grader.csv"))
# (Part3-C) By age group: Test the interested Driving effect(C matrix)------------------
# Derive the distance effect (Near - Far)
df_C_distance_path_all_id=
  df_path_collect_all_id%>%
  select(-X,-region_from)%>%
  filter(matrix_type%in%c("C_Far","C_Medium","C_Near"))%>%
  spread(key=matrix_type,value = value)%>%
  mutate(C_Distance=C_Near-C_Far)

C_models_5_grader=pblapply(interested_region,
                  FUN = function(region){
                    m=lmer(formula = C_Distance~centered_gender_dummy+centered_scan_age_precise+(1|sub_id),
                           data = df_C_distance_path_all_id%>%
                             filter(region_to==region & scan_age_precise>10))
                    return(m)
                  })
C_models_2_grader=pblapply(interested_region,
                           FUN = function(region){
                             m=lmer(formula = C_Distance~centered_gender_dummy+centered_scan_age_precise+(1|sub_id),
                                    data = df_C_distance_path_all_id%>%
                                      filter(region_to==region & scan_age_precise<10))
                             return(m)
                           })
names(C_models_5_grader)=interested_region
names(C_models_2_grader)=interested_region


C_models_summarized_5_grader=lapply(seq_along(C_models_5_grader),
                           FUN = function(model_index){
                             model=C_models_5_grader[[model_index]]
                             model_name=names(C_models_5_grader)[model_index]

                             s=summary(model)
                             coeff=data.frame(s$coefficients)
                             coeff$parameter=rownames(s$coefficients)
                             names(coeff)=c("Estimate","Std_Error","df","t_value","p_value","parameter")
                             coeff$interested_path=model_name
                             return(coeff)
                           })
C_models_summarized_2_grader=lapply(seq_along(C_models_2_grader),
                                    FUN = function(model_index){
                                      model=C_models_2_grader[[model_index]]
                                      model_name=names(C_models_2_grader)[model_index]

                                      s=summary(model)
                                      coeff=data.frame(s$coefficients)
                                      coeff$parameter=rownames(s$coefficients)
                                      names(coeff)=c("Estimate","Std_Error","df","t_value","p_value","parameter")
                                      coeff$interested_path=model_name
                                      return(coeff)
                                    })
df_C_models_summarized_5_grader=do.call("rbind",C_models_summarized_5_grader)
df_C_models_summarized_2_grader=do.call("rbind",C_models_summarized_2_grader)
df_C_models_summarized_5_grader=
  df_C_models_summarized_5_grader%>%
  mutate(significant=ifelse(p_value<0.01,yes = "**",
                            no = ifelse(p_value<0.05,yes = "*",
                                        no = ifelse(p_value<0.1,yes = "+",no = ""))))%>%
  mutate(direction=ifelse(significant!="",
                          yes = ifelse(t_value>0,yes = "pos",no = "neg"),
                          no = ""))
df_C_models_summarized_2_grader=
  df_C_models_summarized_2_grader%>%
  mutate(significant=ifelse(p_value<0.01,yes = "**",
                            no = ifelse(p_value<0.05,yes = "*",
                                        no = ifelse(p_value<0.1,yes = "+",no = ""))))%>%
  mutate(direction=ifelse(significant!="",
                          yes = ifelse(t_value>0,yes = "pos",no = "neg"),
                          no = ""))
write.csv(x = df_C_models_summarized_5_grader,
          file = paste0(path_output,"\\C_Dist_effect_matrix_result_5_grader.csv"))
write.csv(x = df_C_models_summarized_2_grader,
          file = paste0(path_output,"\\C_Dist_effect_matrix_result_2_grader.csv"))

# (Part 4) Plot out scatter plot for all mixed models above---------------------------------------------
library("ggplot2")
library("ggrepel")
library("gridExtra")
# (Part4-A) A matrix scatterplots--------------------------------------
list_collect_plots_A=list()

# Get the y axis limit (to fix the limit for all plots as the overall max and min)
y_max=max(df_path_collect_all_id%>%
            filter(matrix_type=="A")%>%
            pull(value))
y_min=min(df_path_collect_all_id%>%
            filter(matrix_type=="A")%>%
            pull(value))

for(plot in 1:length(interested_path)){
  paths=interested_path[[plot]]
  # Retrive the p values from lmer results (with direction)
  p_labels=
    df_A_models_summarized%>%
    filter(interested_path==paste(paths,collapse = "-"),
           parameter%in%c("(Intercept)","centered_scan_age_precise"))%>%
    mutate(p_label=paste0(round(p_value,3),ifelse(Estimate>0,"+","-")))%>%
    pull(p_label)%>%
    setNames(c("(Intercept)","centered_scan_age_precise"))


  # Plot out the scatter plot with p values
  list_collect_plots_A[[plot]]=
      df_path_collect_all_id%>%
      filter(matrix_type=="A"& region_from==paths[1]& region_to==paths[2])%>%
      mutate(label=paste0(gsub("DF10","",sub_id),"-",gsub("run","",run)))%>%
      # filter(label!="17-4")%>%
      ggplot(aes(x=scan_age_precise,y=value))+
      geom_point(color = "red")+
      geom_smooth(method = "lm")+
      geom_text_repel(aes(label=label),size = 2,segment.size = 0.1)+
      labs(title=paste0("A : ",paste(paths,collapse = "-"),"\n",
                        "p(int.) = ",p_labels["(Intercept)"],"     ",
                        "p(age) = ",p_labels["centered_scan_age_precise"]),
           y="Parameter Estimate (Near-Far)")+
      scale_y_continuous(limits = c(y_min, y_max))
}
g=arrangeGrob(grobs = list_collect_plots_A,ncol = 4)
ggsave(plot = g,
       filename =   paste0(path_output,"\\scatterplot\\A_matrix.pdf"),
       width = 20,height = 30)

# (Part4-B) B matrix scatterplots--------------------------------------
list_collect_plots_B=list()

# Get the y axis limit (to fix the limit for all plots as the overall max and min)
y_max=max(df_B_distance_path_all_id$B_Distance)
y_min=min(df_B_distance_path_all_id$B_Distance)

for(plot in 1:length(interested_path)){
  paths=interested_path[[plot]]
  # Retrive the p values from lmer results (with direction)
  p_labels=
    df_B_models_summarized%>%
    filter(interested_path==paste(paths,collapse = "-"),
           parameter%in%c("(Intercept)","centered_scan_age_precise"))%>%
    mutate(p_label=paste0(round(p_value,3),ifelse(Estimate>0,"+","-")))%>%
    pull(p_label)%>%
    setNames(c("(Intercept)","centered_scan_age_precise"))


  # Plot out the scatter plot with p values
  list_collect_plots_B[[plot]]=
      df_B_distance_path_all_id%>%
        filter(region_from==paths[1]& region_to==paths[2])%>%
        mutate(label=paste0(gsub("DF10","",sub_id),"-",gsub("run","",run)))%>%
        # filter(label!="17-4")%>%
        ggplot(aes(x=scan_age_precise,y=B_Distance))+
        geom_point(color = "red")+
        geom_smooth(method = "lm")+
        geom_text_repel(aes(label=label),size = 2,segment.size = 0.1)+
        labs(title=paste0("B : ",paste(paths,collapse = "-"),"\n",
                          "p(int.) = ",p_labels["(Intercept)"],"     ",
                          "p(age) = ",p_labels["centered_scan_age_precise"]),
                          y="Parameter Estimate (Near-Far)")+
        scale_y_continuous(limits = c(y_min, y_max))
}
g=arrangeGrob(grobs = list_collect_plots_B,ncol = 4)
ggsave(plot = g,
       filename =  paste0(path_output,"\\scatterplot\\B_matrix.pdf"),
       width = 20,height = 30)
# (Part4-C) c matrix scatterplots--------------------------------------
list_collect_plots_C=list()

# Get the y axis limit (to fix the limit for all plots as the overall max and min)
y_max=max(df_C_distance_path_all_id$C_Distance)
y_min=min(df_C_distance_path_all_id$C_Distance)

for(plot in 1:length(interested_region)){
  paths=interested_region[[plot]]
  # Retrive the p values from lmer results (with direction)
  p_labels=
    df_C_models_summarized%>%
    filter(interested_path==paths,
           parameter%in%c("(Intercept)","centered_scan_age_precise"))%>%
    mutate(p_label=paste0(round(p_value,3),ifelse(Estimate>0,"+","-")))%>%
    pull(p_label)%>%
    setNames(c("(Intercept)","centered_scan_age_precise"))


  # Plot out the scatter plot with p values
  list_collect_plots_C[[plot]]=
    df_C_distance_path_all_id%>%
    filter(region_to==paths)%>%
    mutate(label=paste0(gsub("DF10","",sub_id),"-",gsub("run","",run)))%>%
    # filter(label!="17-4")%>%
    ggplot(aes(x=scan_age_precise,y=C_Distance))+
    geom_point(color = "red")+
    geom_smooth(method = "lm")+
    geom_text_repel(aes(label=label),size = 2,segment.size = 0.1)+
    labs(title=paste0("C : ",paste(paths,collapse = "-"),"\n",
                      "p(int.) = ",p_labels["(Intercept)"],"     ",
                      "p(age) = ",p_labels["centered_scan_age_precise"]),
         y="Parameter Estimate (Near-Far)")+
    scale_y_continuous(limits = c(y_min, y_max))
}
g=arrangeGrob(grobs = list_collect_plots_C,ncol = 4)
ggsave(plot = g,
       filename =  paste0(path_output,"\\scatterplot\\C_matrix.pdf"),
       width = 20,height = 15)
