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

# (Part1) Load in demographic information---------------------------------
df_demo=read.csv(file = "D:\\Yun-Shiuan_LAMBDA\\demographic\\tidy_demographic.csv",stringsAsFactors = F)

# (Part2) LMER: Intrinsic connections, Age effect, and Distance Effect------------------
# Read in the data frame (with all path coefficients)
df_path_collect_all_id=read.csv(file = "D:\\Yun-Shiuan_LAMBDA\\DCM\\Part3_second_level_LMER\\collect_all_path_coefficients.csv",
                                stringsAsFactors = F)
# 
# df_path_collect_all_id=
#   df_path_collect_all_id%>%
#   left_join(df_demo%>%select(scan_age_precise,gender,sub_id),
#             by = "sub_id")
# write.csv(x = df_path_collect_all_id,
#           file = "D:\\Yun-Shiuan_LAMBDA\\DCM\\Part3_second_level_LMER\\collect_all_path_coefficients.csv")

interested_path=list(c("L_V1","L_IPS"),c("L_V1","R_IPS"),
                     c("R_V1","L_IPS"),c("R_V1","R_IPS"),
                     c("L_IPS","L_SMA"),c("L_IPS","R_SMA"),
                     c("R_IPS","L_SMA"),c("R_IPS","R_SMA"),
                     c("L_SMA","L_M1"),c("R_SMA","L_M1"),
                     c("R_DLPFC","L_V1"),c("R_DLPFC","R_V1"),
                     c("R_DLPFC","L_IPS"),c("R_DLPFC","R_IPS"),
                     c("R_DLPFC","L_SMA"),c("R_DLPFC","R_SMA"),
                     c("R_DLPFC","L_M1"))
interested_region=c("L_V1","R_V1","L_IPS","R_IPS",
                    "L_SMA","R_SMA","L_M1","R_DLPFC")
# (Part2-A) Test the interested intrinsic connecitons (-and age-related effect) (A matrix)------------------

A_models=pblapply(interested_path,
           FUN = function(paths){
             m=lmer(formula = value~gender+scan_age_precise+(1|sub_id),
                    data = df_path_collect_all_id%>%
                      filter(matrix_type=="A"& region_from==paths[1]& region_to==paths[2]))
             return(m)
         })
names(A_models)=unlist(lapply(interested_path,
                              FUN = function(paths){paste0(paths[1],"-",paths[2])}))
A_models_summarized=lapply(seq_along(A_models),
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
          file = "D:\\Yun-Shiuan_LAMBDA\\DCM\\Part3_second_level_LMER\\A_matrix_result.csv")
# (Part2-B) Test the modulating effect - Distance effect and Age effect (A matrix)------------------
# Derive the distance effect (Near - Far)
df_B_distance_path_all_id=
  df_path_collect_all_id%>%
    select(-X)%>%
    filter(matrix_type%in%c("B_Far","B_Medium","B_Near"))%>%
    spread(key=matrix_type,value = value)%>%
    mutate(B_Distance=B_Near-B_Far)

B_models=pblapply(interested_path,
           FUN = function(paths){
             m=lmer(formula = B_Distance~gender+scan_age_precise+(1|sub_id),
                    data = df_B_distance_path_all_id%>%
                      filter(region_from==paths[1]& region_to==paths[2]))
             return(m)
          })
names(B_models)=unlist(lapply(interested_path,
                              FUN = function(paths){paste0(paths[1],"-",paths[2])}))
B_models_summarized=lapply(seq_along(B_models),
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
          file = "D:\\Yun-Shiuan_LAMBDA\\DCM\\Part3_second_level_LMER\\B_Dist_effect_matrix_result.csv")
# (Part2-C) Test the Driving effect - Distance effect and Age effect (A matrix)------------------

# Derive the distance effect (Near - Far)
df_C_distance_path_all_id=
  df_path_collect_all_id%>%
  select(-X,-region_from)%>%
  filter(matrix_type%in%c("C_Far","C_Medium","C_Near"))%>%
  spread(key=matrix_type,value = value)%>%
  mutate(C_Distance=C_Near-C_Far)

C_models=pblapply(interested_region,
                  FUN = function(region){
                    m=lmer(formula = C_Distance~gender+scan_age_precise+(1|sub_id),
                           data = df_C_distance_path_all_id%>%
                             filter(region_to==region))
                    return(m)
                  })
names(C_models)=interested_region

C_models_summarized=lapply(seq_along(C_models),
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
          file = "D:\\Yun-Shiuan_LAMBDA\\DCM\\Part3_second_level_LMER\\C_Dist_effect_matrix_result.csv")


# (Part3) Split LMER results by age group
# (Part3-A) By age group: Test the interested intrinsic connecitons(A matrix)------------------
A_models_5_grader=pblapply(interested_path,
                  FUN = function(paths){
                    m=lmer(formula = value~gender+scan_age_precise+(1|sub_id),
                           data = df_path_collect_all_id%>%
                             filter(matrix_type=="A"& region_from==paths[1]& region_to==paths[2] & scan_age_precise>10))
                    return(m)
                  })

A_models_2_grader=pblapply(interested_path,
                           FUN = function(paths){
                             m=lmer(formula = value~gender+scan_age_precise+(1|sub_id),
                                    data = df_path_collect_all_id%>%
                                      filter(matrix_type=="A"& region_from==paths[1]& region_to==paths[2] & scan_age_precise<=10))
                             return(m)
                           })
names(A_models_5_grader)=unlist(lapply(interested_path,
                              FUN = function(paths){paste0(paths[1],"-",paths[2])}))
names(A_models_2_grader)=unlist(lapply(interested_path,
                              UN = function(paths){paste0(paths[1],"-",paths[2])}))

A_models_summarized_5_grader=lapply(seq_along(A_models_5_grader),
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
A_models_summarized_2_grader=lapply(seq_along(A_models_2_grader),
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
          file = "D:\\Yun-Shiuan_LAMBDA\\DCM\\Part3_second_level_LMER\\A_matrix_result_5_grader.csv")

write.csv(x = df_A_models_summarized_2_grader,
          file = "D:\\Yun-Shiuan_LAMBDA\\DCM\\Part3_second_level_LMER\\A_matrix_result_2_grader.csv")