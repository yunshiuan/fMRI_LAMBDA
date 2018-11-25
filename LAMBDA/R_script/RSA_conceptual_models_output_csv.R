#Representational Similarity Analysis(RSA)[Visualize Conceptual model]----------------------------------
#Note that all functions were moved to a separate R script: see the constant FILE_RSA_FUNCTIONS
# Note that I only retain the discrete version since contunuous approach is not used
############################################################################
# Re-define the scale: set w from 0.5 to 1.0 (Trial 19, 2018/6/7)
library("dplyr")
library("tidyr")
library("stringr")
library("reshape2")
library("pbapply")
library("pbmcapply")

# General setup==========================================)======================================================
#Helper functions-------------------------------------
# Load in RSA related functions
# Include dissimilarity definitions and helper functions
FILE_RSA_FUNCTIONS="D:\\GoogleDrive\\Lambda_code\\R_script\\R_functions\\RSA_conceptual_models_functions.R"
# FILE_RSA_FUNCTIONS="/study3/devfracs/DCM_YunShiuan/Lambda_code/R_script/R_functions/RSA_conceptual_models_functions.R"
source(FILE_RSA_FUNCTIONS)
RDM_to_csv=function(df,dis_measure_function,parameters,file_name){
  RDM=df_to_RDM(df = df,
                dist_function = dis_measure_function,
                parameters = parameters)
  write.table(RDM,sep = ",",
              file = file_name,
              col.names = FALSE,row.names = FALSE)
}
#Decalre constants----------------
#Parameters
AMOUNT_TRIALS_TOTAL=18
AMOUNT_FORMAT=3
AMOUNT_TRIALS_EACH=AMOUNT_TRIALS_TOTAL/AMOUNT_FORMAT
RDM_MIN=0
RDM_MAX=1
TRIAL_VERSION="trial_19"
LOCATION_INFORMATION=F # Location matters in some versions but not in others
WEIGHTS=seq(from=0.5,to = 0.99,by=0.005) # The weight range for mixed dist measure and file names
names(WEIGHTS)=as.character(1:99)
#Paths
PATH_ROOT="D:\\Yun-Shiuan_LAMBDA"
PATH_OUTPUT_CSV=file.path(PATH_ROOT,"RSA",TRIAL_VERSION,"Part1_neural_and_conceptual_RDMs","conceptual_model_RDMs")
dir.create(PATH_OUTPUT_CSV,recursive = T)

#Files
FILE_BEH_REF=file.path(PATH_ROOT,"EprimeData_raw","For_RSA","df_for_RDM_construction.csv")

# dataframe for RDM derivation----------------------------
df_trials_discrete_24=preprocess_df_discrete_unique_trails_24(FILE_BEH_REF) # Sort the trials by RSA ID (ID: 1~24)derived from discrete distance (N/M/F)
df_trials_discrete_36=preprocess_df_discrete_repeated_trails_36(FILE_BEH_REF) # Sort the trials by RSA ID (ID: 1~36)derived from discrete magnitude (N/M/F)
df_trials_discrete_18=preprocess_df_discrete_repeated_trails_18(FILE_BEH_REF) # Sort the trials by RSA ID derived from discrete magnitude (N/M/F)
df_trials_target=df_trials_discrete_18

# Discrete magnitude RSA ID==========================================)======================================================
# File names of the output csv===========================================)==========
# Pure format and magnitude measure (parameters unused) and file names-------------------)-----------------------------------------------
list_pure_dist_names=c(#Format models
                      "dist_format_without_location_no_distance",
                      "dist_same_format_pair_no_distance",
                      "dist_null",
                      # "dist_all_the_same",
                       #Pure magnitude models
                       "dist_only_abs_dist_discrete",
                       "dist_only_signed_dist_discrete")
if(LOCATION_INFORMATION){
  list_pure_dist_names=append(list_pure_dist_names,"dist_format_with_location_no_distance")
}
#collect measures into a list
list_pure_dist_functions=mget(list_pure_dist_names)
list_pure_dist_file_names=file.path(PATH_OUTPUT_CSV,paste0("RDM_",list_pure_dist_names,".csv"))
# magnitude hybrid measure (parameters used)-------------------------------------
list_mixed_dist_names=c("dist_format_without_location_abs_dist_discrete",
                        "dist_format_without_location_signed_dist_discrete",
                        "dist_same_format_pair_abs_distance",
                        "dist_same_format_pair_signed_distance"
                         #Partial hybrid models
                        # "dist_format_without_location_LL_abs_dist_discrete",
                        # "dist_format_without_location_LL_signed_dist_discrete",
                        # "dist_same_format_pair_LL_abs_dist_discrete",
                        # "dist_same_format_pair_LL_signed_dist_discrete"
                        )

if(LOCATION_INFORMATION){
  list_mixed_dist_file_names=append(list_mixed_dist_names,
                                    c("dist_format_with_location_abs_dist_discrete",
                                      "dist_format_with_location_signed_dist_discrete"))
}
#generate file names by all possible combinations (mixed dist x WEIGHTS)
list_mixed_dist_file_names=outer(list_mixed_dist_names,names(WEIGHTS),paste,sep="_")
list_mixed_dist_file_names= # To retain the shape of Nfunctions x NWEIGHTS (instead of using paste0())
  apply(list_mixed_dist_file_names,MARGIN = c(1,2),
        FUN = function(x){file.path(PATH_OUTPUT_CSV,paste0("RDM_",x,".csv"))})
#collect measures into a list
list_mixed_dist_functions=mget(list_mixed_dist_names)
# Pure difficulty measure------------------------------------------------------)-------------------------------------------------------------------------------------------------------
#This fundamentally differ from pure magnitude/format measures since different age groups need to be treat separately
list_pure_diff_func_names=c("dist_diff_acc")
list_age_group=c("2_grade","5_grade","adult")
list_pure_diff_file_names=matrix(file.path(PATH_OUTPUT_CSV,
                                           paste0("RDM_",list_pure_diff_func_names,"_",list_age_group,".csv")))
#collect measures into a list
list_pure_diff_functions=mget(list_pure_diff_func_names)
# Difficulty hybrid measure (parameters used)--------------------------------------------------------
list_mixed_diff_names=c("dist_format_without_location_diff_acc",
                        "dist_same_format_pair_diff_acc")
list_mixed_diff_file_names=t(outer(list_mixed_diff_names,list_age_group,paste,sep="_"))
list_mixed_diff_file_names=outer(list_mixed_diff_file_names,names(WEIGHTS),paste,sep="_")
list_mixed_diff_file_names= # To retain the shape of Nagegroup x Nfunc x NWEIGHTS (instead of using paste0())
  apply(list_mixed_diff_file_names,MARGIN = c(1,2,3),
        FUN = function(x){file.path(PATH_OUTPUT_CSV,paste0("RDM_",x,".csv"))})
#collect measures into a list
list_mixed_diff_functions=mget(list_mixed_diff_names)

# Pure log magnitude measure------------------------------------------------------)-------------------------------------------------------------------------------------------------------
#This fundamentally differ from pure magnitude/format measures since there are different versions of magnitude
# list_log_version=c("v1","v2","v3")
list_log_version=c("v2")

list_log_parameters=lapply(setNames(list_log_version,
                                rep("parameters",times=length(list_log_version))),
                      FUN = function(version){
                        parameters=list(log_version=version)})
list_log_dist_func_names=c("dist_only_abs_dist_discrete_log",
                           "dist_only_signed_dist_discrete_log")
list_log_dist_file_names=t(outer(list_log_dist_func_names,list_log_version,paste,sep="_"))
list_log_dist_file_names= # To retain the shape of Nagegroup x Nfunc x NWEIGHTS (instead of using paste0())
  apply(list_log_dist_file_names,MARGIN = c(1,2),
        FUN = function(x){file.path(PATH_OUTPUT_CSV,paste0("RDM_",x,".csv"))})

#collect measures into a list
list_log_dist_functions=mget(list_log_dist_func_names)
# Log magnitude hybrid measure (parameters used)--------------------------------------------------------
list_mixed_dist_log_names=c("dist_format_without_location_abs_dist_discrete_log",
                        "dist_format_without_location_signed_dist_discrete_log",
                        "dist_same_format_pair_abs_distance_log",
                        "dist_same_format_pair_signed_distance_log")#,
                        
                        # #Partial hybrid models
                        # "dist_format_without_location_LL_abs_dist_discrete_log",
                        # "dist_format_without_location_LL_signed_dist_discrete_log",
                        # "dist_same_format_pair_LL_abs_dist_discrete_log",
                        # "dist_same_format_pair_LL_signed_dist_discrete_log"
                        #)
list_mixed_dist_log_file_names=t(outer(list_mixed_dist_log_names,
                                       list_log_version,paste,sep="_"))
list_mixed_dist_log_file_names=outer(list_mixed_dist_log_file_names,names(WEIGHTS),paste,sep="_")
list_mixed_dist_log_file_names= # To retain the shape of Nagegroup x Nfunc x NWEIGHTS (instead of using paste0())
  apply(list_mixed_dist_log_file_names,MARGIN = c(1,2,3),
        FUN = function(x){file.path(PATH_OUTPUT_CSV,paste0("RDM_",x,".csv"))})
#collect measures into a list
list_mixed_dist_log_functions=mget(list_mixed_dist_log_names)

# Output csv=====================================================)------------------------------------
# Pure dist measures--------------------------
pblapply(1:length(list_pure_dist_functions),
       FUN = function(dist_func_index){
         dist_function=list_pure_dist_functions[[dist_func_index]]
         csv_name=list_pure_dist_file_names[[dist_func_index]]

         RDM_to_csv(df = df_trials_target,
                    dis_measure_function = dist_function,
                    file_name =csv_name)
       })
# Magnitude hybrid measure--------------------------
pblapply(1:length(list_mixed_dist_functions),
         FUN = function(dist_func_index){
           dist_function=list_mixed_dist_functions[[dist_func_index]]

           pblapply(1:length(WEIGHTS),
                    FUN = function(weight_index){
                      weight=WEIGHTS[weight_index]
                      csv_name=list_mixed_dist_file_names[dist_func_index,weight_index]
                      RDM_to_csv(df = df_trials_target,
                                 dis_measure_function = dist_function,
                                 file_name =csv_name,
                                 parameters =weight )
                    })
         })
# Pure diff measures-----------
# pblapply(1:length(list_pure_diff_functions),
#          FUN = function(dist_func_index){
#            dist_function=list_pure_diff_functions[[dist_func_index]]
#            pblapply(1:length(list_age_group),
#              FUN=function(age_group_index){
#                csv_name=list_pure_diff_file_names[age_group_index,dist_func_index]
#                parameters=list(age_group=list_age_group[[age_group_index]])
#                
#                RDM_to_csv(df = df_trials_target,
#                           dis_measure_function = dist_function,
#                           file_name = csv_name,
#                           parameters = parameters)
#            })
#         })
# Difficulty hybrid measure--------------------------------------------------------
# pblapply(1:length(list_mixed_diff_functions),
#          FUN = function(dist_func_index){
#            dist_function=list_mixed_diff_functions[[dist_func_index]]
#            pblapply(1:length(list_age_group),
#                     FUN=function(age_group_index){
#                       pblapply(1:length(WEIGHTS),
#                             FUN = function(weight_index){
#                               weight=WEIGHTS[weight_index]
#                               csv_name=list_mixed_diff_file_names[age_group_index,dist_func_index,weight_index]
#                               parameters=list(age_group=list_age_group[[age_group_index]],
#                                               w2=weight)
#                               RDM_to_csv(df = df_trials_target,
#                                          dis_measure_function = dist_function,
#                                          file_name = csv_name,
#                                          parameters = parameters)
#                           })
#                   })
#          })

# Pure log distance measures-----------
pblapply(1:length(list_log_dist_functions),
         FUN = function(dist_func_index){
           dist_function=list_log_dist_functions[[dist_func_index]]
           pblapply(1:length(list_log_version),
                    FUN=function(log_version_index){
                      csv_name=list_log_dist_file_names[log_version_index,dist_func_index]
                      
                      parameters=list_log_parameters[[log_version_index]]
                      #Note that version 3 of log transformation is not supported by the signed distance
                      tryCatch({
                        RDM_to_csv(df = df_trials_target,
                                   dis_measure_function = dist_function,
                                   file_name = csv_name,
                                   parameters = parameters)
                      },error=function(e){
                        print(conditionMessage(e))
                      })
                   
                    })
         })
# Log magnitude hybrid measure--------------------------------------------------------
pblapply(1:length(list_mixed_dist_log_functions),
           FUN = function(dist_func_index){
             dist_function=list_mixed_dist_log_functions[[dist_func_index]]
             pblapply(1:length(list_log_version),
                      FUN=function(log_version_index){
                        pblapply(1:length(WEIGHTS),
                                 FUN = function(weight_index){
                                   weight=WEIGHTS[weight_index]
                                   csv_name=list_mixed_dist_log_file_names[log_version_index,dist_func_index,weight_index]
                                   parameters=list_log_parameters[[log_version_index]]
                                   parameters$w2=weight
                                
                                   #Note that version 3 of log transformation is not supported by the signed distance
                                   tryCatch({
                                     if(!file.exists(csv_name)){#Run if not yet generated
                                       RDM_to_csv(df = df_trials_target,
                                                  dis_measure_function = dist_function,
                                                  file_name = csv_name,
                                                  parameters = parameters)
                                       cat(paste0("Done:","\n",csv_name))
                                     }else{
                                       cat(paste0("Already done:","\n",csv_name))
                                     }
                                   },error=function(e){
                                     cat(paste0(conditionMessage(e),":\n",csv_name))
                                   })
                                 })
                      })
           })
sum(length(list_mixed_dist_log_file_names),length(list_log_dist_file_names),
    length(list_mixed_dist_file_names),length(list_pure_dist_file_names))