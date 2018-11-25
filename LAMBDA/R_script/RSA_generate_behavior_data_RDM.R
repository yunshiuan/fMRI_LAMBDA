# Generate difficulty RDM based on RT and accuracy measures-----------------------------
#Purpose:
#(1)Add average(across all valid subject runs) RT/ACC as trial's features and add it to the "reference Run_tidy.csv".
# Make sure to add it back so that I could still utilize the functions in "RSA_conceptual_models_functions.R", and
# keep using the scripts of visualization and csv generation.
#Note:
#(1)Features to be included
#(1-1)2nd grader's RT/ACC
#(1-2)5th grader's RT/ACC
#(1-3)Adult's RT/ACC
#(1-4)Collapsed RT/ACC
library(dplyr)
library(tidyr)
library(stringr)
library(pbapply)
#Declare contsant----------------------------------------------
#Path------
PATH_BEH_FILES_CHILD="/study3/devfracs/DCM_YunShiuan/EprimeData_raw"
PATH_BEH_FILES_ADULT="/study3/devfracs/DCM_YunShiuan/Adult/EprimeData_raw"
#File------
#Valid subject run info
FILE_VALID_RUNS_CHILD="/study3/devfracs/DCM_YunShiuan/Run_inclusion_info/inclusive_runs_indexes.csv"
FILE_VALID_RUNS_ADULT="/study3/devfracs/DCM_YunShiuan/Adult/Run_inclusion_info/inclusive_runs_indexes.csv"
FILE_DEMOGRAPHIC_CHILD="/study3/devfracs/DCM_YunShiuan/demographic/tidy_demographic.csv"
FILE_DEMOGRAPHIC_ADULT="/study3/devfracs/DCM_YunShiuan/Adult/demographic/tidy_demographic.csv"
FILE_RSA_FUNCTIONS="/study3/devfracs/DCM_YunShiuan/Lambda_code/R_script/R_functions/RSA_conceptual_models_functions.R"
FILE_BEH_REFERENCE="/study3/devfracs/DCM_YunShiuan/EprimeData_raw/For_RSA/df_for_RDM_construction.csv"
#Parameters-------------------------
RT_CUTOFF=300 #The mimimun valid RT
#Generate averaged(across subject runs) difficulty-related trial features ================)==========================
#Set up and preprocess data----------------------------------------------------------
source(FILE_RSA_FUNCTIONS)
# Get valid subject run info.
df_valid_sub_run_child=
  read.csv(FILE_VALID_RUNS_CHILD,stringsAsFactors = F,header = T)%>%
  select(-X)
df_valid_sub_run_adult=
  read.csv(FILE_VALID_RUNS_ADULT,stringsAsFactors = F,header = T)%>%
  select(-X)
df_valid_sub_run=
  df_valid_sub_run_child%>%
  bind_rows(df_valid_sub_run_adult)
# Get all dfs for runs
df_csv_child=data.frame(csv=list.files(path = PATH_BEH_FILES_CHILD,pattern = "tidy.csv",
                        recursive=T,full.names=T),stringsAsFactors = F)
df_csv_adult=data.frame(csv=list.files(path = PATH_BEH_FILES_ADULT,pattern = "tidy.csv",
                        recursive=T,full.names=T),stringsAsFactors = F)
df_csv=
  df_csv_child%>%
  bind_rows(df_csv_adult)%>%
  mutate(sub_id=unlist(str_extract_all(string=csv,pattern="(df|XFC)\\d+(?=_Run)")),
         run_num=paste0("r",str_extract_all(string=csv,pattern="(?<=Run)\\d(?=_tidy.csv)")))%>%
  right_join(df_valid_sub_run,by = c("sub_id","run_num"))
#Add in age group information to the "df_csv" ----------------------------------
df_demographic_child=
  read.csv(FILE_DEMOGRAPHIC_CHILD,header = T,stringsAsFactors = F)%>%
  mutate(sub_id=tolower(sub_id))%>%#To match up with the folder names
  select(sub_id,grade)%>%
  mutate(grade=as.character(grade))
df_demographic_adult=
  read.csv(FILE_DEMOGRAPHIC_ADULT,header = T,stringsAsFactors = F)%>%
  select(sub_id,grade)%>%
  mutate(grade=as.character(grade))
df_demographic=
  df_demographic_child%>%
  bind_rows(df_demographic_adult)
df_csv=
  df_csv%>%
  left_join(df_demographic,by = "sub_id")
#Load in and process behavior data-------------------------------
list_all_sub_run=
  pblapply(X = df_csv$csv,
           FUN = function(csv){
             df=read.csv(csv,stringsAsFactors = F,header = T)
             df%>%
               #Filter out the null trials
               filter(!is.na(trial_id_RSA_discrete_18),
                      !is.null(trial_id_RSA_discrete_18))%>%
               select(fraccomp_resp_RT,fraccomp_resp_acc,
                      trial_id_RSA_discrete_18,sub_id)%>%
               #Mark those missing responses
               mutate(fraccomp_resp_RT=ifelse(fraccomp_resp_RT==0,yes = NA,no = fraccomp_resp_RT),
                      fraccomp_resp_acc=ifelse(fraccomp_resp_RT==0,yes = NA,no = fraccomp_resp_acc))
     })
#Collapse all dfs into a big df and derive averaged RT and Acc-----------------------------------
df_all_sub=
  do.call("rbind",list_all_sub_run)
df_all_sub=
  df_all_sub%>%
    left_join(df_demographic%>%
               mutate(sub_id=as.numeric(unlist(str_extract_all(sub_id,"\\d+")))),
              by = "sub_id")
df_beh_mean_all_sub=
  df_all_sub%>%
    group_by(trial_id_RSA_discrete_18)%>%
    summarise(ACC_mean_all=mean(fraccomp_resp_acc,na.rm=T))
df_beh_mean_all_sub=#Add in valid RT
  df_all_sub%>%
  filter(fraccomp_resp_acc==1,
         fraccomp_resp_RT>=RT_CUTOFF)%>%
  group_by(trial_id_RSA_discrete_18)%>%
  summarise(RT_mean_all=mean(fraccomp_resp_RT,na.rm = T))%>%
  right_join(df_beh_mean_all_sub,by = "trial_id_RSA_discrete_18")

df_beh_mean_2_grade=
  df_all_sub%>%
  filter(grade=="2")%>%
  group_by(trial_id_RSA_discrete_18)%>%
  summarise(ACC_mean_2_grade=mean(fraccomp_resp_acc,na.rm=T))
df_beh_mean_2_grade=
  df_all_sub%>%
  filter(grade=="2")%>%
  filter(fraccomp_resp_acc==1,
         fraccomp_resp_RT>=RT_CUTOFF)%>%
  group_by(trial_id_RSA_discrete_18)%>%
  summarise(RT_mean_2_grade=mean(fraccomp_resp_RT,na.rm = T))%>%
  right_join(df_beh_mean_2_grade,by = "trial_id_RSA_discrete_18")

df_beh_mean_5_grade=
  df_all_sub%>%
  filter(grade=="5")%>%
  group_by(trial_id_RSA_discrete_18)%>%
  summarise(ACC_mean_5_grade=mean(fraccomp_resp_acc,na.rm=T))
df_beh_mean_5_grade=
  df_all_sub%>%
  filter(grade=="5")%>%
  filter(fraccomp_resp_acc==1,
         fraccomp_resp_RT>=RT_CUTOFF)%>%
  group_by(trial_id_RSA_discrete_18)%>%
  summarise(RT_mean_5_grade=mean(fraccomp_resp_RT,na.rm = T))%>%
  right_join(df_beh_mean_5_grade,by = "trial_id_RSA_discrete_18")

df_beh_mean_adult=
  df_all_sub%>%
  filter(grade=="Adult")%>%
  group_by(trial_id_RSA_discrete_18)%>%
  summarise(ACC_mean_adult=mean(fraccomp_resp_acc,na.rm=T))
df_beh_mean_adult=
  df_all_sub%>%
  filter(grade=="Adult")%>%
  filter(fraccomp_resp_acc==1,
         fraccomp_resp_RT>=RT_CUTOFF)%>%
  group_by(trial_id_RSA_discrete_18)%>%
  summarise(RT_mean_adult=mean(fraccomp_resp_RT,na.rm = T))%>%
  right_join(df_beh_mean_adult,by = "trial_id_RSA_discrete_18")
#Add the feature back to the reference tidy csv--------------------------------------------------
df_for_RDM_construction=
  read.csv(FILE_BEH_REFERENCE,header = T,stringsAsFactors = F)%>%
  left_join(df_beh_mean_all_sub,by = "trial_id_RSA_discrete_18")%>%
  left_join(df_beh_mean_2_grade,by = "trial_id_RSA_discrete_18")%>%
  left_join(df_beh_mean_5_grade,by = "trial_id_RSA_discrete_18")%>%
  left_join(df_beh_mean_adult,by = "trial_id_RSA_discrete_18")%>%
  mutate_at(.vars=vars(matches("(RT|ACC)_mean")),
            .fun=funs(round(.,4)))%>%
  write.csv(x=.,file=FILE_BEH_REFERENCE)


# #(Incomplete) Generate ifficulty RDMs per each subject (average across runs)================)==========================
# list_all_sub_run=
#   pblapply(X = unique(df_csv$sub_id),
#            FUN = function(sub) {
#              sub_csv=
#                df_csv%>%
#                filter(sub_id==sub)%>%
#                select(csv)%>%
#                pull()
#
#              #List of all the valid csv of the given subject (so that it could be collapsed afterwards)
#              list_df_sub=
#                lapply(
#                  X = sub_csv,
#                  FUN = function(csv){
#                    df=read.csv(csv,stringsAsFactors = F,header = T)
#                    df%>%
#                      #Filter out the null trials
#                      filter(!is.na(trial_id_RSA_discrete_18),
#                             !is.null(trial_id_RSA_discrete_18))%>%
#                      #Mark those missing responses
#                      mutate(fraccomp_resp_RT=ifelse(fraccomp_resp_RT==0,yes = NA,no = fraccomp_resp_RT),
#                             fraccomp_resp_acc=ifelse(fraccomp_resp_RT==0,yes = NA,no = fraccomp_resp_acc))
#                  })
#              df_sub=do.call("rbind",list_df_sub) #Collapse the runs' df into a single df per subject
#            })
