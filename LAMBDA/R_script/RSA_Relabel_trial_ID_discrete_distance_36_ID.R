#===========================================================================================
# Re-label trial ID [Trial 3: discrete distance with 36 IDs  (Trial 4 in the RSA README.txt)] : combining format(with location) and distance(with sign) information
# Note: the old "trial ID" only accounts for distance,
# but not format(i.e.,trials with same trial ID might differ in terms of format)
# 36 IDs:
# Note Each trial in a run has an ID, so that there'll be 36 IDs
# (note there're only 24 unique trials: 4 formats x 2 locations x 3 distance levels)
# [FL: 6 trials = 3 distance levels x 2 direction]
# [LF: 6 trials = 3 distance levels x 2 direction]
# [LL: 12 trials = 3 distance levels x 2 direction x 2 repetition]
# [FF: 12 trials = 3 distance levels x 2 direction x 2 repetition]

# The motivation for having 36 trials is to ensure beta estimate of each condition is subjected to same amount of noise.
# (Having more repetition in a condition entails a smaller degree of random noise.)

# Note for the ID rule:
# Each id should be unique (thus one need to combine all csv in order to see the overall unique trials)
# Note the order of the id is determined by
# (1) Format in the following order: FF-> FL-> LF-> LL
# (2) Within the same format, order trials by absolute distance (discrete scale: 0~1)
# (3) (make sure to distinguish trials that share the same distance but with opposite sign)
# (4) To account for the repetition within a run(36 trials but only 24 are unique),
# I distinguish the repeated trials (which occur within FF & LL) by absolute continuous distance.

# Note: Continuous distance approach has its problem:
# (1) RSA toolbox only allows all subject run having the same model RDS (every run should have the same RSA IDs)
# (2) Discrete approach could ensure every subject run to have the same RSA IDs, although with the tradeoff that distance information becomes blurred.
#===========================================================================================
library("dplyr")
library("tidyr")
library("stringr")
library("reshape2")
library("pbapply")
# save(list = c("df_all_trials"),file = "D:\\GoogleDrive\\Lambda_code\\R_script\\RData\\all_beh_data_tidy.RData")
#Constants-----------------------
run_data_path="D:\\Yun-Shiuan_LAMBDA\\EprimeData_raw"
valid_run_file="D:\\Yun-Shiuan_LAMBDA\\Run_inclusion_info\\inclusive_runs_indexes.csv"
#Labeling------------------------
#Step 1: Read in run data
temp=list.files(path = run_data_path,pattern = "tidy.csv",full.names = T,recursive = T)
df_list=pblapply(X = temp,
                 FUN = function(csv){
                   df=read.csv(csv,header = T,stringsAsFactors = F)%>%
                     select(-X,-X.1)%>%
                     filter(!is.na(trial_id))
                 })
df_all_trials=do.call(what = "rbind",args = df_list)
# Step 2: Add trial ID  for RSA (combining format(with location) and distance(with sign) information)
# Each id should be unique (thus one need to combine all csv in order to see the overall unique trials)
# Note the order of the id is determined by
# (1) Format in the following order: FF-> FL-> LF->LL
# (2) Within the same format, order trials by absolute discrete distance
# (make sure to distinguish trials that share the same distance but with opposite sign)

a=
  df_all_trials%>%
  filter(trial_id!="NULL")%>%
  mutate(fraccomp_sti_format = toupper(fraccomp_sti_type))%>%
  mutate(fraccomp_sti_format = as.factor(fraccomp_sti_format))%>%
  mutate(fraccomp_sti_dis_discrete_abs=ifelse(fraccomp_sti_dis_type=="Near",yes = 1,
                                              no = ifelse(fraccomp_sti_dis_type=="Medium",yes = 2,
                                                          no = ifelse(fraccomp_sti_dis_type=="Far",yes = 3,no = NA))),
         fraccomp_sti_dis_discrete=(fraccomp_sti_dis>0)*(fraccomp_sti_dis_discrete_abs)+
           (fraccomp_sti_dis<0)*(-fraccomp_sti_dis_discrete_abs))%>%
  arrange(sub_id,run_num,fraccomp_sti_format,abs(fraccomp_sti_dis_discrete),fraccomp_sti_dis_discrete,abs(fraccomp_sti_dis))

# Generate trial ID based on 1)format, 2)absolute discrete distance, 3)signed discrete distance,
#4) (only within a run to assign different ID to repeated trials)absoluate continuous distance
# Note that the dataframe a is already sorted in the desired order
a=a%>%
  mutate(trial_id_RSA_discrete_36=rep(1:36,length.out=nrow(a)))
df_all_trials=a
# save(list = c("df_all_trials"),file = "D:\\GoogleDrive\\Lambda_code\\R_script\\RData\\all_beh_data_tidy.RData")


# Add the "valid run" info
valid_runs=
  read.csv(valid_run_file,header = T,stringsAsFactors = F)%>%
  select(-X)%>%
  mutate(sub.id=gsub(pattern="^df","",sub.id))
valid_runs_collapsed=paste0(valid_runs$sub.id,"-",valid_runs$run.num)

df_all_trials=
  df_all_trials%>%
  mutate(id_run=paste0(sub_id,"-",run_num))%>%
  mutate(valid_run=ifelse(id_run%in%valid_runs_collapsed,yes = 1,no = 0),
         trial_id=as.numeric(trial_id))%>%
  select(-id_run)%>%
  arrange(sub_id,run_num,trial_id_RSA_discrete_36)

# See the properties of each unique trial------------------------
trial_info=
  a%>%
  select(trial_id_RSA_discrete_36,fraccomp_sti_format,fraccomp_sti_dis_discrete)%>%
  unique()

trial_info_with_sub_id=
  a%>%
  select(trial_id_RSA_discrete_36,fraccomp_sti_format,fraccomp_sti_dis_discrete,sub_id)%>%
  unique()

collect_sub_id_unique_amount=c()
collect_sub_id_list_collapsed=c()

for(trial_id in 1:length(unique(trial_info_with_sub_id$trial_id_RSA_discrete_36))){
  sub_id_list=
    trial_info_with_sub_id%>%
    filter(trial_id_RSA_discrete_36==trial_id)%>%
    select(sub_id)%>%
    pull()%>%
    gsub("^10","",x = .)
  sub_id_unique_amount=length(sub_id_list)
  sub_id_list_collapsed=paste(sub_id_list,collapse = ";")

  collect_sub_id_unique_amount=append(collect_sub_id_unique_amount,sub_id_unique_amount)
  collect_sub_id_list_collapsed=append(collect_sub_id_list_collapsed,sub_id_list_collapsed)
}

trial_info$sub_id_unique_amount=collect_sub_id_unique_amount
trial_info$sub_id_list_collapsed=collect_sub_id_list_collapsed
# Check if every subject run has 36 ID
a%>%
  group_by(sub_id,run_num)%>%
  summarise(n_trials_unique=length(unique(trial_id_RSA_discrete_36)))%>%
  View()
# Check the stimulus properties of each RSA ID
df_all_trials%>%
  group_by(RSA_ID=trial_id_RSA_discrete_36)%>%
  summarise(num_total_runs=n(),
            format=unique(fraccomp_sti_format),
            dis_discrete=mean(fraccomp_sti_dis_discrete),
            dis_discrete_abs=mean(fraccomp_sti_dis_discrete_abs),
            dis_abs=mean(abs(fraccomp_sti_dis))
            )%>%View()

# Write back to the tidy.csv----------------------------------
temp=list.files(path = run_data_path,pattern = "tidy.csv",full.names = T,recursive = T)
df_list=pblapply(X = temp,
                 FUN = function(csv){
                   csv_sub_id=as.integer(unlist(str_extract_all(string = csv,pattern = "(?<=/df)\\d+(?=_Run)")))
                   csv_run_num=paste0("r",unlist(str_extract_all(string = csv,pattern = "(?<=Run)\\d+(?=_tidy)")))


                   df=read.csv(csv,header = T,stringsAsFactors = F)%>%
                     mutate(trial_id=as.numeric(trial_id))
                   # Skip if already joined the trial_id_RSA info.
                   if("trial_id_RSA_discrete_36" %in% names(df)){
                     print(paste0("Already done: id:",csv_sub_id," ; run: ", csv_run_num))
                   }else{
                     df_all_trials_selected=
                       df_all_trials%>%
                         select(sub_id,run_num,trial_id_RSA_discrete_36,trial_id,fraccomp_sti_dis_discrete)%>%
                         filter(sub_id==csv_sub_id,run_num==csv_run_num)

                     df=df%>%
                       select(-X.1)%>%
                       left_join(df_all_trials_selected,by = c("sub_id","run_num","trial_id","fraccomp_sti_dis_discrete"))

                     write.csv(df,file = csv)
                   }
                 })

