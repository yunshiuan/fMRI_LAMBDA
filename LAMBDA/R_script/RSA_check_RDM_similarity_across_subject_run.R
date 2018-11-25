#Representational Similarity Analysis(RSA) Check RDM similarity across subject run----------------------------------
library("dplyr")
library("tidyr")
library("stringr")
library("reshape2")
source("D:\\GoogleDrive\\Lambda_code\\R_script\\mclapply.hack.R")
#[Part1]Generate each subject run's RDMs====================================================
#Preparation================================================================================
# Functions for RSA - RDM generation----------
load(file_RSA_RDM_functions)
#Constants----------------------
dist_functions_list=c(dist_format_without_location_no_distance,dist_format_with_location_no_distance,
                      dist_only_abs_dist,dist_format_without_location_abs_dist,dist_format_with_location_abs_dist,
                      dist_only_signed_dist,dist_format_without_location_signed_dist,dist_format_with_location_signed_dist)
dist_functions_names=c("dist_format_without_location_no_distance","dist_format_with_location_no_distance",
                       "dist_only_abs_dist","dist_format_without_location_abs_dist","dist_format_with_location_abs_dist",
                       "dist_only_signed_dist","dist_format_without_location_signed_dist","dist_format_with_location_signed_dist")
weight_dist=0.5 #Onle compare the mixed model with weight=0.5 to reduce comparison amount
file_baseline_df="D:\\Yun-Shiuan_LAMBDA\\EprimeData_raw\\DevFracs1002\\df1002_Run2_tidy.csv"
file_valid_run="D:\\GoogleDrive\\Lambda_code\\R_script\\RData\\valid_runs.RData"
path_tidy_csv="D:\\Yun-Shiuan_LAMBDA\\EprimeData_raw"
# file_RSA_RDM_functions="/Users/vimchiz/googledrive/Lambda_code/R_script/RData/RSA_RDM_functions.RData"
file_RSA_RDM_functions="D:\\GoogleDrive\\Lambda_code\\R_script\\RData\\RSA_RDM_functions.RData"
# file_collect_all_RDMs_all_sub_run="/Users/vimchiz/googledrive/Lambda_code/R_script/RData/collect_all_RDMs_all_sub_run.RData"
file_collect_all_RDMs_all_sub_run_valid="D:\\GoogleDrive\\Lambda_code\\R_script\\RData\\collect_all_RDMs_all_sub_run_valid.RData"
file_baseline_RDM_matrices="D:\\GoogleDrive\\Lambda_code\\R_script\\RData\\baseline_RDM_matrices.RData"
#Helper functions----------------

# Function for preprocessing data frame into RSA required format-----------
preprocess_df=function(csv){
  processed_df=
    read.csv(csv)%>%
    filter(!is.na(trial_id))%>%
    select(trial_id_RSA_discrete,fraccomp_sti_dis_discrete,
           fraccomp_sti_type,fraccomp_sti_dis_type,
           fraccomp_sti_value_left,fraccomp_sti_value_right,fraccomp_sti_dis)%>%
    mutate(fraccomp_sti_dis_abs=abs(fraccomp_sti_dis),
           fraccomp_sti_dis_discrete_abs=abs(fraccomp_sti_dis_discrete))%>%
    mutate(F_left=grepl("^(f|F)",fraccomp_sti_type),
           F_right=grepl("(f|F)$",fraccomp_sti_type),
           L_left=grepl("^(l|L)",fraccomp_sti_type),
           L_right=grepl("(l|L)$",fraccomp_sti_type))%>%
    mutate(F_amount=F_left+F_right,
           L_amount=L_left+L_right)%>%
    arrange(trial_id_RSA_discrete)%>%
    mutate(fraccomp_sti_dis_abs_scaled=scale_RSA(fraccomp_sti_dis_abs,0,2),
           fraccomp_sti_dis_scaled=scale_RSA(fraccomp_sti_dis,0,2),
           fraccomp_sti_dis_discrete_scaled=scale_RSA(fraccomp_sti_dis_discrete,0,2),
           fraccomp_sti_dis_discrete_abs_scaled=scale_RSA(fraccomp_sti_dis_discrete_abs,0,2))
  return(processed_df)
}
# Function for generate all RDMs for a given subject run
generate_RDMs=function(df){
  RDM_matrices=list()
  for (dist_fun_index in 1:length(dist_functions_list)){
    RDM_matrix=df_to_RDM(df = df,
                         dist_function = dist_functions_list[[dist_fun_index]],
                         parameters = weight_dist)
    RDM_matrices[[dist_fun_index]]=RDM_matrix
  }
  return(RDM_matrices)
}


#Valid runs---------------------
load(file_valid_run) #valid_runs
valid_runs=
  valid_runs%>%
    mutate(sub_id_run_num=paste0("df",sub.id,"_",gsub(x = run.num,pattern="r",replacement="Run")))
#Baseline matrix-----------------
#Set df1002 run1 as the baseline matrix (to be compared to)
load(file_baseline_RDM_matrices)#baseline_RDM_matrices
# df_baseline=
#   preprocess_df(file_baseline_df)
# baseline_RDM_matrices=
#   generate_RDMs(df_baseline)
# save(baseline_RDM_matrices,file = file_baseline_RDM_matrices)

#Generate each subject run's RDMs===============================================================
temp=list.files(path = path_tidy_csv,pattern = "tidy.csv",recursive = T,full.names = T)
temp_inclusion_index=unlist(str_extract_all(string = temp,pattern = "df\\d+_Run\\d(?=_tidy)")) %in% valid_runs$sub_id_run_num
temp=temp[temp_inclusion_index]
# Only run when needed
# collect_all_RDMs_al_sub_run=mclapply.hack(temp,
#                               function(csv){
#                                 df=preprocess_df(csv)
#                                 print(unlist(str_extract_all(string = csv,pattern = "df\\d+_Run\\d+(?=_tidy)")))
#                                 list_matrices=generate_RDMs(df)
#                                 return(list_matrices)
#                               })

#[Part2]Quantify the RDM difference across subject run=================================================
load(file_collect_all_RDMs_all_sub_run_valid)#collect_all_RDMs_all_sub_run_valid
names(collect_all_RDMs_all_sub_run_valid)=valid_runs$sub_id_run_num
#The absolute difference between matrices----------------------------------------
diff_RDM_abs_diff=function(RDM_baseline,RDM_other){
  num=sum(abs(RDM_other-RDM_baseline))
  dem=sum(abs(RDM_baseline))
  return((num/dem)*100)
}
#The spatial correlation between matrices-----------------------------------------
diff_RDM_spatial_cor=function(RDM_baseline,RDM_other){
  cor=cor.test(RDM_baseline,RDM_other,method = "pearson")
  return(list(r=cor$estimate,p=cor$p.value))
}
#Compute matrix differences-------------------------------------------------------
df_collect_RDM_differences=data.frame(matrix(ncol = 5,nrow = 0))
names(df_collect_RDM_differences)=c("abs_diff","corr_r","corr_p","sub_run","RDM_type")
for (RDM_type in 1:length(dist_functions_list)){
  RDM_baseline=baseline_RDM_matrices[[RDM_type]]
  for(RDM_sub_run_index in 1:length(collect_all_RDMs_all_sub_run_valid)){
    RDM_other=collect_all_RDMs_all_sub_run_valid[[RDM_sub_run_index]][[RDM_type]]

    abs_diff=diff_RDM_abs_diff(RDM_baseline = RDM_baseline,RDM_other = RDM_other)
    spatial_cor=diff_RDM_spatial_cor(RDM_baseline = RDM_baseline,RDM_other = RDM_other)


    df_temp=data.frame(abs_diff,spatial_cor$r,spatial_cor$p,
                       names(collect_all_RDMs_all_sub_run_valid)[RDM_sub_run_index],dist_functions_names[RDM_type])
    names(df_temp)=names(df_collect_RDM_differences)
    df_collect_RDM_differences=
      rbind(df_collect_RDM_differences,df_temp)
    print(paste0("sub_run_id: ",RDM_sub_run_index," RDM_type:",RDM_type))
  }
}

#Conclusion:----------------------------------------------------------------------
tidy_df_collect_RDM_differences=
  df_collect_RDM_differences%>%
  mutate(
    abs_diff=round(abs_diff,2),
    corr_r=round(corr_r,2),
    corr_p=round(corr_p,3)
  )
View(tidy_df_collect_RDM_differences)
#(1)Some RDMs differ quite a lot (r=.032~1.0,med=0.34, esp. for signed dist.; Nor)
tidy_df_collect_RDM_differences%>%
  filter(RDM_type=="dist_only_signed_dist")%>%
  pull(corr_r)%>%
  summary()
#(2)Abs distance is less of a problem(r=0.85~1.0,med=0.97)
tidy_df_collect_RDM_differences%>%
  filter(RDM_type=="dist_only_abs_dist")%>%
  pull(corr_r)%>%
  summary()
#(3)Perfect match: Format with location (no distance) (r=1)
tidy_df_collect_RDM_differences%>%
  filter(RDM_type=="dist_format_with_location_no_distance")%>%
  pull(corr_r)%>%
  summary()
#(4)Perfect match: Format without location (no distance) (r=1)
tidy_df_collect_RDM_differences%>%
  filter(RDM_type=="dist_format_without_location_no_distance")%>%
  pull(corr_r)%>%
  summary()
# Next step:
#Relabel the RSA ID to ensure every run has the same 1~36.
#Convert the continuous distance(0~1) to descrete distance(N/M/F)
#Remain to be checked:
##(1)if every run has the same 1~36.
##(2)after relabel, the RDMs should have perfect match for all runs
