#Generate inclusive run indexes to tidy format(column: ID, run) - From txt to csv
library("dplyr")
library("stringr")
PATH_INCLUSION="/study3/devfracs/DCM_YunShiuan/Run_inclusion_info"
FILE_INCLUSION_RUN_RAW=file.path(PATH_INCLUSION,"DevFracs_5thGraders_BinNote-New_RFX_TAL_8mm_smoothed.txt")
FILE_INCLUSION_RUN_PROCESSED=file.path(PATH_INCLUSION,"inclusive_runs_indexes_new_5_grade_April_2.csv")

run_inclusion=read.table(file = FILE_INCLUSION_RUN_RAW,
                         skip = 11,blank.lines.skip = T,stringsAsFactors = F)

run_inclusion%>%
  select(V1)%>%
  mutate(sub_id=paste0("df",
                      str_extract_all(string=V1,pattern="(?<=DevFracs)\\d+(?=_)")),
         run_num=paste0("r",str_extract_all(string=V1,pattern="(?<=_Run)\\d+(?=_)")))%>%
  select(-V1)%>%
  write.csv(x=.,file=FILE_INCLUSION_RUN_PROCESSED)

#Merge the two age group's inclusion list---------------------------------------------------------
FILE_INCLUSION_RUN_PROCESSED_2_GRADE=file.path(PATH_INCLUSION,"inclusive_runs_indexes_new_2_grade_April_2.csv")
FILE_INCLUSION_RUN_PROCESSED_5_GRADE=file.path(PATH_INCLUSION,"inclusive_runs_indexes_new_5_grade_April_2.csv")
FILE_INCLUSION_RUN_PROCESSED_MERGED=file.path(PATH_INCLUSION,"inclusive_runs_indexes_new_April_2.csv")

run_inclusion_2_grade=read.csv(FILE_INCLUSION_RUN_PROCESSED_2_GRADE,header = T,stringsAsFactors = F)
run_inclusion_5_grade=read.csv(FILE_INCLUSION_RUN_PROCESSED_5_GRADE,header = T,stringsAsFactors = F)

run_inclusion_2_grade%>%
  bind_rows(run_inclusion_5_grade)%>%
  arrange(sub_id)%>%
  select(sub_id,run_num)%>%
  write.csv(x=.,file=FILE_INCLUSION_RUN_PROCESSED_MERGED)
