#Generate inclusive run indexes to tidy format(column: ID, run) - with manual help
#The old list excludes all subjects with only 3 runs
#The new list includes them back
library("dplyr")
library("stringr")

PATH_INCLUSION="D:\\Yun-Shiuan_LAMBDA\\Run_inclusion_info"
FILE_INCLUSION_RUN_OLD=file.path(PATH_INCLUSION,"inclusive_runs_indexes_new_April_3.csv")
FILE_INCLUSION_RUN_NEW=file.path(PATH_INCLUSION,"inclusive_runs_indexes_new_April_5.csv")
DF1001=c("r1","r2","r6")
#Changing from inclusive_runs_index_April_2 to inclusive_runs_index_April_3
#The manually specified valid runs (subjects with only 3 valid runs)
# DF_INCLUSION_3_RUNS=data.frame(matrix(data=c("df1001","r1",
#                                              "df1001","r2",
#                                              "df1001","r6",
#                                              "df1028","r1",
#                                              "df1028","r5",
#                                              "df1028","r6",
#                                              "df1046","r1",
#                                              "df1046","r2",
#                                              "df1046","r3",
#                                              "df1051","r1",
#                                              "df1051","r3",
#                                              "df1051","r5"),byrow = T,ncol = 2),
#                                stringsAsFactors = F)
#Changing from inclusive_runs_index_April_3 to inclusive_runs_index_April_5
DF_INCLUSION_3_RUNS=data.frame(matrix(data=c("df1018","r1",
                                             "df1018","r2",
                                             "df1018","r4"),byrow = T,ncol = 2),
                               stringsAsFactors = F)

names(DF_INCLUSION_3_RUNS)=c("sub_id","run_num")
run_inclusion_old=read.csv(FILE_INCLUSION_RUN_OLD,header = T,stringsAsFactors = F)
run_inclusion_old%>%
  bind_rows(DF_INCLUSION_3_RUNS)%>%
  arrange(sub_id)%>%
  select(sub_id,run_num)%>%
  write.csv(x=.,file=FILE_INCLUSION_RUN_NEW)
