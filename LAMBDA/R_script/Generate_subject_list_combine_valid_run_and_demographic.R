#Generate the valid subject list with demographic information-----------------------------------
library("dplyr")
#Constants----------------------------------------------------------------------------------------
#Path
PATH_ROOT="D:\\Yun-Shiuan_LAMBDA"
PATH_OUTPUT_CSV=file.path(PATH_ROOT,"Run_inclusion_info")
  
#subject list
FILE_SUBJECT_LIST_CHILD=file.path(PATH_ROOT,"Run_inclusion_info","inclusive_runs_indexes_new_June_11.csv")
FILE_SUBJECT_LIST_ADULT=file.path(PATH_ROOT,"Adult","Run_inclusion_info","inclusive_runs_indexes.csv")
FILE_DEMOGRAPHIC_CHILD=file.path(PATH_ROOT,"demographic","tidy_demographic.csv")
FILE_DEMOGRAPHIC_ADULT=file.path(PATH_ROOT,"Adult","demographic","tidy_demographic.csv")
#Read in demographic information----------------------------------------------------------------
df_subject_list_child=read.csv(FILE_SUBJECT_LIST_CHILD,stringsAsFactors = F,header = T)
df_subject_list_adult=read.csv(FILE_SUBJECT_LIST_ADULT,stringsAsFactors = F,header = T)
df_demographic_child=read.csv(FILE_DEMOGRAPHIC_CHILD,stringsAsFactors = F,header = T)
df_demographic_adult=read.csv(FILE_DEMOGRAPHIC_ADULT,stringsAsFactors = F,header = T)
df_demographic=
  df_demographic_child%>%
  select(X,sub_id,grade)%>%
  rbind.data.frame(df_demographic_adult)%>%
  select(-X)
df_subject_list=
  df_subject_list_child%>%
  rbind.data.frame(df_subject_list_adult)%>%
  select(-X)
df_subject_list_with_demo=
  df_subject_list%>%
  mutate(sub_id=toupper(sub_id))%>%
  left_join(df_demographic,by = c("sub_id"))
df_subject_list_with_demo=
  df_subject_list_with_demo%>%
  mutate(sub_id=ifelse(grepl(x = sub_id,
                             pattern="DF"),
                       yes=tolower(sub_id),
                       no = sub_id),
         grade_character=case_when(grade=="2" ~ "grade_2",
                                   grade=="5" ~ "grade_5",
                                   grade=="Adult" ~ "adult"),
         grade_numeric=ifelse(grade=="Adult",
                              yes = 0,no = grade))%>%
  select(-grade)
write.csv(x = df_subject_list_with_demo,
          file = file.path(PATH_OUTPUT_CSV,
                           "inclusive_runs_indexes_with_demographic_June_11.csv"))