#Generate the list of valid runs for the adults---------------------------------
#Transform inclusive run indexes to tidy format(column: ID, run)
library(dplyr)
library(tidyr)
library(stringr)
#Declare constants-----------------------
FILE_RAW="D:\\Yun-Shiuan_LAMBDA\\Adult\\Run_inclusion_info\\FracComp_Data_ParticipantTracking.csv"
FILE_TIDY_OUTPUT="D:\\Yun-Shiuan_LAMBDA\\Adult\\Run_inclusion_info\\inclusive_runs_indexes.csv"
NUM_MAX_RUN=6
#Clean the data into desired format-----------------------
df=read.csv(FILE_RAW,header = T,stringsAsFactors = F)
df%>%
  select(sub_id=ID,
         exclusion_decision=Exclusion.Decision,
         exclusion_note=Exclusion.Notes)%>%
  mutate(exclusion_amount=str_extract_all(exclusion_decision,pattern = "(^\\d(?=/))|(EXCLUDE)"),
         excluded_run=str_extract_all(exclusion_note,"(?<=Run ).*(?= Excluded)"))%>%
  mutate(excluded_run=gsub(x=excluded_run,pattern="( |&)",replacement=""))%>%
  mutate(sub_id=gsub(x=sub_id,pattern="_",replacement=""))%>%
  # Initialize all runs as true (and then being excluded below)
  mutate(r1=1,r2=1,r3=1,r4=1,r5=1,r6=1)%>%
  #Deal with those fully excluded
  mutate(r1=ifelse(exclusion_amount=="EXCLUDE",yes = 0,no = r1),
         r2=ifelse(exclusion_amount=="EXCLUDE",yes = 0,no = r2),
         r3=ifelse(exclusion_amount=="EXCLUDE",yes = 0,no = r3),
         r4=ifelse(exclusion_amount=="EXCLUDE",yes = 0,no = r4),
         r5=ifelse(exclusion_amount=="EXCLUDE",yes = 0,no = r5),
         r6=ifelse(exclusion_amount=="EXCLUDE",yes = 0,no = r6)
  )%>%
#Deal with those who have some excluded runs
  mutate(r1=ifelse(grepl(excluded_run,pattern = "1"),yes = 0,no = r1),
         r2=ifelse(grepl(excluded_run,pattern = "2"),yes = 0,no = r2),
         r3=ifelse(grepl(excluded_run,pattern = "3"),yes = 0,no = r3),
         r4=ifelse(grepl(excluded_run,pattern = "4"),yes = 0,no = r4),
         r5=ifelse(grepl(excluded_run,pattern = "5"),yes = 0,no = r5),
         r6=ifelse(grepl(excluded_run,pattern = "6"),yes = 0,no = r6))%>%
  select(sub_id,matches("^r\\d$"))%>%
  gather(run_num,include,matches("^r\\d$"))%>%
  arrange(sub_id)%>%
  filter(include==1)%>%
  select(-include)%>%
  write.csv(x=.,
            file = FILE_TIDY_OUTPUT)
