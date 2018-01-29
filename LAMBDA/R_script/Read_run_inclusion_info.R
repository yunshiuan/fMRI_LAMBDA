
#(PART2) Transform inclusive run indexes to tidy format(column: ID, run)

setwd("D:\\Users\ychuang26\\GoogleDrive\\Lambda_code\\R_script")
library("dplyr")
library("tidyr")
library("stringr")
path="D:\\Yun-Shiuan_LAMBDA\\Run_inclusion_info\\runs_included_xfc_child.txt"
run_inclusion=read.table(file = path,skip = 11,blank.lines.skip = T,stringsAsFactors = F)

run_inclusion%>%
  select(V1)%>%
  mutate(sub.id=paste0("df",
                      str_extract_all(string=V1,pattern="(?<=DevFracs)\\d+(?=_)")),
         run.num=paste0("r",str_extract_all(string=V1,pattern="(?<=_Run)\\d+(?=_)")))%>%
  select(-V1)%>%
  write.csv(x=.,file="D:\\Yun-Shiuan_LAMBDA\\Run_inclusion_info\\inclusive_runs_indexes.csv")
         
                  