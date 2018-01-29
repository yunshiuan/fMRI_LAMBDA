#(PART3) Extract interested demographics for second-level analysis (e.g., gender, age, )
setwd("D:\\GoogleDrive\\Lambda_code\\R_script")
library("dplyr")
library("tidyr")
library("stringr")
path="D:\\Yun-Shiuan_LAMBDA\\demographic\\demographics_scans.csv"
df=read.csv(path,header = T,stringsAsFactors = F)
df%>%
  filter(Imaging.Code!="")%>% # filter out those without scans
  rename(birth_date=Birth.date,
         test_date=Test.date,
         scan_date=Scan.Date,
         sub_id=Imaging.Code,
         age=Age..years.,
         gender=Sex,
         ethinicity=Ethnicity,
         race=Race,
         event_name=Event.Name
         )%>%
  mutate_at(vars(birth_date,test_date,scan_date),
            funs(as.Date(.,format='%m/%d/%Y')))%>%
  mutate(scan_age_precise=round((as.numeric(scan_date-birth_date)/365.25),2),
         gender_dummy=ifelse(gender=="Male",1,
                       ifelse(gender=="Female",-1,"Unknown")),
         grade=ifelse(age>9,yes = 5,no = 2))%>%
  select(-Subject.ID,-event_name,-Form,-X,-X.1,-X.2,-X.3)%>%
  write.csv(x = .,file = "D:\\Yun-Shiuan_LAMBDA\\demographic\\tidy_demographic.csv")
  
       