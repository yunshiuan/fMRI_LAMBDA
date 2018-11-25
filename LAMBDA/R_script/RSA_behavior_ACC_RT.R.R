#RSA: Check if the behavior data seems valid(i.e., Adults > 5th Grade > 2nd Grade)
library(dplyr)
library(tidyr)
library(stringr)
library(pbapply)
#Declare contsant----------------------------------------------
#Path------
PATH_ROOT="D:\\Yun-Shiuan_LAMBDA"
PATH_BEH_FILES_CHILD=file.path(PATH_ROOT,"EprimeData_raw")
PATH_BEH_FILES_ADULT=file.path(PATH_ROOT,"Adult","EprimeData_raw")
#File------
#Valid subject run info
FILE_VALID_RUNS_CHILD=file.path(PATH_ROOT,"Run_inclusion_info","inclusive_runs_indexes_new_April_5.csv")
FILE_VALID_RUNS_ADULT=file.path(PATH_ROOT,"Adult","Run_inclusion_info","inclusive_runs_indexes.csv")
FILE_DEMOGRAPHIC_CHILD=file.path(PATH_ROOT,"demographic","tidy_demographic.csv")
FILE_DEMOGRAPHIC_ADULT=file.path(PATH_ROOT,"Adult","demographic","tidy_demographic.csv")
FILE_RSA_FUNCTIONS="D:\\GoogleDrive\\Lambda_code\\R_script\\R_functions\\RSA_conceptual_models_functions.R"
# FILE_BEH_REFERENCE=file.path(PATH_ROOT,"EprimeData_raw","For_RSA","df_for_RDM_construction.csv")

#Parameters-------------------------
RT_CUTOFF=300 #The mimimun valid RT
PAR_FORMAT_COLOR=c("#1a1aff","#ff4dff","#ff0000")

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
                      trial_id_RSA_discrete_18,sub_id,
                      fraccomp_sti_dis_type,fraccomp_sti_type)%>%
               #Retrieve the format and magnitude information
               mutate(format=tolower(fraccomp_sti_type),
                      format=case_when(format=="f - f" ~ "FF",
                                       format %in% c("f - l","l - f") ~ "CN",
                                       format=="l - l" ~ "LL"))%>%
               rename(magnitude=fraccomp_sti_dis_type)%>%
               #Mark those missing responses
               mutate(fraccomp_resp_RT=ifelse(fraccomp_resp_RT==0,yes = NA,no = fraccomp_resp_RT),
                      fraccomp_resp_acc=ifelse(fraccomp_resp_RT==0,yes = NA,no = fraccomp_resp_acc))
           })
#Collapse all dfs into a big df and derive averaged RT and Acc (group by Age Group x Magnitude x Format)-----------------------------------
df_all_sub=
  do.call("rbind",list_all_sub_run)
df_all_sub=
  df_all_sub%>%
  left_join(df_demographic%>%
              mutate(sub_id=as.numeric(unlist(str_extract_all(sub_id,"\\d+")))),
            by = "sub_id")
#Add ACC
df_beh_mean_all_sub=
  df_all_sub%>%
  group_by(format,magnitude,grade)%>%
  summarise(ACC_mean=mean(fraccomp_resp_acc,na.rm=T),
            n=n(),
            ACC_se=sd(fraccomp_resp_acc,na.rm=T)/sqrt(n))
#Add in valid RT
df_beh_mean_all_sub=
  df_all_sub%>%
  filter(fraccomp_resp_acc==1,
         fraccomp_resp_RT>=RT_CUTOFF)%>%
  group_by(format,magnitude,grade)%>%
  summarise(RT_mean=mean(fraccomp_resp_RT,na.rm = T),
            n=n(),
            RT_se=sd(fraccomp_resp_RT,na.rm=T)/sqrt(n))%>%
  select(-n)%>%
  right_join(df_beh_mean_all_sub,by = c("format","magnitude","grade"))
#Factorize the variable to faciliate visualization
df_beh_mean_all_sub=
  df_beh_mean_all_sub%>%
  as.data.frame()%>%
  mutate(grade=factor(grade,levels = c("2","5","Adult"),
                      labels = c("Second Graders","Fifth Graders","Adults")),
         magnitude=factor(magnitude,levels=c("Near","Medium","Far")),
         format=factor(format,levels=c("FF","CN","LL")))
#Plot the average ACC & RT by age group plot-----------------------------------------------------------------------------------
g_RT=
  df_beh_mean_all_sub%>%
    ggplot(aes(x=magnitude,y=RT_mean,group=format,color=format))+
    geom_line()+
    scale_color_manual(values = PAR_FORMAT_COLOR)+
    geom_point()+
    geom_errorbar(aes(ymax=RT_mean+RT_se,ymin=RT_mean-RT_se),width=0.2)+
    facet_grid(.~grade)+
    theme_bw()
g_ACC=
  df_beh_mean_all_sub%>%
    ggplot(aes(x=magnitude,ACC_mean,group=format,color=format))+
    geom_line()+
    scale_color_manual(values = PAR_FORMAT_COLOR)+
    geom_point()+
    geom_errorbar(aes(ymax=ACC_mean+ACC_se,ymin=ACC_mean-ACC_se),width=0.2)+
    facet_grid(.~grade)+
    theme_bw()
grid.arrange(g_RT,g_ACC)

#Collapse all dfs into a big df and derive averaged RT and Acc (group by Age Group x Magnitude x Format x Sub_id)-----------------------------------

#Add ACC
df_beh_ind=
  df_all_sub%>%
  group_by(format,magnitude,grade,sub_id)%>%
  summarise(ACC_mean=mean(fraccomp_resp_acc,na.rm=T),
            n=n(),
            ACC_se=sd(fraccomp_resp_acc,na.rm=T)/sqrt(n))
#Add in valid RT
df_beh_ind=
  df_all_sub%>%
  filter(fraccomp_resp_acc==1,
         fraccomp_resp_RT>=RT_CUTOFF)%>%
  group_by(format,magnitude,grade,sub_id)%>%
  summarise(RT_mean=mean(fraccomp_resp_RT,na.rm = T),
            n=n(),
            RT_se=sd(fraccomp_resp_RT,na.rm=T)/sqrt(n))%>%
  select(-n)%>%
  right_join(df_beh_ind,by = c("format","magnitude","grade","sub_id"))
#Factorize the variable to faciliate visualization
df_beh_ind=
  df_beh_ind%>%
  as.data.frame()%>%
  mutate(grade=factor(grade,levels = c("2","5","Adult"),
                      labels = c("Second Graders","Fifth Graders","Adults")),
         magnitude=factor(magnitude,levels=c("Near","Medium","Far")),
         format=factor(format,levels=c("FF","CN","LL")),
         sub_id=str_extract(string=sub_id,pattern="(?<=^(10)|(3))\\d+"))

#Plot the individual ACC & RT by age group plot-----------------------------------------------------------------------------------
#Factorize the variable to faciliate visualization
pos <- position_jitter(width = 0.5)
#devtools::install_github("slowkow/ggrepel") #Allow pos argument for ggrepel
g_RT=
  df_beh_ind%>%
  ggplot(aes(x=magnitude,y=RT_mean,group=format,color=format))+
  # geom_line()+
  scale_color_manual(values = PAR_FORMAT_COLOR)+
  geom_point()+
  # geom_errorbar(aes(ymax=RT_mean+RT_se,ymin=RT_mean-RT_se),width=0.2)+
  facet_grid(format~grade)+
  theme_bw()+
  # geom_text(aes(label=sub_id),position=pos)+
  # geom_jitter()+
  ggrepel::geom_text_repel(aes(label=sub_id))
g_ACC=
  df_beh_ind%>%
  ggplot(aes(x=magnitude,y=ACC_mean,group=format,color=format))+
  # geom_line()+
  scale_color_manual(values = PAR_FORMAT_COLOR)+
  geom_point()+
  # geom_errorbar(aes(ymax=RT_mean+RT_se,ymin=RT_mean-RT_se),width=0.2)+
  geom_jitter()+
  facet_grid(.~grade)+
  theme_bw()
  # ggrepel::geom_text_repel(aes(label=sub_id))
grid.arrange(g_RT,g_ACC)