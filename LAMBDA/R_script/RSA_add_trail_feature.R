#Add features of conditions for constructing RDMs------------------------------
#E.g., log-transformed discrete distance
# Note that only absolute distance could be log-transformed, since log-transformation only allows positive input.

#For log-transformed discrete distance, there're different version:
#[Version 1] Based on the assumption that value is represented logarithmically while the distance linearlly
# i.e., log(value 1)-log(value 2) [Fechner's law]
#1) Compute the log(value 1)-log(value 2) as the "logarithmic distance"
#2) Average the log distance within each RSA ID batch (across all 6 runs)
#(Note: Raw distance and new log distance is no longer monotonically correspondent.
# This may require re-defining the RSA ID, i.e., use the new log distance to lump trials together for t value estimation)
#(Note: This version allow both absolute and signed distance)
#Absolute distance [Assume 1:3 and 3:1 is the same]: abs(log(value 1)-log(value 2)) 
#Signed distance [Assume 1:3 and 3:1 is different]: log(value 1)-log(value 2)

#[Version 2] Based on the assumption that value is represented  linearlly while the distance logarithmically
# i.e., log(abs(value 1-value 2))
#1) Compute the log(value 1-value 2) as the "logarithmic distance"
#2) Average the log distance within each RSA ID batch (across all 6 runs)
#(Note: This version allow both absolute and signed distance)
#Absolute distance [Assume 1:3 and 3:1 is the same]: log(abs(value 1-value 2))
#Signed distance [Assume 1:3 and 3:1 is different]: sign(rescale_min_to_zero(log(abs(value 1-value 2))))

#[Version 3] Based on the assumption that value is represented logarithmically and so does the distance
# i.e., log[log(value 1)-log(value 2)]
#1) Compute the [log(value 1)-log(value 2)] as the "logarithmic distance"
#2) Average the log distance within each RSA ID batch (across all 6 runs)
#(Note: Raw distance and new log distance is no longer monotonically correspondent.
# This may require re-defining the RSA ID, i.e., use the new log distance to lump trials together for t value estimation)

library(dplyr)
library(tidyr)
library(stringr)
library(pbapply)
#Declare contsant=================================================
#Path------
PATH_ROOT="D:\\Yun-Shiuan_LAMBDA"
PATH_RSA_ROOT=file.path(PATH_ROOT,"RSA")
PATH_BEH_FILES_CHILD=file.path(PATH_ROOT,"EprimeData_raw")
PATH_BEH_FILES_ADULT=file.path(PATH_ROOT,"Adult","EprimeData_raw")
# PATH_BEH_DATA=file.path(PATH_ROOT,"EprimeData_raw","DevFracs1002") #This could be any subject, as long as there''re 6 runs
#File------
#Valid subject run info
FILE_RSA_FUNCTIONS="D:\\GoogleDrive\\Lambda_code\\R_script\\R_functions\\RSA_conceptual_models_functions.R"
FILE_BEH_REFERENCE=file.path(PATH_ROOT,"EprimeData_raw","For_RSA","df_for_RDM_construction.csv") # The file where new var is added
FILE_VALID_RUNS_CHILD=file.path(PATH_ROOT,"Run_inclusion_info","inclusive_runs_indexes.csv")
FILE_VALID_RUNS_ADULT=file.path(PATH_ROOT,"Adult","Run_inclusion_info","inclusive_runs_indexes.csv")
FILE_DEMOGRAPHIC_CHILD=file.path(PATH_ROOT,"demographic","tidy_demographic.csv")
FILE_DEMOGRAPHIC_ADULT=file.path(PATH_ROOT,"Adult","demographic","tidy_demographic.csv")

#Parameters-------------------------
VAR_RSA_ID_VERSION="trial_id_RSA_discrete_18"
#Helper funtion=====================================================
source(FILE_RSA_FUNCTIONS)
#Get valid subject run info------------------------------------
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

#Generate log-transformed discrete distance for each RSA ID (18)---------------------------)---------------------------

# list_beh_csv=list.files(path = PATH_BEH_DATA,
#                         pattern = "Run\\d_tidy\\.csv",
#                         full.names = T)
list_beh_csv=df_csv_child
list_new_log_measures=
  lapply(1:length(list_beh_csv$csv),
       FUN = function(csv_index){
         csv=list_beh_csv$csv[csv_index]
         df=read.csv(csv,header = T,stringsAsFactors = F)
         df=
           df%>%
             #Filter out the null trials
             filter(!is.na(trial_id_RSA_discrete_18),
                    !is.null(trial_id_RSA_discrete_18))%>%
             select(!!VAR_RSA_ID_VERSION,run_num,fraccomp_sti_dis_type,
                    fraccomp_sti_dis,fraccomp_sti_value_left,fraccomp_sti_value_right)%>%
             mutate(fraccomp_sti_dis_abs=abs(fraccomp_sti_dis),#Transfer to abs. distance before log-transformation
                    sign=ifelse(fraccomp_sti_dis>0,yes = 1,no = -1),#The sign information
                    #Version 1
                    fraccomp_sti_dis_abs_logv1=abs(log(fraccomp_sti_value_right)-log(fraccomp_sti_value_left)),
                    fraccomp_sti_dis_signed_logv1=log(fraccomp_sti_value_right)-log(fraccomp_sti_value_left),
                    #Version 2
                    fraccomp_sti_dis_abs_logv2=log(fraccomp_sti_dis_abs),
                    fraccomp_sti_dis_signed_logv2=scale_RSA(log(fraccomp_sti_dis_abs),new_min = 0,new_max = 1)*sign,
                    #version 3
                    fraccomp_sti_dis_abs_logv3=log(abs(log(fraccomp_sti_value_right)-log(fraccomp_sti_value_left)))
                    )%>%
             select(-sign)
         df$sub_id=list_beh_csv$sub_id[csv_index]
         return(df)
       })
# Average across runs
df_new_log_measures_merged=do.call(what = "rbind.data.frame",
                                   args = list_new_log_measures)
df_new_log_measures_merged_avg=
  df_new_log_measures_merged%>%
    group_by(trial_id_RSA_discrete_18)%>%
    summarise(fraccomp_sti_dis_abs_logv1_mean=mean(fraccomp_sti_dis_abs_logv1),
              fraccomp_sti_dis_signed_logv1_mean=mean(fraccomp_sti_dis_signed_logv1),
              
              fraccomp_sti_dis_abs_logv2_mean=mean(fraccomp_sti_dis_abs_logv2),
              fraccomp_sti_dis_signed_logv2_mean=mean(fraccomp_sti_dis_signed_logv2),
              
              fraccomp_sti_dis_abs_logv3_mean=mean(fraccomp_sti_dis_abs_logv3),
              fraccomp_sti_dis_abs_mean=mean(fraccomp_sti_dis_abs),#For referece only
              
              fraccomp_sti_dis_abs_logv1_se=sd(fraccomp_sti_dis_abs_logv1),
              fraccomp_sti_dis_signed_logv1_se=sd(fraccomp_sti_dis_signed_logv1),
              
              fraccomp_sti_dis_abs_logv2_se=sd(fraccomp_sti_dis_abs_logv2),
              fraccomp_sti_dis_signed_logv2_se=sd(fraccomp_sti_dis_signed_logv2),
              
              fraccomp_sti_dis_abs_logv3_se=sd(fraccomp_sti_dis_abs_logv3),
              fraccomp_sti_dis_abs_se=sd(fraccomp_sti_dis_abs),
              n_rt=sqrt(12) #Each ID repeated 12 times (6 runs)
              # n_rt=sqrt(n())# SQRT of Number of trials (for deriving se)
              )%>%
    as.data.frame()%>% #Disable tibble display
    mutate_at(.vars=vars(matches("_se$")),
               .funs=funs(./n_rt))

df_new_log_measures_merged_avg=
  df_new_log_measures_merged_avg%>%
  mutate(dist_discrete_abs=rep(c(1:3),each=2,times=3),
         format=rep(c("FF","CN","LL"),each=6),
         sign=rep(c("L>R","R>L"),times=9))

#Write back to the reference beh csv
df_beh_reference=read.csv(FILE_BEH_REFERENCE,header = T,stringsAsFactors = F)
df_beh_reference%>%
  select(-matches("(^X.\\d*$)|(^X$)|log|abs_mean"))%>%
  left_join(df_new_log_measures_merged_avg%>%
              select(!!VAR_RSA_ID_VERSION,
                     fraccomp_sti_dis_abs_logv1_mean,
                     fraccomp_sti_dis_signed_logv1_mean,
                     fraccomp_sti_dis_abs_logv2_mean,
                     fraccomp_sti_dis_signed_logv2_mean,
                     fraccomp_sti_dis_abs_logv3_mean,
                     fraccomp_sti_dis_abs_mean),
            by = VAR_RSA_ID_VERSION)%>%
  write.csv(.,file = FILE_BEH_REFERENCE)
#Visualization of th log-transformation===============================================
# Scatter plot (raw abs. continuous distance vs. new abs. log continuous distance)------------------------------------
#Note that raw distance and new log distance is no longer monotonically correspondent.
ggplot_log_types=
  df_new_log_measures_merged%>%
  gather(log_type,log_dist_abs,matches("fraccomp_sti_dis_abs_logv\\d"))%>%
  rename(raw_dist_abs=fraccomp_sti_dis_abs)%>%
  ggplot(aes(x=raw_dist_abs,y=log_dist_abs,color=log_type))+
  geom_point()+
  theme_classic()+
  labs(x="Raw Abs. Difference Magnitude",
       y="Log-Transformed Abs. Difference Magnitude")+
  scale_color_discrete(name="Log-Transformation Type",
                       labels=c("Log fraction value, \nLinear difference mag. (v1)\n log(v1)-log(v2)",
                                "Linear fracion value, \nLog difference mag. (v2)\n log( |v1-v2| )",
                                "Log fracion value, \nLog difference mag. (v3) \n log( |log(v1)-log(v2)| )"))+
  theme(axis.title = element_text(size=16),
        axis.text= element_text(size=14),
        legend.title = element_text(size=13,face = "bold"),
        legend.text = element_text(size=12),
        legend.key.height =unit(3,"line"))
#Add the lines to indicate the raw discrete distance levels
df_cutoff_lines=
  df_new_log_measures_merged%>%
    filter(!is.na(trial_id_RSA_discrete_18),
           !is.null(trial_id_RSA_discrete_18))%>%
    select(fraccomp_sti_dis_type,fraccomp_sti_dis)%>%
    mutate(fraccomp_sti_dis_abs=abs(fraccomp_sti_dis))%>%
    group_by(fraccomp_sti_dis_type)%>%
    summarise(min=min(fraccomp_sti_dis_abs),
              max=max(fraccomp_sti_dis_abs))
line_N_M=mean(c(as.numeric(df_cutoff_lines[df_cutoff_lines$fraccomp_sti_dis_type=="Near","max"]),
                as.numeric(df_cutoff_lines[df_cutoff_lines$fraccomp_sti_dis_type=="Medium","min"])))
line_M_F=mean(c(as.numeric(df_cutoff_lines[df_cutoff_lines$fraccomp_sti_dis_type=="Medium","max"]),
                as.numeric(df_cutoff_lines[df_cutoff_lines$fraccomp_sti_dis_type=="Far","min"])))
ggplot_log_types+
  geom_vline(xintercept = c(line_N_M,line_M_F),linetype="dashed")
  # ggrepel::geom_text_repel(aes(label=trial_id_RSA_discrete_18))
# Scatter plot (raw signed continuous distance vs. new log continuous distance)------------------------------------
#Note that raw distance and new log distance is no longer monotonically correspondent.
ggplot_log_types=
  df_new_log_measures_merged%>%
  mutate(raw_dist_signed=fraccomp_sti_dis)%>%
  gather(log_type,log_dist,matches("signed_logv\\d|fraccomp_sti_dis$"))%>%
  mutate(log_type=factor(log_type,levels = c("fraccomp_sti_dis_signed_logv1",
                                             "fraccomp_sti_dis_signed_logv2",
                                             "fraccomp_sti_dis")))%>%
  ggplot(aes(x=raw_dist_signed,y=log_dist,color=log_type))+
  geom_point()+
  theme_classic()+
  labs(x="Raw Abs. Difference Magnitude",
       y="Log-Transformed Abs. Difference Magnitude")+
  scale_color_manual(name="Log-Transformation Type",
                     labels=c(
                              "Log fraction value, \nLinear difference mag. (v1)\n Mag = log(v1)-log(v2)",
                              "Linear fracion value, \nLog difference mag. (v2)\n Mag = sign*shifted(log(|v1-v2|))",
                              "Linear fracion value, \nLinear difference mag. (raw) \n Mag = (v1-v2)"
                              ),
                     values = c("#F8766D","#00BA38","black"))+
  theme(axis.title = element_text(size=16),
        axis.text= element_text(size=14),
        legend.title = element_text(size=13,face = "bold"),
        legend.text = element_text(size=12),
        legend.key.height =unit(3,"line"))
#Add the lines to indicate the raw discrete distance levels
df_cutoff_lines=
  df_new_log_measures_merged%>%
  filter(!is.na(trial_id_RSA_discrete_18),
         !is.null(trial_id_RSA_discrete_18))%>%
  select(fraccomp_sti_dis_type,fraccomp_sti_dis)%>%
  mutate(fraccomp_sti_dis_abs=abs(fraccomp_sti_dis))%>%
  group_by(fraccomp_sti_dis_type)%>%
  summarise(min=min(fraccomp_sti_dis_abs),
            max=max(fraccomp_sti_dis_abs))
line_N_M=mean(c(as.numeric(df_cutoff_lines[df_cutoff_lines$fraccomp_sti_dis_type=="Near","max"]),
                as.numeric(df_cutoff_lines[df_cutoff_lines$fraccomp_sti_dis_type=="Medium","min"])))
line_M_F=mean(c(as.numeric(df_cutoff_lines[df_cutoff_lines$fraccomp_sti_dis_type=="Medium","max"]),
                as.numeric(df_cutoff_lines[df_cutoff_lines$fraccomp_sti_dis_type=="Far","min"])))
ggplot_log_types+
  geom_vline(xintercept = c(-line_M_F,-line_N_M,0,line_N_M,line_M_F),linetype="dashed")
# ggrepel::geom_text_repel(aes(label=trial_id_RSA_discrete_18))

# Averaged point plot (raw abs. discrete distance vs. new abs. log discrete distance)------------------------------------
pattern_interested_variables_mean="fraccomp_sti_dis_abs_logv\\d_mean|fraccomp_sti_dis_abs_mean"
pattern_interested_variables_se="fraccomp_sti_dis_abs_logv\\d_se|fraccomp_sti_dis_abs_se"

df_log_type_mean=
  df_new_log_measures_merged_avg%>%
  # filter(trial_id_RSA_discrete_18%%2==0)%>% #Even: (+: R>L)
  gather(log_type_mean,log_dist_abs_mean,matches(pattern_interested_variables_mean))
df_log_type_se=
  df_new_log_measures_merged_avg%>%
  # filter(trial_id_RSA_discrete_18%%2==0)%>% #Even: (+: R>L)
  gather(log_type_se,log_dist_abs_se,matches(pattern_interested_variables_se))
df_log_type=
  df_log_type_mean%>%
  bind_cols(df_log_type_se%>%select(log_type_se,log_dist_abs_se))

df_log_type%>%
  rename(raw_dist_abs=dist_discrete_abs)%>%
  ggplot(aes(x=raw_dist_abs,y=log_dist_abs_mean,fill=log_type_mean))+
  geom_bar(stat = "identity",position=position_dodge())+
  geom_errorbar(aes(ymax=log_dist_abs_mean+1.96*log_dist_abs_se,
                    ymin=log_dist_abs_mean-1.96*log_dist_abs_se),position=position_dodge())+
  theme_classic()+
  labs(x="Absolute Difference Magnitude Discrete Level",
       y="Averaged Log-Transformed Abs. Difference Magnitude")+
  scale_fill_discrete(name="Log-Transformation Type",
                      labels=c("Log fraction value, \nLinear difference mag. (v1)\n Mag.=log(v1)-log(v2)",
                               "Linear fracion value, \nLog difference mag. (v2)\n Mag.=log( |v1-v2| )",
                               "Log fracion value, \nLog difference mag. (v3) \n Mag.=log( |log(v1)-log(v2)| )",
                               "Linear fracion value, \nLinear difference mag. (v2)\n Mag.=|v1-v2| "))+
  scale_x_discrete(limits=c(1,2,3),
                   labels=c("Near","Medium","Far"))+
  theme(axis.title = element_text(size=16),
        axis.text= element_text(size=12),
        legend.title = element_text(size=13,face = "bold"),
        legend.text = element_text(size=12),
        legend.key.height =unit(3,"line"))+
  facet_wrap(sign~format)

# Rescale version: Averaged point plot (raw abs. discrete distance vs. new abs. log discrete distance)-------------
#To compare the degree of linearity
pattern_interested_variables_mean="fraccomp_sti_dis_abs_logv\\d_mean|fraccomp_sti_dis_abs_mean"

df_new_log_measures_merged_avg%>%
  rename(raw_dist_abs=dist_discrete_abs)%>%
  mutate_at(.vars=vars(matches(pattern_interested_variables_mean)),
            .funs=funs(scale_RSA(.,new_min = 1,new_max = 3)))%>%
  gather(log_type_mean,log_dist_abs_mean,matches(pattern_interested_variables_mean))%>%
  ggplot(aes(x=raw_dist_abs,y=log_dist_abs_mean,fill=log_type_mean))+
  geom_bar(stat = "identity",position=position_dodge())+
  theme_classic()+
  labs(x="Absolute Difference Magnitude Discrete Level",
       y="Rescaled Averaged Log-Transformed Abs. Difference Magnitude")+
  scale_fill_discrete(name="Log-Transformation Type",
                      labels=c("Log fraction value, \nLinear difference mag. (v1)\n log(v1)-log(v2)",
                               "Linear fracion value, \nLog difference mag. (v2)\n log( |v1-v2| )",
                               "Log fracion value, \nLog difference mag. (v3) \n log( |log(v1)-log(v2)| )",
                               "Linear fracion value, \nLinear difference mag. (v2)\n Mag.=|v1-v2| "))+
  scale_x_discrete(limits=c(1,2,3),
                   labels=c("Near","Medium","Far"))+
  theme(axis.title = element_text(size=16),
        axis.text= element_text(size=12),
        legend.title = element_text(size=13,face = "bold"),
        legend.text = element_text(size=12),
        legend.key.height =unit(3,"line"))+
  facet_wrap(sign~format)+
  geom_hline(yintercept =2 ,linetype="dashed")

#TO BE DEBBUGED=======================================================================================
#Some messed up trials for adults (left_value-right_value!=fraccomp_sti_dis)-----------------------------
#For now, only use the children's beh files to derive log distance as in above.
df_new_log_measures_merged%>%
  mutate(fraccomp_sti_dis_lr=round(fraccomp_sti_value_right-fraccomp_sti_value_left,2),
         descrepancy=fraccomp_sti_dis_lr-fraccomp_sti_dis,
         consistency=(abs(descrepancy)<=0.0051))%>%
  filter(!consistency)%>%
  write.csv(.,file = "D:\\Yun-Shiuan_LAMBDA\\to_debug\\Inconsistency_between_dist_and_pair_value.csv")

df_new_log_measures_merged%>%
  filter(fraccomp_sti_dis_abs_logv3==-Inf)%>%
  write.csv(.,file = "D:\\Yun-Shiuan_LAMBDA\\to_debug\\Inconsistency_between_dist_and_pair_value_extreme_case.csv")
