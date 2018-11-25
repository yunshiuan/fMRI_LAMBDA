#Transform run info. from Eprime raw xlxs to a tidy format for MATLAB to read

library("dplyr")
library("tidyr")
library("stringr")
library("rio")
library("pbapply")
#Declare constants----------------------------------
PATH_BEH_DATA="D:\\Yun-Shiuan_LAMBDA\\Adult\\EprimeData_raw"
STRING_SUB_PREFIX="XFC3"
#The pattern of run number in the "run_id" variable (this varies across children and adults)
PATTERN_RUN_NUM="(?<=^0)\\d+(?=-|_)"#Children : "(?<=^r)\\d+(?=_xfc)"
NUM_TRIALS=38
SUBJECT_LIST=dir(path = PATH_BEH_DATA,
               pattern=paste0("^",STRING_SUB_PREFIX))
# Note that this variable is named differently between adults and children
# For children: congruent=ComponentCongruency
VAR_NAME_CONGRUENCY="NumCongruency"
VAR_ADJUST_DIST_DIRECTION=1 # -1 for children(original: dist= left-right) and +1 for adults(original: dist= right-left)
#Helper functions
as.frac.numeric=function(x)
              {sapply(x,function(x)eval(parse(text=x)))}#"1/2" to 0.5; works for vector

lapply_with_error <- function(X,FUN,...){
lapply(X, function(x, ...) tryCatch(FUN(x, ...),
                                    error=function(e)
                                      paste0(x," has error:",e)))
}#define lapply with try-catch

#Transform csv xlsx to csv and then output my desired format of run_info.csv file
pblapply(SUBJECT_LIST,function(sub){

    directory_sub=paste0(PATH_BEH_DATA,"\\",sub)
    xlsx=list.files(path = directory_sub,
                    pattern = ".xlsx",
                    full.names = T)

    lapply_with_error(X = xlsx,FUN = function(x){

       csv=gsub("xlsx", "csv", x)

       # convert(in_file=x,out_file = csv)#xlsx to csv first

       df=read.csv(file = csv,header = T,stringsAsFactors = F)
       id.imaging=unlist(str_extract_all(string = csv,pattern = "(?<=Run\\d_)\\d*(?=.csv)"))#Extract sub id from file name
       path_output=unlist(str_extract_all(string = csv,pattern = "^D.*(?=/Run)"))

       #Check if each files have 38 trials
       if (nrow(df)!=NUM_TRIALS){
         print(paste0(x,"Only have",nrow(df),"trials!"))
         }

       df%>%
         select(
                run_id=ExperimentName,
                Date=SessionDate,
                trial_id=TrialList,
                trial_order=TrialList.Sample,

                fix_onset=FixationCross.OnsetTime,
                fix_duration=FixationCross.OnsetToOnsetTime,

                fraccomp_sti_onset=FracComp.OnsetTime,
                fraccomp_sti_type=TypeDescription,
                fraccomp_sti_value_pair=Pair,
                fraccomp_sti_dis=Dist,
                fraccomp_sti_dis_type=Bin,

                fraccomp_ans=FracComp.CRESP,
                fraccomp_resp_onset=FracComp.RTTime,
                fraccomp_resp_RT=FracComp.RT,
                fraccomp_resp=FracComp.RESP,
                fraccomp_resp_acc=FracComp.ACC,


                congruent=(!!VAR_NAME_CONGRUENCY)
                )%>%
         mutate(
                fraccomp_sti_dis=fraccomp_sti_dis*VAR_ADJUST_DIST_DIRECTION, # Set the left one as baseline: dist = right-left
                sub_id=id.imaging,
                run_num=paste0("r",unlist(str_extract_all(string=run_id,
                                       pattern=PATTERN_RUN_NUM,simplify =F))),
                # Note that the position might be incorrect since "fraccomp_sti_value_pair" itself doesn't encode position info.
                # Need to be adjusted by "dist" (with direction info) . See below
                fraccomp_sti_value_left=
                  unlist(unname(as.frac.numeric(str_extract_all(string=fraccomp_sti_value_pair,
                                             pattern="^.*(?=_)",simplify=T)))),
                fraccomp_sti_value_right=
                  unlist(unname(as.frac.numeric(str_extract_all(string=fraccomp_sti_value_pair,
                                             pattern="(?<=_).*",simplify=T)))),
                Date=gsub(str_extract_all(string = df$SessionDate[1],
                                     pattern = "^\\d.*(?=T)",simplify = T),
                          pattern = "-",replacement="_"))%>%
         # Swap the left and the right if misplaced
         mutate(fraccomp_sti_value_left_adjusted=ifelse((fraccomp_sti_value_right-fraccomp_sti_value_left)*fraccomp_sti_dis>0,
                                               yes = fraccomp_sti_value_left,no = fraccomp_sti_value_right),
                fraccomp_sti_value_right_adjusted=ifelse((fraccomp_sti_value_right-fraccomp_sti_value_left)*fraccomp_sti_dis>0,
                                               yes = fraccomp_sti_value_right,no = fraccomp_sti_value_left))%>%
         # Retain the adjusted ones
         select(-fraccomp_sti_value_right,-fraccomp_sti_value_left)%>%
         rename(fraccomp_sti_value_left=fraccomp_sti_value_left_adjusted,
                fraccomp_sti_value_right=fraccomp_sti_value_right_adjusted)%>%
         write.csv(x = .,
                   file=paste0(path_output,
                               "\\",STRING_SUB_PREFIX,id.imaging,"_",
                               str_extract_all(string = x,pattern = "Run\\d(?=)",simplify = T),
                               "_tidy.csv"))

  })
})
#Add in sti duration info
pblapply(SUBJECT_LIST,function(sub){

  directory_sub=paste0(PATH_BEH_DATA,"\\",sub)
  csv_list=list.files(path = directory_sub,
                  pattern = "tidy.csv",
                  full.names = T)
  lapply_with_error(X = csv_list,FUN = function(x){

      df=read.csv(x,stringsAsFactors = F,header = T)
      df=
        df%>%
          mutate(fraccomp_sti_duration=lead(fix_onset)-fraccomp_sti_onset)
      df[nrow(df),"fraccomp_sti_duration"]=round(mean(df$fraccomp_sti_duration,na.rm = T))
      write.csv(df,file = x) # Overwrite the orignal file
  })
})
#Remove extra undesired columns (whenever needed)
# pblapply(SUBJECT_LIST,function(sub){
#
#   directory_sub=paste0(PATH_BEH_DATA,"\\",sub)
#   csv_list=list.files(path = directory_sub,
#                       pattern = "tidy.csv",
#                       full.names = T)
#   lapply_with_error(X = csv_list,FUN = function(x){
#
#     df=read.csv(x,stringsAsFactors = F,header = T)
#     df=
#       df%>%
#         select(-X.1,-X.2)
#
#     write.csv(df,file = x) # Overwrite the orignal file
#   })
# })

#Rename the csv (whenever needed)
pblapply(SUBJECT_LIST,function(sub){

  directory_sub=paste0(PATH_BEH_DATA,"\\",sub)
  csv_list=list.files(path = directory_sub,
                      pattern = "tidy.csv",
                      full.names = T)
  lapply_with_error(X = csv_list,FUN = function(x){
    old_name=x
    new_name=gsub(x = x,pattern = "df",replacement = "XFC")
    file.rename(old_name,new_name)
  })
})
