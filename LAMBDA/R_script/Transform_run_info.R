
#(PART1) Transform run info. from Eprime raw xlxs to tidy format for MATLAB to read

setwd("C:\\Users\ychuang26\\GoogleDrive\\Lambda_code\\R_script")
library("dplyr")
library("tidyr")
library("stringr")
library("rio")
library("pbapply")
      path="D:\\Yun-Shiuan_LAMBDA\\EprimeData_raw"
      subject_list=dir(path = path,
                       pattern="^DevFrac")
      as.frac.numeric=function(x) 
                      {sapply(x,function(x)eval(parse(text=x)))}#"1/2" to 0.5; works for vector
      
      lapply_with_error <- function(X,FUN,...){    
        lapply(X, function(x, ...) tryCatch(FUN(x, ...),
                                            error=function(e) 
                                              paste0(x," has error:",e)))
      }#define lapply with try-catch

#Transform csv xlsx to csv and then output my desired format of run_info.csv file
pblapply(subject_list,function(sub){
  
    directory_sub=paste0(path,"\\",sub)
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
       # if (nrow(df)!=38){
       #   print(paste0(x,nrow(df)))}
       
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

                congruent=ComponentCongruency)%>%
         mutate(sub_id=id.imaging,
                run_num=str_extract_all(string=run_id,
                                       pattern="^r\\d(?=_xfc)",simplify =T),
                fraccomp_sti_value_left=
                  unlist(unname(as.frac.numeric(str_extract_all(string=fraccomp_sti_value_pair,
                                             pattern="^.*(?=_)",simplify=T)))),
                fraccomp_sti_value_right=
                  unlist(unname(as.frac.numeric(str_extract_all(string=fraccomp_sti_value_pair,
                                             pattern="(?<=_).*",simplify=T)))),
                Date=gsub(str_extract_all(string = df$SessionDate[1],
                                     pattern = "^\\d.*(?=T)",simplify = T),
                          pattern = "-",replacement="_")
                )%>%
         write.csv(x = .,
                   file=paste0(path_output,
                               "\\df",id.imaging,"_",
                               str_extract_all(string = x,pattern = "Run\\d(?=)",simplify = T),
                               "_tidy.csv"))
         
  })
})