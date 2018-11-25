#Add the missing column in adults' behavioral data
#E.g., TypeDescription, Pair
library("dplyr")
library("tidyr")
library("stringr")
library("rio")
library("pbapply")
#Declare constants
PATH_BEH_DATA="D:\\Yun-Shiuan_LAMBDA\\Adult\\EprimeData_raw"
STRING_SUB_PREFIX="XFC3"
NUM_TRIALS=38
SUBJECT_LIST=dir(path = PATH_BEH_DATA,
                 pattern=paste0("^",STRING_SUB_PREFIX))
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

    convert(in_file=x,out_file = csv)#xlsx to csv first

    df=read.csv(file = csv,header = T,stringsAsFactors = F)
    id.imaging=unlist(str_extract_all(string = csv,pattern = "(?<=Run\\d_)\\d*(?=.csv)"))#Extract sub id from file name
    path_output=unlist(str_extract_all(string = csv,pattern = "^D.*(?=/Run)"))

    #Check if each files have 38 trials
    if (nrow(df)!=NUM_TRIALS){
      print(paste0(x,"Only have",nrow(df),"trials!"))
    }

    df%>% #for adults(original: dist= right-left)
      mutate(TypeDescription=case_when(Type=="Frac-Frac"& Dist>0 ~ "f - F",
                                       Type=="Frac-Frac"& Dist<0 ~ "F - f",
                                       Type=="Line-Frac"& Dist>0 ~ "f - L",
                                       Type=="Line-Frac"& Dist<0 ~ "l - F",
                                       Type=="Line-Line"& Dist>0 ~ "l - L",
                                       Type=="Line-Line"& Dist<0 ~ "L - l"
                                       ),
             value_left=ifelse(!is.na(N1),yes = paste0(N1,"/",D1),
                                no = gsub(x=unlist(str_extract_all(string=img1,pattern="^\\d+_\\d+(?=_)")),
                                          pattern="_",replacement="/")),
             value_right=ifelse(!is.na(N2),yes = paste0(N2,"/",D2),
                               no = gsub(x=unlist(str_extract_all(string=img2,pattern="^\\d+_\\d+(?=_)")),
                                         pattern="_",replacement="/")),
             Pair=paste0(value_left,"_",value_right)
             )%>%
      select(-value_left,-value_right)%>%
      write.csv(x = .,
                file=csv)

  })
})

