# install.packages("devtools")
#library(devtools)
#install_github("yunshiuan/label4MRI",subdir = "label4MRI") 
library(mni2aal)
library(dplyr)
library(tidyr)
library(stringr)
library(pbapply)

path="D:\\Yun-Shiuan_LAMBDA\\"
result_tables=list.files(pattern = "k15.csv", # Make sure it only include result report csv
                        recursive = T,path = path,
                        full.names = T)
# Make sure it only includeds non-labeled csv
result_tables=grep(x  = result_tables,pattern = "labeled",value = T,invert = T)

error_collection=c()
pblapply(X =result_tables,
         function(csv){
               tryCatch({
                
                df=read.table(csv,header = F,stringsAsFactors = F,
                               sep = ",",col.names=paste0("V",seq_len(15)),fill = T)
                 mni_x=as.numeric(df$V12[3:nrow(df)])
                 mni_y=as.numeric(df$V13[3:nrow(df)])
                 mni_z=as.numeric(df$V14[3:nrow(df)])
                 
                 aal_info=mapply(mni_to_region_name,mni_x,mni_y,mni_z)
                 aal_name=unlist(aal_info[1,])
                 aal_distance=unlist(aal_info[2,])
                 df$V15=c("","aal_name",aal_name)
                 df$V16=c("","aal_distance",round(aal_distance,2))
                 
                 file_name=paste0(str_extract_all(csv,pattern = ".*Result_table(?=/)"),"/",
                                  "labeled_",
                                  str_extract_all(csv,pattern = "(?<=Result_table/).*csv"))
                 
                 write.table(df,file=file_name,
                             sep=",",  col.names=F,row.names = F)
               },
               warning=function(e){
                 error_message=paste0(message(e),"====",csv)
                 print(error_message)
                 error_collection=append(error_collection,
                                         values = paste0(error_message,"====",csv))
                }, 
               error=function(e){
                 error_message=paste0(message(e),"====",csv)
                 print(error_message)
                 error_collection=append(error_collection,
                                         values = paste0(error_message,"====",csv))
               })
         })