#Create the demographic file for the adults
#Include Grade as "Adult" and sub_id
PATH_ADULT_ID="D:\\Yun-Shiuan_LAMBDA\\Adult\\preprocessed_data_with_reslice"
FILE_OUTPUT="D:\\Yun-Shiuan_LAMBDA\\Adult\\demographic\\tidy_demographic.csv"
id_list=list.dirs(PATH_ADULT_ID,full.names = F)
id_list=grep(x = id_list,pattern = "^XFC",value = T)
df=data.frame(sub_id=id_list,
              grade="Adult")
write.csv(df,FILE_OUTPUT)
