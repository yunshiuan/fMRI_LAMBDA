# In order to only upload necessary data onto the server,--------------------------------------------------------------
# first check data size on the local PC and see how much I could reduce the data size.
#Helper function
get_file_size_MB=function(file_name){
  size=file.info(file_name)$size/(1024*1024)
  return(size)
}
#Constant
PATH_ROOT="D:\\Yun-Shiuan_LAMBDA"

all_files=list.files(path = PATH_ROOT,recursive = T)
#system('powershell -noprofile -command "ls -r|measure -s Length"')
