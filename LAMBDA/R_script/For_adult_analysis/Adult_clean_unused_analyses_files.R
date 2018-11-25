#Clean unnecessary files in the adult brain data folder (only keep the raw dicom files)
library(dplyr)
library(tidyr)

PATH_RAW_DATA="D:\\Yun-Shiuan_LAMBDA\\Adult\\raw_data"
PATTERN_UNUSED="analyses"#These are data for old analyses
#Clean those unneeded data
dirs_unused=list.dirs(PATH_RAW_DATA,full.names = T,recursive = T)
dirs_unused=grep(dirs_unused,pattern = PATTERN_UNUSED,value = T)

unlink(dirs_unused,recursive = T)
