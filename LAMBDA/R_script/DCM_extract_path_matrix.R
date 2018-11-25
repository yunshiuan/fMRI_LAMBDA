#===========================================================================================
# Extract path and effect values from DCM.mat results
#===========================================================================================
library("R.matlab")
library("pbapply")
library("dplyr")
library("tidyr")
library("stringr")
# ------------------------------------------------------------------------------------
# (Part1) EXtract and collect path matrix(i.e., A, B, C matrices) from DCM_estimated.mat to a data frame------------------
# ------------------------------------------------------------------------------------
# Constant
# To be adjusted according to need
DCM_dir="D:\\Yun-Shiuan_LAMBDA\\DCM\\Part2_DCM_Specify_and_Estimate_with_reslice_normalized_no_SMA\\full_duration\\resliced_peaks\\Houde_resliced_peak_4mm_C_matrix_only_V1\\B_full_specified"
output_dir="D:\\Yun-Shiuan_LAMBDA\\DCM\\Part3_second_level_LMER\\full_duration\\Houde_resliced_peak_4mm_C_matrix_only_V1"
#Read in the DCM_estimated.mat file
#List all DCM result files
temp=list.files(path = DCM_dir,pattern = "result_estimated_run",full.names = T,recursive = T)
# region_name=c("L_IPS","L_M1","L_SMA","L_V1",
#               "R_IPS","R_SMA","R_V1","R_DLPFC")
region_name=c("L_IPS","L_M1","L_V1",
              "R_IPS","R_V1","R_DLPFC")

list_path_collect_all_id=
  pblapply(X = temp,FUN = function(mat){
  # Read in parameters
  Mat=readMat(mat)
  A=Mat$Ep[[1]]
  B_Near=Mat$Ep[[2]][,,1]
  B_Medium=Mat$Ep[[2]][,,2]
  B_Far=Mat$Ep[[2]][,,3]

  C_Near=Mat$Ep[[3]][,1]
  C_Medium=Mat$Ep[[3]][,2]
  C_Far=Mat$Ep[[3]][,3]

  # Rename (add in the region name information)
  colnames(A)=region_name
  colnames(B_Near)=region_name
  colnames(B_Medium)=region_name
  colnames(B_Far)=region_name

  rownames(A)=region_name
  rownames(B_Near)=region_name
  rownames(B_Medium)=region_name
  rownames(B_Far)=region_name

  names(C_Near)=region_name
  names(C_Medium)=region_name
  names(C_Far)=region_name

  #Reshape the matrices into data frame
  df_A=
    data.frame(A)%>%
    gather(key = "region_from")%>%
    mutate(region_to=rep(region_name,length.out=length(region_name)^2),
           matrix_type="A")

  df_B_Near=
    data.frame(B_Near)%>%
    gather(key = "region_from")%>%
    mutate(region_to=rep(region_name,length.out=length(region_name)^2),
           matrix_type="B_Near")
  df_B_Medium=
    data.frame(B_Medium)%>%
    gather(key = "region_from")%>%
    mutate(region_to=rep(region_name,length.out=length(region_name)^2),
           matrix_type="B_Medium")
  df_B_Far=
    data.frame(B_Far)%>%
    gather(key = "region_from")%>%
    mutate(region_to=rep(region_name,length.out=length(region_name)^2),
           matrix_type="B_Far")

  df_C_Near=
    data.frame(value=C_Near)%>%
    mutate(region_to=rep(region_name,length.out=length(region_name)),
           matrix_type="C_Near",
           region_from=NA)
  df_C_Medium=
    data.frame(value=C_Medium)%>%
    mutate(region_to=rep(region_name,length.out=length(region_name)),
           matrix_type="C_Medium",
           region_from=NA)
  df_C_Far=
    data.frame(value=C_Far)%>%
    mutate(region_to=rep(region_name,length.out=length(region_name)),
           matrix_type="C_Far",
           region_from=NA)

  #Bind sub-dataframes into a big dataframe
  df_path_collect=do.call("rbind",list(df_A,
                                       df_B_Far,df_B_Medium,df_B_Near,
                                       df_C_Far,df_C_Medium,df_C_Near))
  # Add ID and run info
  id=unlist(str_extract_all(string = mat,pattern = "df\\d+(?=/DCM_model_result)"))
  run=unlist(str_extract_all(string = mat,pattern = "(?<=estimated_)run\\d+(?=.mat)"))
  df_path_collect=
    df_path_collect%>%
      mutate(sub_id=toupper(id),
             run=run)
  return(df_path_collect)
  })

# Collapse the list into a big data frame
df_path_collect_all_id=do.call("rbind",list_path_collect_all_id)

# Write the result to csv
write.csv(x = df_path_collect_all_id,
          file = paste0(output_dir,"\\collect_all_path_coefficients.csv"))
