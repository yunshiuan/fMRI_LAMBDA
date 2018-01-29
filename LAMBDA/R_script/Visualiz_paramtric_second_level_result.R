#Visualizd Parametric Secone Level Results
#Visulaize all
##(1)Dist. results, including Dist_all,Dist_LL,Dist_FL,Dist_FF
##(2)B0 results, including B0_LL. B0_FL.B0_FF
library(mni2aal)
library(dplyr)
library(tidyr)
library(stringr)
library(pbapply)
library(ggplot2)
library(ggrepel)

#Read in the demographic file (for age and gneder)
df_demographic=read.csv("D:\\Yun-Shiuan_LAMBDA\\demographic\\tidy_demographic.csv",
                        header = T,stringsAsFactors = F)
df_demographic= # muetate the ID column to match the beta files
  df_demographic%>%
    rename(ID=sub_id)%>%
    mutate(ID=as.integer(str_extract_all(ID,pattern="\\w{2}$",simplify = T)))
#Read in the Dist Betas files
beta_file_path="D:\\Yun-Shiuan_LAMBDA\\Second_level_parametric_new\\First_Betas_extract_from_Second_level_peaks"
setwd(beta_file_path)
temp=list.files(path = beta_file_path,
                pattern="^(Dist_|B0_\\w{2}_from)")

pblapply(X = temp,FUN = function(x){
  filename=x # save the filename for further usage(e.g., retrive the mni coordinates)
  #Create the plot title for ggplot below
  k_size=str_extract_all(string = filename,pattern = "(?<=_k_)\\d+(?=_x)",simplify = T)
  k_size=ifelse(length(k_size)==0,yes = "NA",no = k_size)# set to NA if it's subpeak
  
  mni_x=as.numeric(str_extract_all(string = filename,pattern = "(?<=_x).*\\d+(?=_y)",simplify = T))
  mni_y=as.numeric(str_extract_all(string = filename,pattern = "(?<=_y).*\\d+(?=_z)",simplify = T))
  mni_z=as.numeric(str_extract_all(string = filename,pattern = "(?<=_z).*\\d+(?=_aal)",simplify = T))
  
  aal_info=mni_to_region_name(mni_x,mni_y,mni_z)
  aal_name=aal_info$region
  aal_dist=aal_info$distance
  
  first_level_beta=str_extract_all(filename,pattern = "^.*(?=_from)",simplify = T)
  second_level_peak=str_extract_all(filename,pattern = "(?<=from_).*(?=_pk_)",simplify = T)
  peak_num=str_extract_all(filename,pattern = "(?<=_pk_).*(?=_k_)",simplify = T)
  
  plot_title=paste0("1st = ",first_level_beta,
                    "; 2nd = ",second_level_peak,
                    "; peak order = ",peak_num,
                    "; \n k = ",k_size,
                    "; aal Region = ", aal_name,
                    "; dist. = ",aal_dist
                    )
  output_file_name=paste0(str_extract_all(string = filename,
                                   pattern = "^.*(?=_aal)",simplify = T),
                          "_",aal_name,".pdf")
    
  df=read.csv(x,header = F)
  names(df)=c("Beta","ID")
  g=
    df%>%
      left_join(df_demographic%>%select(ID,gender,scan_age_precise,grade),
                by = "ID")%>%
      ggplot(aes(x=scan_age_precise,y=Beta))+
      geom_point(aes(color=gender))+
      geom_smooth(method = "lm")+
      geom_text_repel(aes(label=ID))+
      labs(x="Age",y="Paramater Estimate (a.u.)",
           title=plot_title)
    
  
  ggsave(plot = g,filename = paste0("D:\\Yun-Shiuan_LAMBDA\\Second_level_parametric_new\\Visuaize_First_Betas_extract_from_Second_level_peaks\\",
                                    output_file_name),device = "pdf",
         height = 6,width = 8)
  })
#Scatterplot (Dist Betas~Age)