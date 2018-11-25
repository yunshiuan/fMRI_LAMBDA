#Functions for Representational Similarity Analysis(RSA)[Visualize Conceptual model]----------------------------------
#Note: this is sourced by "RSA_conceptual_models_visualization" and "RSA_conceptual_models_output_csv"
library("dplyr")
library("tidyr")
library("stringr")
library("reshape2")
library("ggplot2")
#Constants------------------------------------------------
SCALE_MIN=0
SCALE_MAX=1

# #The color bar for white to dark blue
# COLOR_LOW="#ffffff"
# COLOR_HIGH="#002b80"
#The color bar for bright blue to dark blue
COLOR_LOW="#56B1F7"
COLOR_HIGH="#132B43"
#Basic RSA functions
FILE_RSA_BASIC_FUNCTIONS="D:\\GoogleDrive\\Lambda_code\\R_script\\R_functions\\RSA_basic_functions.R"
# Helper functions=========================================)==============================
#Basic RSA functions
source(FILE_RSA_BASIC_FUNCTIONS)
#Function for preprocessing data frame into RSA required format-----------
preprocess_df_discrete=function(csv){
  processed_df=
    read.csv(csv)%>%
    select(-matches("(^X.\\d$)|(X)"))%>%
    filter(!is.na(trial_id))%>%
    mutate(fraccomp_sti_dis_abs=abs(fraccomp_sti_dis),
           sign=ifelse(fraccomp_sti_dis>=0,1,-1),
           fraccomp_sti_dis_mean=fraccomp_sti_dis_abs_mean*sign)%>%
    mutate(F_left=grepl("^(f|F)",fraccomp_sti_type),
           F_right=grepl("(f|F)$",fraccomp_sti_type),
           L_left=grepl("^(l|L)",fraccomp_sti_type),
           L_right=grepl("(l|L)$",fraccomp_sti_type))%>%
    mutate(F_amount=F_left+F_right,
           L_amount=L_left+L_right)%>%
    mutate(#Continuous distance
           fraccomp_sti_dis_abs_scaled=scale_RSA(fraccomp_sti_dis_abs,SCALE_MIN,SCALE_MAX),
           fraccomp_sti_dis_scaled=scale_RSA(fraccomp_sti_dis,SCALE_MIN,SCALE_MAX),
           #Discrete distance
           fraccomp_sti_dis_discrete_scaled=scale_RSA(fraccomp_sti_dis_mean,SCALE_MIN,SCALE_MAX),
           fraccomp_sti_dis_discrete_abs_scaled=scale_RSA(fraccomp_sti_dis_abs_mean,SCALE_MIN,SCALE_MAX))%>%
    mutate_at(vars(matches("RT_mean|ACC_mean")),funs(scale_RSA(.,SCALE_MIN,SCALE_MAX)))%>%
    #Log-distance
    mutate_at(vars(matches("logv\\d")),funs(scale_RSA(.,SCALE_MIN,SCALE_MAX)))
  return(processed_df)
}

preprocess_df_continuous=function(csv){
  processed_df=
    read.csv(csv)%>%
    filter(!is.na(trial_id))%>%
    mutate(fraccomp_sti_dis_abs=abs(fraccomp_sti_dis),
           fraccomp_sti_dis_discrete_abs=abs(fraccomp_sti_dis_discrete))%>%
    mutate(F_left=grepl("^(f|F)",fraccomp_sti_type),
           F_right=grepl("(f|F)$",fraccomp_sti_type),
           L_left=grepl("^(l|L)",fraccomp_sti_type),
           L_right=grepl("(l|L)$",fraccomp_sti_type))%>%
    mutate(F_amount=F_left+F_right,
           L_amount=L_left+L_right)%>%
    arrange(trial_id_RSA)%>%
    mutate(fraccomp_sti_dis_abs_scaled=scale_RSA(fraccomp_sti_dis_abs,0,2),
           fraccomp_sti_dis_scaled=scale_RSA(fraccomp_sti_dis,0,2),
           fraccomp_sti_dis_discrete_scaled=scale_RSA(fraccomp_sti_dis_discrete,0,2),
           fraccomp_sti_dis_discrete_abs_scaled=scale_RSA(fraccomp_sti_dis_discrete_abs,0,2))
  return(processed_df)
}
# 24 trials version: Only retain  unique trails (to match up with neural RDM)
preprocess_df_discrete_unique_trails_24=function(csv){
  processed_df=
    preprocess_df_discrete(csv)%>%
    arrange(trial_id_RSA_discrete)%>%
    distinct(trial_id_RSA_discrete,.keep_all = T)
}
preprocess_df_continuous_unique_trails_24=function(csv){
  processed_df=
    preprocess_df_continuous(csv)%>%
    arrange(trial_id_RSA)%>%
    distinct(trial_id_RSA_discrete,.keep_all = T)
}
# 36 trials version
preprocess_df_discrete_repeated_trails_36=function(csv){
  processed_df=
    preprocess_df_discrete(csv)%>%
    arrange(trial_id_RSA_discrete_36)%>%
    distinct(trial_id_RSA_discrete_36,.keep_all = T)
}
# 18 trials version
preprocess_df_discrete_repeated_trails_18=function(csv){
  processed_df=
    preprocess_df_discrete(csv)%>%
    arrange(trial_id_RSA_discrete_18)%>%
    distinct(trial_id_RSA_discrete_18,.keep_all = T)
}
#Function to convert trial-by-trial long format data frame into dissimilarity matrix (RDM)----------
df_to_RDM=
  function(df,dist_function,parameters){
    matrix_RDM<<-matrix(0,nrow = nrow(df),ncol=nrow(df))

    lapply(1:(nrow(df)-1),FUN = function(trial1_index){
      trial1_content=df[trial1_index,]

        lapply((trial1_index+1):nrow(df),FUN = function(trial2_index){
          trial2_content=df[trial2_index,]
          matrix_RDM[trial1_index,trial2_index]<<-dist_function(trial1_content,trial2_content,parameters)
        })
    })

    matrix_RDM<<-t(matrix_RDM)+matrix_RDM # Symmetry by the diagonal line

    # To rescale the dissimilarity to 0~2 (necessary before mixing models)
    if(length(table(matrix_RDM))>1){
      matrix_RDM=scale_RSA(matrix_RDM,new_min = SCALE_MIN, new_max = SCALE_MAX)
    #Special case: Don't scale if  all the values are the same,
    }else if(length(table(matrix_RDM))==1){}
    
    diag(matrix_RDM)<-0 #Diagonal elements are identical to itself hence disimilarity=0

    return(matrix_RDM)
  }
#Function to perform rank-transformation----------------------------------------
#Function to generate plot from RDM matrix--------------------------------------
#Serve as a helper function for both RDM_model_plot() and RDM_neural_plot()
RDM_plot=function(RDM_matrix,num_trial=18,plot_title,
                  plot_rank_transform_and_scale =  F,RDM_rank_transform_and_scale=F,
                  show_legend = T,show_axis_text = T){
  # Accomodate different version of number of trials (24,36,18) Defaults to n = 18.
  if(num_trial==24){
    intervals=c(0.5,1.5,2.5,3.5)
    batch_size=num_trial/4
    labels=c("FF", "FL","LF", "LL")
    breaks=intervals*batch_size
    
  }else if(num_trial==36){
    intervals=c(0.5,1.25,1.75,2.5)
    batch_size=num_trial/3
    labels=c("FF", "FL","LF", "LL")
    breaks=intervals*batch_size
    
  }else if (num_trial==18){
    
    # Enable tick position to differ from label position
    ticks=c(6.5,12.5)
    breaks=c(3.5,9.5,15.5,ticks)
    labels=c("FF", "CN", "LL","","")
    
    # intervals=c(0.5,1.5,2.5)
    # batch_size=num_trial/3
    # labels=c("FF", "CN", "LL")
    # breaks=intervals*batch_size
  }
  
  #Convert the RDM matrix to long format data frame
  melt_RDM = melt(RDM_matrix)
  
  #Rank transform and scale to 0 ~ 100 (as percentile)
  if(plot_rank_transform_and_scale){
    SCALE_MAX=100
    melt_RDM=
      melt_RDM%>%
      mutate(value=rank(value,ties.method = "average"),
             value=scale_RSA(value,new_min = SCALE_MIN,new_max = SCALE_MAX-0.000001))
  }
  
  if(RDM_rank_transform_and_scale){
    SCALE_MAX=100
    RDM_matrix=rank(RDM_matrix,ties.method = "average")
    RDM_matrix=scale_RSA(RDM_matrix,new_min = SCALE_MIN,new_max = SCALE_MAX-0.000001)
  }
    
  g=melt_RDM%>%
    ggplot(aes(x=Var1,y=Var2,fill=value))+
    geom_tile()+
    scale_y_reverse(breaks = breaks,
                    labels = labels,
                    expand=c(0,0))+
    scale_x_continuous(breaks = breaks,
                       labels = labels,
                       expand=c(0,0))+
    # scale_fill_gradient2(name="Dissimilarity\n",
    #                     limits=c(SCALE_MIN,SCALE_MAX),#Set the limits of the color bar
    #                     midpoint = 0.5,
    #                     low=scales::muted("blue"),high=scales::muted("red"))+# Set the color of the color bar
    scale_fill_gradient(name="Dissimilarity \n [percentile]",
                        limits=c(SCALE_MIN,SCALE_MAX),# Set the limits of the color bar
                        low=COLOR_LOW,high=COLOR_HIGH)+# Set the color of the color bar
    coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE)+ #Fix the x/y axis as having same units
    labs(title=plot_title,
         x="trials",y="trials",fill="dissimilarity")+
    theme_classic()+
    theme(plot.title = element_text(hjust = 0.5),
          legend.title = element_text(size=13),
          legend.text = element_text(size=10),
          legend.key.height =unit(1.5,"line"),
          axis.text=element_text(size=14),
          axis.ticks = element_line(color = c(rep(NA, length(ticks) + 1),
                                              rep("black", length(ticks))),
                                    size=rel(1.5)),
          axis.line=element_blank())
  
  #Remove legend if not needed
  if(!show_legend){
    g=g+theme(legend.position = "none")
  }
  #Remove axix text if not needed
  if(!show_axis_text){
    g=g+theme(axis.title = element_blank(),
              axis.text = element_blank())
  }
  return(list(plot=g,RDM=RDM_matrix))
}
#Function to generate conceptual RDM plot according to different dissimilarity measures--------
RDM_model_plot=
  function(df,dis_measure_function,plot_title,parameters,
           plot_rank_transform_and_scale=F,RDM_rank_transform_and_scale=F,
           show_legend=T,show_axis_text=T){

    # Check the amount of trials and adjust the axis labels accordingly (get the number from the name of the dataframe)
    num_trial=as.numeric(unlist(str_extract_all(deparse(substitute(df)),pattern = "\\d+$")))
    
    #Compute the conceptual RDM
    RDM_model_matrix = df_to_RDM(df,dis_measure_function,
                          parameters = parameters)
    #Plot the RDM and return the RDM matrix
    RDM_plot(RDM_matrix=RDM_model_matrix,
             num_trial=num_trial,
             plot_title = plot_title,
             plot_rank_transform_and_scale = plot_rank_transform_and_scale,
             RDM_rank_transform_and_scale = RDM_rank_transform_and_scale,
             show_legend = show_legend,
             show_axis_text = show_axis_text)
  }
#Function to visualize RDM from csv sheet------------------------------
RDM_csv_plot=function(csv,num_trial=18,plot_title,
                      plot_rank_transform_and_scale =  F,RDM_rank_transform_and_scale=F,
                      show_legend=T,show_axis_text=T){
  # num_trial: Accomodate different version of number of trials (24,36,18) Defaults to n = 18.
  
  # Read in neural RDM
  RDM_matrix = unname(as.matrix(read.csv(csv,header = F)))
  
  #Plot the RDM and return the RDM matrix
  RDM_plot(RDM_matrix=RDM_model_matrix,
           num_trial=num_trial,
           plot_title = plot_title,
           plot_rank_transform_and_scale = plot_rank_transform_and_scale,
           RDM_rank_transform_and_scale = RDM_rank_transform_and_scale,
           show_legend = show_legend,
           show_axis_text = show_axis_text)
}
# Dissimilarity functions===================================)===============================================
# No Distance===============================================)============
#Model: Null model - Doing nothing-------------------------------
dist_null=
  function(trial1_content,trial2_content,parameters){
    return(SCALE_MAX) # Null model
  }
#Model: FCF(Full cross-format similarity) All-the-same model - Highest level of abstraction-------------------------
dist_all_the_same=
  function(trial1_content,trial2_content,parameters){
    return(SCALE_MIN) # All-the-same model
  }
#Model: Format with location information-------------------------------------
dist_format_with_location_no_distance=
  function(trial1_content,trial2_content,parameters){
    trial1_features=
      trial1_content%>%
      select(F_left,F_right,L_left,L_right)

    trial2_features=
      trial2_content%>%
      select(F_left,F_right,L_left,L_right)
    return(sum(abs(trial1_features-trial2_features))/2)
  }
#Model: GCF(Graded cross-format similarity) Format without location information-------------------------------------
dist_format_without_location_no_distance=
  function(trial1_content,trial2_content,parameters){
    trial1_features=
      trial1_content%>%
      select(F_amount,L_amount)

    trial2_features=
      trial2_content%>%
      select(F_amount,L_amount)
    return(sum(abs(trial1_features-trial2_features))/2)
  }
#Model: NCF(No cross-format similarity)Pair with the same format----------------------
dist_same_format_pair_no_distance=
  function(trial1_content,trial2_content,parameters){
    trial1_features=
      trial1_content%>%
      select(F_amount)

    trial2_features=
      trial2_content%>%
      select(F_amount)
    return(trial1_features!=trial2_features) # True only if the format is with different format (1: diff; 0: same format)
  }
#Hybrid Models Definitions===============================)=================
#Hybrid Distance Definition-------------------------------------------------
dist_hybrid=
  function(dist_func_1,dist_func_2,trial1_content,trial2_content,weight_dist_2){
    w=weight_dist_2# The coeffient of how much the 2nd dist. weighs compare to 1st dist.

    first_diff=dist_func_1(trial1_content,trial2_content) #e.g., format diff
    second_diff=dist_func_2(trial1_content,trial2_content) # e.g., distance diff

    return(as.numeric((1-w)*first_diff+(w*second_diff)))
  }
#Hybrid(partial) Dist. in particular batch Definition--------------------------
#Hybrid the first dist. with the second dist. only in sepcific batch(e.g., LL-LL)
dist_partial_hybrid= #Only hybrid two models within the filter
  function(dist_func_1,dist_func_2,trial1_content,trial2_content,weight_dist_2,hybrid_filter){
    if(hybrid_filter==T){ #hybrid the pixel within the filter
      w=weight_dist_2# The coeffient of how much the 2nd dist. weighs compare to 1st dist.

      hybrid_dist=dist_hybrid(dist_func_1 = dist_func_1,
                              dist_func_2 = dist_func_2,
                              weight_dist_2 = weight_dist_2,
                              trial1_content =  trial1_content,
                              trial2_content =  trial2_content)

      return(hybrid_dist)
    }else{ #Do not hybrid the pixel outside the filter
      unhybrid_dist=dist_func_1(trial1_content =  trial1_content,
                                trial2_content =  trial2_content)
      return(unhybrid_dist)
    }
  }
#Additive Hybrid Defitition-------------------------------------------------
#The main purpose of this is for LL-batch distance-weighted models.
#Since smaller distance should only make conditions in the LL batch more similar,
#larger distance between two LL conditions should not be more dissimilar than LL-CN.
#That is, hybriding LL and format model in the LL batch should only increase the similarities.
dist_addictive_hybrid=
  function(dist_func_1,dist_func_2,trial1_content,trial2_content,weight_dist_2){
    w=weight_dist_2# The coeffient of how much the 2nd dist. weighs compare to 1st dist.

    first_diff=dist_func_1(trial1_content,trial2_content) #e.g., format diff
    second_diff=dist_func_2(trial1_content,trial2_content) # e.g., distance diff
    #The larger the w, the smaller the disimilarity (due to increasing importance of closer distance).
    #Shift the distance RDM by -2 so that the max becomes 0 (won't increase dissimilarity).
    return(as.numeric(first_diff+w*(second_diff-2)))
  }
#Additive Hybrid(partial) Dist. in particular batch Definition--------------------------
#Hybrid the first dist. with the second dist. only in sepcific batch(e.g., LL-LL)
dist_addictive_partial_hybrid= #Only hybrid two models within the filter
  function(dist_func_1,dist_func_2,trial1_content,trial2_content,weight_dist_2,hybrid_filter){
    if(hybrid_filter==T){ #hybrid the pixel within the filter
      w=weight_dist_2# The coeffient of how much the 2nd dist. weighs compare to 1st dist.

      hybrid_dist=dist_addictive_hybrid(dist_func_1 = dist_func_1,
                                        dist_func_2 = dist_func_2,
                                        weight_dist_2 = weight_dist_2,
                                        trial1_content =  trial1_content,
                                        trial2_content =  trial2_content)

      return(hybrid_dist)
    }else{ #Do not hybrid the pixel outside the filter
      unhybrid_dist=dist_func_1(trial1_content =  trial1_content,
                                trial2_content =  trial2_content)
      return(unhybrid_dist)
    }
  }

#Hybrid Difficulty Definition-------------------------------------------------
#This type of hybrid is different fundamentally from the "dist_hybrid"
# in a sense that there are two parameters instead of one
# ,i.e., weight and age group(for difficulty).
diff_hybrid=
  function(dist_func_other,dist_func_difficulty,trial1_content,trial2_content,weight_dist_2,age_group){
    w=weight_dist_2# The coeffient of how much the 2nd dist. weighs compare to 1st dist.
    parameters=list(age_group=age_group)
    first_diff=dist_func_other(trial1_content,trial2_content) #e.g., format diff
    second_diff=dist_func_difficulty(trial1_content,trial2_content,parameters) # e.g., difficulty diff
    
    return(as.numeric((1-w)*first_diff+(w*second_diff)))
  }

#Hybrid Log Magnitude Definition-------------------------------------------------
#This type of hybrid is different fundamentally from the "dist_hybrid"
# in a sense that there are two parameters instead of one
# ,i.e., weight and log version.
dist_hybrid_log=
  function(dist_func_1,dist_func_2,trial1_content,trial2_content,parameters){
    w=parameters$w2# The coeffient of how much the 2nd dist. weighs compare to 1st dist.
    first_diff=dist_func_1(trial1_content,trial2_content) #e.g., format diff
    #The log version information is stored in the parameters
    second_diff=dist_func_2(trial1_content,trial2_content,parameters=parameters) # e.g., difficulty diff
    
    return(as.numeric((1-w)*first_diff+(w*second_diff)))
  }
#Hybrid(partial) Log Dist. in particular batch Definition--------------------------
#Hybrid the first format dist. with the second log magnitude dist. only in sepcific batch(e.g., LL-LL)
dist_partial_hybrid_log= #Only hybrid two models within the filter
  function(dist_func_1,dist_func_2,trial1_content,trial2_content,parameters,hybrid_filter){
    if(hybrid_filter==T){ #hybrid the pixel within the filter
      hybrid_dist=dist_hybrid_log(dist_func_1 = dist_func_1,
                                  dist_func_2 = dist_func_2,
                                  trial1_content = trial1_content,
                                  trial2_content = trial2_content,
                                  parameters = parameters)
      return(hybrid_dist)
    }else{ #Do not hybrid the pixel outside the filter
      unhybrid_dist=dist_func_1(trial1_content =  trial1_content,
                                trial2_content =  trial2_content)
      return(unhybrid_dist)
    }
  }


# Continuous distance========================================)=======
#Model: absolute distance information (no format or position)-------------------------------------
dist_only_abs_dist=
  function(trial1_content,trial2_content,parameters){
    trial1_features=
      trial1_content%>%
      select(fraccomp_sti_dis_abs_scaled)

    trial2_features=
      trial2_content%>%
      select(fraccomp_sti_dis_abs_scaled)
    dist_diff=abs(trial1_features["fraccomp_sti_dis_abs_scaled"]-trial2_features["fraccomp_sti_dis_abs_scaled"])
    return(as.numeric(dist_diff))
  }
#Model: signed distance information (no format nor position)-------------------------------------
dist_only_signed_dist=
  function(trial1_content,trial2_content,parameters){
    trial1_features=
      trial1_content%>%
      select(fraccomp_sti_dis_scaled)

    trial2_features=
      trial2_content%>%
      select(fraccomp_sti_dis_scaled)
    dist_diff=abs(trial1_features["fraccomp_sti_dis_scaled"]-trial2_features["fraccomp_sti_dis_scaled"])
    return(as.numeric(dist_diff))
  }

# Hybrid models weighted by continuous distance--------------------------------------
#Model: Format with location information(with absolute distance information)-------------------------------------
#w is The relative weight of distance difference compared to format difference
dist_format_with_location_abs_dist=
  function(trial1_content,trial2_content,parameters){
    hybrid_dist=dist_hybrid(dist_func_1 = dist_format_with_location_no_distance,
                            dist_func_2 = dist_only_abs_dist,
                            weight_dist_2 = parameters,
                            trial1_content =  trial1_content,
                            trial2_content =  trial2_content)
    return(hybrid_dist)
  }
#Model: GCF: Format without location information(with absolute distance information)-------------------------------------
#w is The relative weight of distance difference compared to format difference
dist_format_without_location_abs_dist=
  function(trial1_content,trial2_content,parameters){
    hybrid_dist=dist_hybrid(dist_func_1 = dist_format_without_location_no_distance,
                            dist_func_2 = dist_only_abs_dist,
                            weight_dist_2 = parameters,
                            trial1_content =  trial1_content,
                            trial2_content =  trial2_content)
    return(hybrid_dist)
  }
#Model: Format with location information(with signed distance information)-------------------------------------
#w is The relative weight of distance difference compared to format difference
dist_format_with_location_signed_dist=
  function(trial1_content,trial2_content,parameters){
    hybrid_dist=dist_hybrid(dist_func_1 = dist_format_with_location_no_distance,
                            dist_func_2 = dist_only_signed_dist,
                            weight_dist_2 = parameters,
                            trial1_content =  trial1_content,
                            trial2_content =  trial2_content)
    return(hybrid_dist)
  }
#Model: GCF: Format without location information(with signed distance information)-------------------------------------
#w is The relative weight of distance difference compared to format difference
dist_format_without_location_signed_dist=
  function(trial1_content,trial2_content,parameters){
    hybrid_dist=dist_hybrid(dist_func_1 = dist_format_without_location_no_distance,
                            dist_func_2 = dist_only_signed_dist,
                            weight_dist_2 = parameters,
                            trial1_content =  trial1_content,
                            trial2_content =  trial2_content)
    return(hybrid_dist)
  }
# Discrete distance========================================)=======
#Model: absolute distance information (no format or position)-------------------------------------
dist_only_abs_dist_discrete=
  function(trial1_content,trial2_content,parameters){
    trial1_features=
      trial1_content%>%
      select(fraccomp_sti_dis_discrete_abs_scaled)

    trial2_features=
      trial2_content%>%
      select(fraccomp_sti_dis_discrete_abs_scaled)
    dist_diff=abs(trial1_features["fraccomp_sti_dis_discrete_abs_scaled"]-trial2_features["fraccomp_sti_dis_discrete_abs_scaled"])
    return(as.numeric(dist_diff))
  }
#Model: signed distance information (no format or position)-------------------------------------
dist_only_signed_dist_discrete=
  function(trial1_content,trial2_content,parameters){
    trial1_features=
      trial1_content%>%
      select(fraccomp_sti_dis_discrete_scaled)

    trial2_features=
      trial2_content%>%
      select(fraccomp_sti_dis_discrete_scaled)
    dist_diff=abs(trial1_features["fraccomp_sti_dis_discrete_scaled"]-trial2_features["fraccomp_sti_dis_discrete_scaled"])
    return(as.numeric(dist_diff))
  }

# Hybrid models weighted by discrete distance---------------------)----------------------
#Model: Format with location information(with absolute distance information)-------------------------------------
#w is The relative weight of distance difference compared to format difference
dist_format_with_location_abs_dist_discrete=
  function(trial1_content,trial2_content,parameters){
    hybrid_dist=dist_hybrid(dist_func_1 = dist_format_with_location_no_distance,
                            dist_func_2 = dist_only_abs_dist_discrete,
                            weight_dist_2 = parameters,
                            trial1_content =  trial1_content,
                            trial2_content =  trial2_content)
    return(hybrid_dist)
  }
#Model: GCF: Format without location information(with absolute distance information)-------------------------------------
#w is The relative weight of distance difference compared to format difference
dist_format_without_location_abs_dist_discrete=
  function(trial1_content,trial2_content,parameters){
    hybrid_dist=dist_hybrid(dist_func_1 = dist_format_without_location_no_distance,
                            dist_func_2 = dist_only_abs_dist_discrete,
                            weight_dist_2 = parameters,
                            trial1_content =  trial1_content,
                            trial2_content =  trial2_content)
    return(hybrid_dist)
  }
#Model: Format with location information(with signed distance information)-------------------------------------
#w is The relative weight of distance difference compared to format difference
dist_format_with_location_signed_dist_discrete=
  function(trial1_content,trial2_content,parameters){
    hybrid_dist=dist_hybrid(dist_func_1 = dist_format_with_location_no_distance,
                            dist_func_2 = dist_only_signed_dist_discrete,
                            weight_dist_2 = parameters,
                            trial1_content =  trial1_content,
                            trial2_content =  trial2_content)
    return(hybrid_dist)
  }
#Model: GCF: Format without location information(with signed distance information)-------------------------------------
#w is The relative weight of distance difference compared to format difference
dist_format_without_location_signed_dist_discrete=
  function(trial1_content,trial2_content,parameters){
    hybrid_dist=dist_hybrid(dist_func_1 = dist_format_without_location_no_distance,
                            dist_func_2 = dist_only_signed_dist_discrete,
                            weight_dist_2 = parameters,
                            trial1_content =  trial1_content,
                            trial2_content =  trial2_content)
    return(hybrid_dist)
  }
#Model: NCF: Pair with the same format (with absolute distance information)-------------------------------------
#w is The relative weight of distance difference compared to format difference
dist_same_format_pair_abs_distance=
  function(trial1_content,trial2_content,parameters){
    hybrid_dist=dist_hybrid(dist_func_1 = dist_same_format_pair_no_distance,
                            dist_func_2 = dist_only_abs_dist_discrete,
                            weight_dist_2 = parameters,
                            trial1_content =  trial1_content,
                            trial2_content =  trial2_content)
    return(hybrid_dist)
  }
#Model: NCF: Pair with the same format (with signed distance information)-------------------------------------
#w is The relative weight of distance difference compared to format difference
dist_same_format_pair_signed_distance=
  function(trial1_content,trial2_content,parameters){
    hybrid_dist=dist_hybrid(dist_func_1 = dist_same_format_pair_no_distance,
                            dist_func_2 = dist_only_signed_dist_discrete,
                            weight_dist_2 = parameters,
                            trial1_content =  trial1_content,
                            trial2_content =  trial2_content)
    return(hybrid_dist)
  }

# Partially Hybrid distance-weighted models-----------------------)--------------------------------------------
# Partially Hybrid within the LL-LL batch----------------------------------
dist_partial_hybrid_in_LL=
  function(dist_func_1,dist_func_2,trial1_content,trial2_content,weight_dist_2,hybrid_filter){
    boolean_both_LL=prod(as.numeric(trial1_content[,c("L_left","L_right")]))*
        prod(as.numeric(trial2_content[,c("L_left","L_right")]))

    hybrid_dist=dist_partial_hybrid(dist_func_1 = dist_func_1,
                                              dist_func_2 = dist_func_2,
                                              trial1_content = trial1_content,
                                              trial2_content = trial2_content,
                                              weight_dist_2 =  weight_dist_2,
                                              hybrid_filter = boolean_both_LL)
    return(hybrid_dist)
  }
#Model: GCF: Format without location information with LL-absolute distance information---------------------------
dist_format_without_location_LL_abs_dist_discrete=
  function(trial1_content,trial2_content,parameters){
    hybrid_dist=dist_partial_hybrid_in_LL(dist_func_1 = dist_format_without_location_no_distance,
                                                    dist_func_2 = dist_only_abs_dist_discrete,
                                                    weight_dist_2 = parameters,
                                                    trial1_content =  trial1_content,
                                                    trial2_content =  trial2_content)
    return(hybrid_dist)
  }
#Model: GCF: Format without location information with LL-signed distance information---------------------------
dist_format_without_location_LL_signed_dist_discrete=
  function(trial1_content,trial2_content,parameters){
    hybrid_dist=dist_partial_hybrid_in_LL(dist_func_1 = dist_format_without_location_no_distance,
                                                    dist_func_2 = dist_only_signed_dist_discrete,
                                                    weight_dist_2 = parameters,
                                                    trial1_content =  trial1_content,
                                                    trial2_content =  trial2_content)
    return(hybrid_dist)
  }
#Model: NCF: Pair with the same format with LL-absolute distance information---------------------------
dist_same_format_pair_LL_abs_dist_discrete=
  function(trial1_content,trial2_content,parameters){
    hybrid_dist=dist_partial_hybrid_in_LL(dist_func_1 = dist_same_format_pair_no_distance,
                                                    dist_func_2 = dist_only_abs_dist_discrete,
                                                    weight_dist_2 = parameters,
                                                    trial1_content =  trial1_content,
                                                    trial2_content =  trial2_content)
    return(hybrid_dist)
  }
#Model: NCF: Pair with the same format with LL-signed distance information---------------------------
dist_same_format_pair_LL_signed_dist_discrete=
  function(trial1_content,trial2_content,parameters){
    hybrid_dist=dist_partial_hybrid_in_LL(dist_func_1 = dist_same_format_pair_no_distance,
                                                    dist_func_2 = dist_only_signed_dist_discrete,
                                                    weight_dist_2 = parameters,
                                                    trial1_content =  trial1_content,
                                                    trial2_content =  trial2_content)
    return(hybrid_dist)
  }

# Difficulty====================================================================)=======
#Model: Difficulty - Accuracy--------------------------------------
dist_diff_acc=
  function(trial1_content,trial2_content,parameters){
    if(!"age_group" %in% names(parameters)){
      stop("Age group must be specified in the parameters as parameters$age_group.")
    }
    age_group=parameters$age_group #all/2_grade/5_grade/adult
    ACC_var=case_when(age_group == "all" ~ "ACC_mean_all",
                      age_group == "2_grade" ~ "ACC_mean_2_grade",
                      age_group == "5_grade" ~ "ACC_mean_5_grade",
                      age_group == "adult" ~ "ACC_mean_adult")
    if(is.na(ACC_var)){
      stop("Age group should be one of the option: all, 2_grade, 5_grade, or adult")
    }
    
    trial1_features=
      trial1_content%>%
      select_(ACC_var)
    
    trial2_features=
      trial2_content%>%
      select_(ACC_var)
    dist_diff=abs(trial1_features[ACC_var]-trial2_features[ACC_var])
    return(as.numeric(dist_diff))
  }
#Model: Difficulty - RT--------------------------------------------
dist_diff_RT=
  function(trial1_content,trial2_content,parameters){
    if(!"age_group" %in% names(parameters)){
      stop("Age group must be specified in the parameters as parameters$age_group.")
    }
    age_group=parameters$age_group #all/2_grade/5_grade/adult
    RT_var=case_when(age_group == "all" ~ "RT_mean_all",
                     age_group == "2_grade" ~ "RT_mean_2_grade",
                     age_group == "5_grade" ~ "RT_mean_5_grade",
                     age_group == "adult" ~ "RT_mean_adult")
    if(is.na(RT_var)){
      stop("Age group should be one of the option: all, 2_grade, 5_grade, or adult")
    }
    trial1_features=
      trial1_content%>%
      select_(RT_var)
    
    trial2_features=
      trial2_content%>%
      select_(RT_var)
    dist_diff=abs(trial1_features[RT_var]-trial2_features[RT_var])
    return(as.numeric(dist_diff))
  }
# Hybrid models weighted by difficulty(ACC)-----------------------)---------------
#Model: GCF: Format without location information(with difficulty information)-------------------------------------
#w is The relative weight of distance difference compared to format difference
dist_format_without_location_diff_acc=
  function(trial1_content,trial2_content,parameters){
    hybrid_dist=diff_hybrid(dist_func_other = dist_format_without_location_no_distance,
                            dist_func_difficulty= dist_diff_acc,
                            weight_dist_2 = parameters$w2,
                            age_group = parameters$age_group,
                            trial1_content =  trial1_content,
                            trial2_content =  trial2_content)
    return(hybrid_dist)
  }
#Model: NCF: Pair with the same format (with difficulty information)-------------------------------------
#w is The relative weight of distance difference compared to format difference
dist_same_format_pair_diff_acc=
  function(trial1_content,trial2_content,parameters){
    hybrid_dist=diff_hybrid(dist_func_other = dist_same_format_pair_no_distance,
                            dist_func_difficulty= dist_diff_acc,
                            weight_dist_2 = parameters$w2,
                            age_group = parameters$age_group,
                            trial1_content =  trial1_content,
                            trial2_content =  trial2_content)
    return(hybrid_dist)
  }

# #Not used=================================)========================================:
# Log-transformed discrete distance======================================================)=======
#Model: absolute log-transformed distance information (no format or position)-------------------------------------
dist_only_abs_dist_discrete_log=
  function(trial1_content,trial2_content,parameters){
    #Decide the version of log-transformation 
    log_version=parameters$log_version
    if(log_version=="v1"){
      vars="fraccomp_sti_dis_abs_logv1_mean"
    }else if(log_version=="v2"){
      vars="fraccomp_sti_dis_abs_logv2_mean"
    }else if(log_version=="v3"){
      vars="fraccomp_sti_dis_abs_logv3_mean"
    }else{
      stop("Version not supported.")
    }
    trial1_features=
      trial1_content%>%
      select(!!vars)

    trial2_features=
      trial2_content%>%
      select(!!vars)
    dist_diff=abs(trial1_features[vars]-trial2_features[vars])
    return(as.numeric(dist_diff))
  }
#Model: signed log-transformed distance information (no format or position)-------------------------------------
dist_only_signed_dist_discrete_log=
  function(trial1_content,trial2_content,parameters){
    log_version=parameters$log_version
    if(log_version=="v1"){
      vars="fraccomp_sti_dis_signed_logv1_mean"
    }else if(log_version=="v2"){
      vars="fraccomp_sti_dis_signed_logv2_mean"
    }else{
      stop("Version not supported.")
    }
    trial1_features=
      trial1_content%>%
      select(!!vars)
    
    trial2_features=
      trial2_content%>%
      select(!!vars)
    dist_diff=abs(trial1_features[vars]-trial2_features[vars])
    return(as.numeric(dist_diff))
  }
# Hybrid models weighted by log-transformed discrete distance---------------------)----------------------
#Model: GCF: Format without location information(with absolute log-transformed distance information)-------------------------------------
#w is The relative weight of distance difference compared to format difference
dist_format_without_location_abs_dist_discrete_log=
  function(trial1_content,trial2_content,parameters){
    hybrid_dist=dist_hybrid_log(dist_func_1 = dist_format_without_location_no_distance,
                                dist_func_2 = dist_only_abs_dist_discrete_log,
                                parameters = parameters,
                                trial1_content =  trial1_content,
                                trial2_content =  trial2_content)
    return(hybrid_dist)
  }
#Model: GCF: Format without location information(with signed log-transformed distance information)-------------------------------------
#w is The relative weight of distance difference compared to format difference
dist_format_without_location_signed_dist_discrete_log=
  function(trial1_content,trial2_content,parameters){
    hybrid_dist=dist_hybrid_log(dist_func_1 = dist_format_without_location_no_distance,
                                dist_func_2 = dist_only_signed_dist_discrete_log,
                                parameters = parameters,
                                trial1_content =  trial1_content,
                                trial2_content =  trial2_content)
    return(hybrid_dist)
  }
#Model: NCF: Pair with the same format (with absolute log-transformed distance information)-------------------------------------
#w is The relative weight of distance difference compared to format difference
dist_same_format_pair_abs_distance_log=
  function(trial1_content,trial2_content,parameters){
    hybrid_dist=dist_hybrid_log(dist_func_1 = dist_same_format_pair_no_distance,
                                dist_func_2 = dist_only_abs_dist_discrete_log,
                                parameters = parameters,
                                trial1_content =  trial1_content,
                                trial2_content =  trial2_content)
    return(hybrid_dist)
  }
#Model: NCF: Pair with the same format (with absolute log-transformed distance information)-------------------------------------
#w is The relative weight of distance difference compared to format difference
dist_same_format_pair_signed_distance_log=
  function(trial1_content,trial2_content,parameters){
    hybrid_dist=dist_hybrid_log(dist_func_1 = dist_same_format_pair_no_distance,
                                dist_func_2 = dist_only_signed_dist_discrete_log,
                                parameters = parameters,
                                trial1_content =  trial1_content,
                                trial2_content =  trial2_content)
    return(hybrid_dist)
  }
# Partially Hybrid log-transformed distance-weighted models-----------------------)--------------------------------------------
# Partially Hybrid by log magnitude within the LL-LL batch----------------------------------
dist_partial_hybrid_in_LL_log=
  function(dist_func_1,dist_func_2,trial1_content,trial2_content,weight_dist_2,parameters = parameters){
    boolean_both_LL=prod(as.numeric(trial1_content[,c("L_left","L_right")]))*
      prod(as.numeric(trial2_content[,c("L_left","L_right")]))
    
    hybrid_dist=dist_partial_hybrid_log(dist_func_1 = dist_func_1,
                                    dist_func_2 = dist_func_2,
                                    trial1_content = trial1_content,
                                    trial2_content = trial2_content,
                                    parameters = parameters,
                                    hybrid_filter = boolean_both_LL)
    return(hybrid_dist)
  }

#Model: GCF: Format without location information with LL-log-transformed absolute distance information---------------------------
dist_format_without_location_LL_abs_dist_discrete_log=
  function(trial1_content,trial2_content,parameters){
    hybrid_dist=dist_partial_hybrid_in_LL_log(dist_func_1 = dist_format_without_location_no_distance,
                                              dist_func_2 = dist_only_abs_dist_discrete_log,
                                              trial1_content =  trial1_content,
                                              trial2_content =  trial2_content,
                                              parameters = parameters)
    return(hybrid_dist)
  }
#Model: GCF: Format without location information with LL-log-transformed signed distance information---------------------------
dist_format_without_location_LL_signed_dist_discrete_log=
  function(trial1_content,trial2_content,parameters){
    hybrid_dist=dist_partial_hybrid_in_LL_log(dist_func_1 = dist_format_without_location_no_distance,
                                          dist_func_2 = dist_only_signed_dist_discrete_log,
                                          trial1_content =  trial1_content,
                                          trial2_content =  trial2_content,
                                          parameters = parameters)
    return(hybrid_dist)
  }
#Model: NCF: Pair with the same format with LL-log-transformed absolute distance information---------------------------
dist_same_format_pair_LL_abs_dist_discrete_log=
  function(trial1_content,trial2_content,parameters){
    hybrid_dist=dist_partial_hybrid_in_LL_log(dist_func_1 = dist_same_format_pair_no_distance,
                                          dist_func_2 = dist_only_abs_dist_discrete_log,
                                          trial1_content =  trial1_content,
                                          trial2_content =  trial2_content,
                                          parameters = parameters)
    return(hybrid_dist)
  }
#Model: NCF: Pair with the same format with LL-log-transformed signed distance information---------------------------
dist_same_format_pair_LL_signed_dist_discrete_log=
  function(trial1_content,trial2_content,parameters){
    hybrid_dist=dist_partial_hybrid_in_LL_log(dist_func_1 = dist_same_format_pair_no_distance,
                                          dist_func_2 = dist_only_signed_dist_discrete_log,
                                          trial1_content =  trial1_content,
                                          trial2_content =  trial2_content,
                                          parameters = parameters)
    return(hybrid_dist)
  }

# # (Not Used)Addictive Partially Hybrid distance-weighted models-----------------------)--------------------------------------------
# # Can't resolve the diagonal cells value problems.
# # That is, if the LL batch is addictively more similar, the unscaled dissimilarity value will be smaller than 0.
# # The negative values are smaller than diagonal cells in other batches after scaling. This does not make sense.
# # If one forces the diagonal cells of other batches to become 0 after scaling,
# # those cells will have smaller dissimalarity values than off-diagonal cells in the same batches,
# # which again, does not make sense.
# # Addictive Partially Hybrid within the LL-LL batch----------------------------------
# dist_addictive_partial_hybrid_in_LL=
#   function(dist_func_1,dist_func_2,trial1_content,trial2_content,weight_dist_2,hybrid_filter){
#     boolean_both_LL=prod(as.numeric(trial1_content[,c("L_left","L_right")]))*
#       prod(as.numeric(trial2_content[,c("L_left","L_right")]))
#
#     hybrid_dist=dist_addictive_partial_hybrid(dist_func_1 = dist_func_1,
#                                               dist_func_2 = dist_func_2,
#                                               trial1_content = trial1_content,
#                                               trial2_content = trial2_content,
#                                               weight_dist_2 =  weight_dist_2,
#                                               hybrid_filter = boolean_both_LL)
#     return(hybrid_dist)
#   }
# #Model: GCF: Format without location information with LL-absolute distance information---------------------------
# dist_format_without_location_LL_abs_dist_discrete=
#   function(trial1_content,trial2_content,parameters){
#     hybrid_dist=dist_addictive_partial_hybrid_in_LL(dist_func_1 = dist_format_without_location_no_distance,
#                                                     dist_func_2 = dist_only_abs_dist_discrete,
#                                                     weight_dist_2 = parameters,
#                                                     trial1_content =  trial1_content,
#                                                     trial2_content =  trial2_content)
#     return(hybrid_dist)
#   }
# #Model: GCF: Format without location information with LL-signed distance information---------------------------
# dist_format_without_location_LL_signed_dist_discrete=
#   function(trial1_content,trial2_content,parameters){
#     hybrid_dist=dist_addictive_partial_hybrid_in_LL(dist_func_1 = dist_format_without_location_no_distance,
#                                                     dist_func_2 = dist_only_signed_dist_discrete,
#                                                     weight_dist_2 = parameters,
#                                                     trial1_content =  trial1_content,
#                                                     trial2_content =  trial2_content)
#     return(hybrid_dist)
#   }
# #Model: NCF: Pair with the same format with LL-absolute distance information---------------------------
# dist_same_format_pair_LL_abs_dist_discrete=
#   function(trial1_content,trial2_content,parameters){
#     hybrid_dist=dist_addictive_partial_hybrid_in_LL(dist_func_1 = dist_same_format_pair_no_distance,
#                                                     dist_func_2 = dist_only_abs_dist_discrete,
#                                                     weight_dist_2 = parameters,
#                                                     trial1_content =  trial1_content,
#                                                     trial2_content =  trial2_content)
#     return(hybrid_dist)
#   }
# #Model: NCF: Pair with the same format with LL-signed distance information---------------------------
# dist_same_format_pair_LL_signed_dist_discrete=
#   function(trial1_content,trial2_content,parameters){
#     hybrid_dist=dist_addictive_partial_hybrid_in_LL(dist_func_1 = dist_same_format_pair_no_distance,
#                                                     dist_func_2 = dist_only_signed_dist_discrete,
#                                                     weight_dist_2 = parameters,
#                                                     trial1_content =  trial1_content,
#                                                     trial2_content =  trial2_content)
#     return(hybrid_dist)
#   }
