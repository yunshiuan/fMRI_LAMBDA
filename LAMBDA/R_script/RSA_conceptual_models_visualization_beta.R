#Representational Similarity Analysis(RSA)[Visualize Conceptual model]----------------------------------
#Beta version:
#Try to generalize the df_to_RDM()
library("dplyr")
library("tidyr")
library("stringr")
library("ggplot2")
library("reshape2")
# save.image("D:\\GoogleDrive\\Lambda_code\\R_script\\RData\\RSA.RData")
# load(file = "/Users/vimchiz/googledrive/Lambda_code/R_script/RData/RSA_model_visualization.RData")
#(Step1): Model Specification---------------------------------------------------------------------------
#Helper functions and constants----------------
AMOUNT_TRIALS_TOTAL=36
AMOUNT_FORMAT=3
AMOUNT_TRIALS_EACH=AMOUNT_TRIALS_TOTAL/AMOUNT_FORMAT
RDM_MIN=0
RDM_MAX=2
#Function to convert trial-by-trial long format data frame into dissimilarity matrix (RDM)----------
df_to_RDM=
  function(df,dist_function,parameters){
    matrix_RDM<<-matrix(0,nrow = nrow(df),ncol=nrow(df))

    lapply(1:nrow(df),FUN = function(trial1_index){
      trial1_content=df[trial1_index,]

      if((trial1_index+1)<nrow(df)){

        lapply((trial1_index+1):nrow(df),FUN = function(trial2_index){
          trial2_content=df[trial2_index,]
          matrix_RDM[trial1_index,trial2_index]<<-dist_function(trial1_content,trial2_content,parameters)
        })
      }
    })

    matrix_RDM<<-t(matrix_RDM)+matrix_RDM # Symmetry by the diagonal line
    diag(matrix_RDM)<-0 #Diagonal elements are identical to itself hence disimilarity=0

    # To rescale the dissimilarity to 0~2 (necessary when distance is intorduced)
    return(matrix_RDM,new_min = 0,new_max = 2)
  }

# Unused attempt to generalize the df_to_RDM
dist_functions=c(dist_format_with_location_no_distance,dist_only_abs_dist)
df_to_RDM=
  function(df,dist_functions,weights=1,parameters){
    # Function to convert trial-by-trial long format data frame into dissimilarity matrix (RDM)
    # This support multiple dissimilarity functions to be linearly combined (e.g., format + distance).
    # Args:
    #   dist_functions: A vector of dissimilarity functions which could be linearly mixed
    #   weights: A vector determining the weights of each dissimilarity function. The vector
    #            should be summed up one. If user not provides the weights, this will by default set
    #            to c(1,0,0...), which makes the function only relies on the first disimilarity function.
    # Return:
    #   The representational dissimilarity matrix (RDM) derived from the dissimilarity functions.

    matrix_RDM<<-matrix(0,nrow = nrow(df),ncol=nrow(df))

    lapply(1:nrow(df),FUN = function(trial1_index){
      trial1_content=df[trial1_index,]

      if((trial1_index+1)<nrow(df)){

        lapply((trial1_index+1):nrow(df),FUN = function(trial2_index){
          trial2_content=df[trial2_index,]

          #Derive the weighted sum based on dist_functions and weights
          dist_weighted_sum=0
          for(dist_fun_index in 1:length(dist_functions)){
            dist_fun=dist_functions[[dist_fun_index]]
            weight=weights[dist_fun_index]
            dist_weighted_sum=dist_weighted_sum + weight*dist_fun(trial1_content,trial2_content,parameters)
          }
          matrix_RDM[trial1_index,trial2_index]<<-dist_function(trial1_content,trial2_content,parameters)
        })
      }
    })

    matrix_RDM<<-t(matrix_RDM)+matrix_RDM # Symmetry by the diagonal line
    diag(matrix_RDM)<-0 #Diagonal elements are identical to itself hence disimilarity=0

    # To rescale the dissimilarity to 0~2 (necessary when distance is intorduced)
    return(scale_RDM(matrix_RDM,new_min = 0,new_max = 2))
  }
#Function to generate plot according to different dissimilarity measures--------
RDM_model_plot=
  function(dis_measure_function,plot_title,weights=1,parameters){
    RDM_model = df_to_RDM(df_trials,dis_measure_function,
                          weights,
                          parameters = parameters)
    melt_RDM =melt(RDM_model)

    g=melt_RDM%>%
      ggplot(aes(x=Var1,y=Var2,fill=value))+
      geom_tile()+
      scale_y_reverse(breaks = c(AMOUNT_TRIALS_EACH/2,(2.5*AMOUNT_TRIALS_EACH)/2, # 2.5 & 3.5 to adjust the position of the axis ticks
                                 (3.5*AMOUNT_TRIALS_EACH)/2,(5*AMOUNT_TRIALS_EACH)/2),
                      labels = c("FF", "FL","LF", "LL"))+
      theme(plot.title = element_text(hjust = 0.5))+
      scale_x_continuous(breaks = c(AMOUNT_TRIALS_EACH/2,(2.5*AMOUNT_TRIALS_EACH)/2,
                                    (3.5*AMOUNT_TRIALS_EACH)/2,(5*AMOUNT_TRIALS_EACH)/2),
                         labels = c("FF", "FL","LF", "LL"))+
      labs(title=plot_title,
           x="trials",y="trials",fill="dissimilarity")+
      theme(axis.text=element_text(size=14))
    return(list(plot=g,RDM=RDM_model))
  }
#Function to scale RDM with any unit to the scale of 0~2---------------------------
scale_RSA=
  function(variable,new_min,new_max){
    #variable: could be a vector, a matrix , or an dataframe

    old_min=min(variable)
    old_max=max(variable)
    old_range=old_max-old_min
    new_range=new_max-new_min

    shift=new_min-old_min
    scaling=new_range/old_range

    scaled_variable=(variable+shift)*(scaling)
    return(scaled_variable)
  }
#Preprocess the data-------------------
df_trials=
  read.csv("D:\\Yun-Shiuan_LAMBDA\\EprimeData_raw\\DevFracs1002\\df1002_Run2_tidy.csv")%>%
  filter(!is.na(trial_id))%>%
  select(fraccomp_sti_type,fraccomp_sti_dis_type,
         fraccomp_sti_value_left,fraccomp_sti_value_right,fraccomp_sti_dis)%>%
  mutate(fraccomp_sti_dis_abs=abs(fraccomp_sti_dis))%>%
  mutate(F_left=grepl("^(f|F)",fraccomp_sti_type),
         F_right=grepl("(f|F)$",fraccomp_sti_type),
         L_left=grepl("^(l|L)",fraccomp_sti_type),
         L_right=grepl("(l|L)$",fraccomp_sti_type))%>%
  mutate(F_amount=F_left+F_right,
         L_amount=L_left+L_right)%>%
  arrange(desc(F_amount),F_left,fraccomp_sti_dis_abs)%>%
  mutate(fraccomp_sti_dis_abs_scaled=scale_RSA(fraccomp_sti_dis_abs,0,2),
         fraccomp_sti_dis_scaled=scale_RSA(fraccomp_sti_dis,0,2))


#Model-1: Format with location information(no distance information)-------------------------------------
dist_format_with_location_no_distance=
  function(trial1_content,trial2_content,weights,parameters){
    trial1_features=
      trial1_content%>%
      select(F_left,F_right,L_left,L_right)

    trial2_features=
      trial2_content%>%
      select(F_left,F_right,L_left,L_right)
    return(sum(abs(trial1_features-trial2_features))/2)
  }

plot_title="Format with location information \n (no distance information)"
RDM_model1_format_with_location =
  RDM_model_plot(dis_measure_function = dist_format_with_location_no_distance,
                 plot_title = plot_title )
print(RDM_model1_format_with_location$plot)
#Model-2: Format without location information(no distance information)-------------------------------------
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
plot_title="Format without location information \n (no distance information)"
RDM_model2_format_without_location =
  RDM_model_plot(dis_measure_function = dist_format_without_location_no_distance,
                 plot_title = plot_title )
print(RDM_model2_format_without_location$plot)
#Model-3: Format with location information(with absolute distance information)-------------------------------------
#w is The relative weight of distance difference compared to format difference
dist_format_with_location_abs_dist=
  function(trial1_content,trial2_content,parameters){
    w=parameters[1]# The coeffient of how much distance diff. weigh compare to format diff.
    trial1_features=
      trial1_content%>%
      select(F_left,F_right,L_left,L_right,fraccomp_sti_dis_abs)

    trial2_features=
      trial2_content%>%
      select(F_left,F_right,L_left,L_right,fraccomp_sti_dis_abs)
    format_diff=sum(abs(trial1_features-trial2_features))/2#Scale 0~2 (no need to scale)
    dist_diff=abs(trial1_features["fraccomp_sti_dis_abs"]-trial2_features["fraccomp_sti_dis_abs"])
    return(as.numeric((1-w)*format_diff+(w*dist_diff)))
  }
#w=0.2-------------------------------------------------------
w=0.2 # The relative weight of distance difference compared to format difference
plot_title = paste0("Format with location information \n (absolute dist.; weight(dist) = ",w,")")
RDM_model3_format_with_location_abs_dist_k0.2=
  RDM_model_plot(dis_measure_function = dist_format_with_location_abs_dist,
                 plot_title = plot_title,
                 parameters = w )
print(RDM_model3_format_with_location_abs_dist_k0.2$plot)

#w=0.5-------------------------------------------------------
w=0.5 # The relative weight of distance difference compared to format difference
plot_title = paste0("Format with location information \n (absolute dist.; weight(dist) = ",w,")")
RDM_model3_format_with_location_abs_dist_k0.5=
  RDM_model_plot(dis_measure_function = dist_format_with_location_abs_dist,
                 plot_title = plot_title,
                 parameters = w )
print(RDM_model3_format_with_location_abs_dist_k0.5$plot)

#w=0.8-------------------------------------------------------
w=0.8 # The relative weight of distance difference compared to format difference
plot_title = paste0("Format with location information \n (absolute dist.; weight(dist) = ",w,")")
RDM_model3_format_with_location_abs_dist_k0.8=
  RDM_model_plot(dis_measure_function = dist_format_with_location_abs_dist,
                 plot_title = plot_title,
                 parameters = w )
print(RDM_model3_format_with_location_abs_dist_k0.8$plot)

#Model-4: Format without location information(with absolute distance information)-------------------------------------
#w is The relative weight of distance difference compared to format difference
dist_format_without_location_abs_dist=
  function(trial1_content,trial2_content,parameters){
    w=parameters[1]# The coeffient of how much distance diff. weigh compare to format diff.
    trial1_features=
      trial1_content%>%
      select(F_amount,L_amount,fraccomp_sti_dis_abs)

    trial2_features=
      trial2_content%>%
      select(F_amount,L_amount,fraccomp_sti_dis_abs)
    format_diff=sum(abs(trial1_features[c("F_amount","L_amount")]-trial2_features[c("F_amount","L_amount")]))/2
    dist_diff=abs(trial1_features["fraccomp_sti_dis_abs"]-trial2_features["fraccomp_sti_dis_abs"])
    return(as.numeric((1-w)*format_diff+(w*dist_diff)))
  }
#w=0.2-------------------------------------------------------
w=0.2 # The relative weight of distance difference compared to format difference
plot_title=paste0("Format without location information \n (absolute dist.; weight(dist) = ",w,")")
RDM_model4_format_without_location_abs_dist_k0.2=
  RDM_model_plot(dis_measure_function = dist_format_without_location_abs_dist,
                 plot_title = plot_title,
                 parameters = w)
print(RDM_model4_format_without_location_abs_dist_k0.2$plot)
#w=0.5-------------------------------------------------------
w=0.5 # The relative weight of distance difference compared to format difference
plot_title=paste0("Format without location information \n (absolute dist.; weight(dist) = ",w,")")
RDM_model4_format_without_location_abs_dist_k0.5=
  RDM_model_plot(dis_measure_function = dist_format_without_location_abs_dist,
                 plot_title = plot_title,
                 parameters = w)
print(RDM_model4_format_without_location_abs_dist_k0.5$plot)
#w=0.8-------------------------------------------------------
w=0.8 # The relative weight of distance difference compared to format difference
plot_title=paste0("Format without location information \n (absolute dist.; weight(dist) = ",w,")")
RDM_model4_format_without_location_abs_dist_k3=
  RDM_model_plot(dis_measure_function = dist_format_without_location_abs_dist,
                 plot_title = plot_title,
                 parameters = w)
print(RDM_model4_format_without_location_abs_dist_k3$plot)

#Model-5: absolute distance information (no format or position)-------------------------------------
dist_only_abs_dist=
  function(trial1_content,trial2_content,parameters){
    trial1_features=
      trial1_content%>%
      select(fraccomp_sti_dis_abs)

    trial2_features=
      trial2_content%>%
      select(fraccomp_sti_dis_abs)
    dist_diff=abs(trial1_features["fraccomp_sti_dis_abs"]-trial2_features["fraccomp_sti_dis_abs"])
    return(as.numeric(dist_diff))
  }

plot_title="Absolute distance \n (No format nor position)"
RDM_model5_dist_only_abs_dist=RDM_model_plot(dis_measure_function = dist_only_abs_dist,
                                             plot_title = plot_title)
print(RDM_model5_dist_only_abs_dist$plot)
#Model-6: signed distance information (no format nor position)-------------------------------------
dist_only_signed_dist=
  function(trial1_content,trial2_content,parameters){
    trial1_features=
      trial1_content%>%
      select(fraccomp_sti_dis)

    trial2_features=
      trial2_content%>%
      select(fraccomp_sti_dis)
    dist_diff=abs(trial1_features["fraccomp_sti_dis"]-trial2_features["fraccomp_sti_dis"])
    return(as.numeric(dist_diff))
  }
plot_title="Signed distance \n (No format nor position)"
RDM_model6_dist_only_signed_dist=RDM_model_plot(dis_measure_function = dist_only_signed_dist,
                                                plot_title = plot_title)
print(RDM_model6_dist_only_signed_dist$plot)

#Model-7: Format with location information(with absolute distance information)-------------------------------------
#w is The relative weight of distance difference compared to format difference
dist_format_with_location_signed_dist=
  function(trial1_content,trial2_content,parameters){
    w=parameters[1]# The coeffient of how much distance diff. weigh compare to format diff.
    trial1_features=
      trial1_content%>%
      select(F_left,F_right,L_left,L_right,fraccomp_sti_dis)

    trial2_features=
      trial2_content%>%
      select(F_left,F_right,L_left,L_right,fraccomp_sti_dis)
    format_diff=sum(abs(trial1_features-trial2_features))/2
    dist_diff=abs(trial1_features["fraccomp_sti_dis"]-trial2_features["fraccomp_sti_dis"])
    return(as.numeric((1-w)*format_diff+(w*dist_diff)))
  }
#w=0.2-------------------------------------------------------
w=0.2 # The relative weight of distance difference compared to format difference
plot_title=paste0("Format with location information \n (signed dist.; weight(dist) = ",w,")")
RDM_model_format_with_location_signed_dist_k0.2=
  RDM_model_plot(dis_measure_function = dist_format_with_location_signed_dist,
                 plot_title = plot_title,parameters = w)
print(RDM_model_format_with_location_signed_dist_k0.2$plot)
#w=0.5-------------------------------------------------------
w=0.5 # The relative weight of distance difference compared to format difference
plot_title=paste0("Format with location information \n (signed dist.; weight(dist) = ",w,")")
RDM_model_format_with_location_signed_dist_k0.5=
  RDM_model_plot(dis_measure_function = dist_format_with_location_signed_dist,
                 plot_title = plot_title,parameters = w)
print(RDM_model_format_with_location_signed_dist_k0.5$plot)
#w=0.8-------------------------------------------------------
w=0.8 # The relative weight of distance difference compared to format difference
plot_title=paste0("Format with location information \n (signed dist.; weight(dist) = ",w,")")
RDM_model_format_with_location_signed_dist_k0.8=
  RDM_model_plot(dis_measure_function = dist_format_with_location_signed_dist,
                 plot_title = plot_title,parameters = w)
print(RDM_model_format_with_location_signed_dist_k0.8$plot)

#Model-8: Format without location information(with absolute distance information)-------------------------------------
#w is The relative weight of distance difference compared to format difference
dist_format_without_location_signed_dist=
  function(trial1_content,trial2_content,parameters){
    w=parameters[1]# The coeffient of how much distance diff. weigh compare to format diff.
    trial1_features=
      trial1_content%>%
      select(F_amount,L_amount,fraccomp_sti_dis)

    trial2_features=
      trial2_content%>%
      select(F_amount,L_amount,fraccomp_sti_dis)
    format_diff=sum(abs(trial1_features[c("F_amount","L_amount")]-trial2_features[c("F_amount","L_amount")]))/2
    dist_diff=abs(trial1_features["fraccomp_sti_dis"]-trial2_features["fraccomp_sti_dis"])
    return(as.numeric((1-w)*format_diff+(w*dist_diff)))
  }
#w=0.2-------------------------------------------------------
w=0.2 # The relative weight of distance difference compared to format difference
plot_title=paste0("Format without location information \n (absolute dist.; weight(dist) = ",w,")")
RDM_model_format_without_location_signed_dist_k0.2=
  RDM_model_plot(dis_measure_function = dist_format_without_location_signed_dist,
                 plot_title = plot_title,parameters = w)
print(RDM_model_format_without_location_signed_dist_k0.2$plot)
#w=0.5-------------------------------------------------------
w=0.5 # The relative weight of distance difference compared to format difference
plot_title=paste0("Format without location information \n (absolute dist.; weight(dist) = ",w,")")
RDM_model_format_without_location_signed_dist_k0.5=
  RDM_model_plot(dis_measure_function = dist_format_without_location_signed_dist,
                 plot_title = plot_title,parameters = w)
print(RDM_model_format_without_location_signed_dist_k0.5$plot)
#w=0.8-------------------------------------------------------
w=0.8 # The relative weight of distance difference compared to format difference
plot_title=paste0("Format without location information \n (absolute dist.; weight(dist) = ",w,")")
RDM_model_format_without_location_signed_dist_k0.8=
  RDM_model_plot(dis_measure_function = dist_format_without_location_signed_dist,
                 plot_title = plot_title,parameters = w)
print(RDM_model_format_without_location_signed_dist_k0.8$plot)
