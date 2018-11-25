#Representational Similarity Analysis(RSA)[Visualize Conceptual model]----------------------------------
#Note that all functions were moved to a separate R script: see the constant FILE_RSA_FUNCTIONS
library("dplyr")
library("tidyr")
library("stringr")
library("ggplot2")
library("reshape2")
#NOTE:
# - Continuous: Regard distance as continuous scale (ID: 1~144).
#It is impratical due to the limit of RSA toolbox, which requires all runs to have the same set of trials.
# - Discrete with 24 trials:  Regard distance as discrete scale (ID: 1~24)
#24 unique trials = 4 formats x 2 direction x 3 distance levels.
# - Discrete with 36 trials:
# Use RSA ID 1~36 so that each ID refers to a trial in a run.
# Note Each trial in a run has an ID, so that there'll be 36 IDs
# 	      (note there're only 24 unique trials: 4 formats x 2 locations x 3 distance levels)
# [FL: 6 trials = 3 distance levels x 2 direction]
# [LF: 6 trials = 3 distance levels x 2 direction]
# [LL: 12 trials = 3 distance levels x 2 direction x 2 repetition]
# [FF: 12 trials = 3 distance levels x 2 direction x 2 repetition]
#
# The motivation for having 36 trials is to ensure beta estimate of each condition
# is subjected to same amount of noise.
# (Having more repetition in a condition entails a smaller degree of random noise.)
############################################################################
# Re-define the scale: set w from 0.5 to 1.0 (Trial 19, 2018/6/7)


#Declare constants----------------
#Parameters
# AMOUNT_TRIALS_TOTAL=36 # Need to be modified depends on using 24ID(each for a unique trial) or 36ID(each for a trial in a run) version
# AMOUNT_FORMAT=3# Need to be modified depends on using 24ID(each for a unique trial) or 36ID(each for a trial in a run) version
# AMOUNT_TRIALS_EACH=AMOUNT_TRIALS_TOTAL/AMOUNT_FORMAT
RDM_MIN=0
RDM_MAX=1
#Paths
PATH_ROOT="D:\\Yun-Shiuan_LAMBDA"

#Files
FILE_RSA_FUNCTIONS="D:\\GoogleDrive\\Lambda_code\\R_script\\R_functions\\RSA_conceptual_models_functions.R"
# FILE_RSA_FUNCTIONS="/study3/devfracs/DCM_YunShiuan/Lambda_code/R_script/R_functions/RSA_conceptual_models_functions.R"
FILE_BEH_DATA=file.path(PATH_ROOT,"EprimeData_raw","For_RSA","df_for_RDM_construction.csv")

# FILE_BEH_DATA="/study3/devfracs/DCM_YunShiuan/EprimeData_raw/For_RSA/df_for_RDM_construction.csv"
#For hybrid models: w is The relative weight of distance difference compared to format difference
# LIST_W=c(0.2,0.5,0.8) #(Before re-define the scale)
WEIGHTS=seq(from=0.5,to = 0.99,by=0.005) # The weight range for mixed dist measure and file names
names(WEIGHTS)=as.character(1:99)
LIST_W=WEIGHTS[c(20,40,60,80)]
# Load in RSA related functions--------------------
# Include dissimilarity definitions and helper functions
source(FILE_RSA_FUNCTIONS)
# dataframe for RDM derivation----------------------------
df_trials_continuous_24=preprocess_df_continuous_unique_trails_24(FILE_BEH_DATA) # Sort the trials by RSA ID derived from continuous distance (0~1)
df_trials_discrete_24=preprocess_df_discrete_unique_trails_24(FILE_BEH_DATA) # Sort the trials by RSA ID derived from discrete distance (N/M/F)
df_trials_discrete_36=preprocess_df_discrete_repeated_trails_36(FILE_BEH_DATA) # Sort the trials by RSA ID derived from discrete distance (N/M/F)
df_trials_discrete_18=preprocess_df_discrete_repeated_trails_18(FILE_BEH_DATA) # Sort the trials by RSA ID derived from discrete distance (N/M/F)

# Continuous distance RSA ID=======================================)=========================================================
# NOTE: For the sake of retrieval convenience , do NOT wrap the models below in a loop.
#Model-1: Format with location information(no distance information)-------------------------------------
plot_title="Format with location information \n (no distance information)"
RDM_model1_format_with_location =
  RDM_model_plot(df = df_trials_continuous_24,
                 dis_measure_function = dist_format_with_location_no_distance,
                 plot_title = plot_title )
print(RDM_model1_format_with_location$plot)
#Model-2: Format without location information(no distance information)-------------------------------------
plot_title="Format without location information \n (no distance information)"
RDM_model2_format_without_location =
  RDM_model_plot(df = df_trials_continuous_24,
                 dis_measure_function = dist_format_without_location_no_distance,
                 plot_title = plot_title )
print(RDM_model2_format_without_location$plot)
#Model-5: absolute distance information (no format or position)-------------------------------------
plot_title="Absolute distance \n (No format or position)"
RDM_model5_dist_only_abs_dist=
  RDM_model_plot(df = df_trials_continuous_24,
                 dis_measure_function = dist_only_abs_dist,
                 plot_title = plot_title)
print(RDM_model5_dist_only_abs_dist$plot)

#Model-3: Format with location information(with absolute distance information)-------------------------------------
#w is The relative weight of distance difference compared to format difference
#w=0.2-------------------------------------------------------
w=0.2 # The relative weight of distance difference compared to format difference
plot_title = paste0("Format with location information \n (absolute dist.; weight(dist) = ",w,")")
RDM_model3_format_with_location_abs_dist_k0.2=
  RDM_model_plot(df = df_trials_continuous_24,
                 dis_measure_function = dist_format_with_location_abs_dist,
                 plot_title = plot_title,
                 parameters = w )
print(RDM_model3_format_with_location_abs_dist_k0.2$plot)
#w=0.5-------------------------------------------------------
w=0.5 # The relative weight of distance difference compared to format difference
plot_title = paste0("Format with location information \n (absolute dist.; weight(dist) = ",w,")")
RDM_model3_format_with_location_abs_dist_k0.5=
  RDM_model_plot(df = df_trials_continuous_24,
                 dis_measure_function = dist_format_with_location_abs_dist,
                 plot_title = plot_title,
                 parameters = w )
print(RDM_model3_format_with_location_abs_dist_k0.5$plot)
#w=0.8-------------------------------------------------------
w=0.8 # The relative weight of distance difference compared to format difference
plot_title = paste0("Format with location information \n (absolute dist.; weight(dist) = ",w,")")
RDM_model3_format_with_location_abs_dist_k0.8=
  RDM_model_plot(df = df_trials_continuous_24,
                 dis_measure_function = dist_format_with_location_abs_dist,
                 plot_title = plot_title,
                 parameters = w )
print(RDM_model3_format_with_location_abs_dist_k0.8$plot)

#Model-4: Format without location information(with absolute distance information)-------------------------------------
#w is The relative weight of distance difference compared to format difference
#w=0.2-------------------------------------------------------
w=0.2 # The relative weight of distance difference compared to format difference
plot_title=paste0("Format without location information \n (absolute dist.; weight(dist) = ",w,")")
RDM_model4_format_without_location_abs_dist_k0.2=
  RDM_model_plot(df = df_trials_continuous_24,
                 dis_measure_function = dist_format_without_location_abs_dist,
                 plot_title = plot_title,
                 parameters = w)
print(RDM_model4_format_without_location_abs_dist_k0.2$plot)
#w=0.5-------------------------------------------------------
w=0.5 # The relative weight of distance difference compared to format difference
plot_title=paste0("Format without location information \n (absolute dist.; weight(dist) = ",w,")")
RDM_model4_format_without_location_abs_dist_k0.5=
  RDM_model_plot(df = df_trials_continuous_24,
                 dis_measure_function = dist_format_without_location_abs_dist,
                 plot_title = plot_title,
                 parameters = w)
print(RDM_model4_format_without_location_abs_dist_k0.5$plot)
#w=0.8-------------------------------------------------------
w=0.8 # The relative weight of distance difference compared to format difference
plot_title=paste0("Format without location information \n (absolute dist.; weight(dist) = ",w,")")
RDM_model4_format_without_location_abs_dist_k3=
  RDM_model_plot(df = df_trials_continuous_24,
                 dis_measure_function = dist_format_without_location_abs_dist,
                 plot_title = plot_title,
                 parameters = w)
print(RDM_model4_format_without_location_abs_dist_k3$plot)
#Model-6: signed distance information (no format nor position)-------------------------------------
plot_title="Signed distance \n (No format or position)"
RDM_model6_dist_only_signed_dist=
  RDM_model_plot(df = df_trials_continuous_24,
                 dis_measure_function = dist_only_signed_dist,
                 plot_title = plot_title)
print(RDM_model6_dist_only_signed_dist$plot)

#Model-7: Format with location information(with signed distance information)-------------------------------------
#w is The relative weight of distance difference compared to format difference

#w=0.2-------------------------------------------------------
w=0.2 # The relative weight of distance difference compared to format difference
plot_title=paste0("Format with location information \n (signed dist.; weight(dist) = ",w,")")
RDM_model_format_with_location_signed_dist_k0.2=
  RDM_model_plot(df = df_trials_continuous_24,
                 dis_measure_function = dist_format_with_location_signed_dist,
                 plot_title = plot_title,parameters = w)
print(RDM_model_format_with_location_signed_dist_k0.2$plot)
#w=0.5-------------------------------------------------------
w=0.5 # The relative weight of distance difference compared to format difference
plot_title=paste0("Format with location information \n (signed dist.; weight(dist) = ",w,")")
RDM_model_format_with_location_signed_dist_k0.5=
  RDM_model_plot(df = df_trials_continuous_24,
                 dis_measure_function = dist_format_with_location_signed_dist,
                 plot_title = plot_title,parameters = w)
print(RDM_model_format_with_location_signed_dist_k0.5$plot)
#w=0.8-------------------------------------------------------
w=0.8 # The relative weight of distance difference compared to format difference
plot_title=paste0("Format with location information \n (signed dist.; weight(dist) = ",w,")")
RDM_model_format_with_location_signed_dist_k0.8=
  RDM_model_plot(df = df_trials_continuous_24,
                 dis_measure_function = dist_format_with_location_signed_dist,
                 plot_title = plot_title,parameters = w)
print(RDM_model_format_with_location_signed_dist_k0.8$plot)

#Model-8: Format without location information(with signed distance information)-------------------------------------
#w is The relative weight of distance difference compared to format difference

#w=0.2-------------------------------------------------------
w=0.2 # The relative weight of distance difference compared to format difference
plot_title=paste0("Format without location information \n (signed dist.; weight(dist) = ",w,")")
RDM_model_format_without_location_signed_dist_k0.2=
  RDM_model_plot(df = df_trials_continuous_24,
                 dis_measure_function = dist_format_without_location_signed_dist,
                 plot_title = plot_title,parameters = w)
print(RDM_model_format_without_location_signed_dist_k0.2$plot)
#w=0.5-------------------------------------------------------
w=0.5 # The relative weight of distance difference compared to format difference
plot_title=paste0("Format without location information \n (signed dist.; weight(dist) = ",w,")")
RDM_model_format_without_location_signed_dist_k0.5=
  RDM_model_plot(df = df_trials_continuous_24,
                 dis_measure_function = dist_format_without_location_signed_dist,
                 plot_title = plot_title,parameters = w)
print(RDM_model_format_without_location_signed_dist_k0.5$plot)
#w=0.8-------------------------------------------------------
w=0.8 # The relative weight of distance difference compared to format difference
plot_title=paste0("Format without location information \n (signed dist.; weight(dist) = ",w,")")
RDM_model_format_without_location_signed_dist_k0.8=
  RDM_model_plot(df = df_trials_continuous_24,
                 dis_measure_function = dist_format_without_location_signed_dist,
                 plot_title = plot_title,parameters = w)
print(RDM_model_format_without_location_signed_dist_k0.8$plot)

# Discrete distance RSA ID (unique ID: 1~24)==========================================)======================================================
#Model-1: Format with location information(no distance information)-------------------------------------
plot_title="Format with location information \n (no distance information)"
RDM_model1_format_with_location =
  RDM_model_plot(df = df_trials_discrete_24,
                 dis_measure_function = dist_format_with_location_no_distance,
                 plot_title = plot_title )
print(RDM_model1_format_with_location$plot)
#Model-2: Format without location information(no distance information)-------------------------------------
plot_title="Format without location information \n (no distance information)"
RDM_model2_format_without_location =
  RDM_model_plot(df = df_trials_discrete_24,
                 dis_measure_function = dist_format_without_location_no_distance,
                 plot_title = plot_title )
print(RDM_model2_format_without_location$plot)
#Model-5: absolute distance information (no format or position)-------------------------------------
plot_title="Absolute distance \n (No format or position)"
RDM_model5_dist_only_abs_dist_discrete=
  RDM_model_plot(df = df_trials_discrete_24,
                 dis_measure_function = dist_only_abs_dist_discrete,
                 plot_title = plot_title)
print(RDM_model5_dist_only_abs_dist_discrete$plot)

#Model-3: Format with location information(with absolute distance information)-------------------------------------
#w is The relative weight of distance difference compared to format difference
#w=0.2-------------------------------------------------------
w=0.2 # The relative weight of distance difference compared to format difference
plot_title = paste0("Format with location information \n (absolute dist.; weight(dist) = ",w,")")
RDM_model3_format_with_location_abs_dist_k0.2=
  RDM_model_plot(df = df_trials_discrete_24,
                 dis_measure_function = dist_format_with_location_abs_dist_discrete,
                 plot_title = plot_title,
                 parameters = w )
print(RDM_model3_format_with_location_abs_dist_k0.2$plot)
#w=0.5-------------------------------------------------------
w=0.5 # The relative weight of distance difference compared to format difference
plot_title = paste0("Format with location information \n (absolute dist.; weight(dist) = ",w,")")
RDM_model3_format_with_location_abs_dist_k0.5=
  RDM_model_plot(df = df_trials_discrete_24,
                 dis_measure_function = dist_format_with_location_abs_dist_discrete,
                 plot_title = plot_title,
                 parameters = w )
print(RDM_model3_format_with_location_abs_dist_k0.5$plot)

#w=0.8-------------------------------------------------------
w=0.8 # The relative weight of distance difference compared to format difference
plot_title = paste0("Format with location information \n (absolute dist.; weight(dist) = ",w,")")
RDM_model3_format_with_location_abs_dist_k0.8=
  RDM_model_plot(df = df_trials_discrete_24,
                 dis_measure_function = dist_format_with_location_abs_dist_discrete,
                 plot_title = plot_title,
                 parameters = w )
print(RDM_model3_format_with_location_abs_dist_k0.8$plot)

#Model-4: Format without location information(with absolute distance information)-------------------------------------
#w is The relative weight of distance difference compared to format difference
#w=0.2-------------------------------------------------------
w=0.2 # The relative weight of distance difference compared to format difference
plot_title=paste0("Format without location information \n (absolute dist.; weight(dist) = ",w,")")
RDM_model4_format_without_location_abs_dist_k0.2=
  RDM_model_plot(df = df_trials_discrete_24,
                 dis_measure_function = dist_format_without_location_abs_dist_discrete,
                 plot_title = plot_title,
                 parameters = w)
print(RDM_model4_format_without_location_abs_dist_k0.2$plot)
#w=0.5-------------------------------------------------------
w=0.5 # The relative weight of distance difference compared to format difference
plot_title=paste0("Format without location information \n (absolute dist.; weight(dist) = ",w,")")
RDM_model4_format_without_location_abs_dist_k0.5=
  RDM_model_plot(df = df_trials_discrete_24,
                 dis_measure_function = dist_format_without_location_abs_dist_discrete,
                 plot_title = plot_title,
                 parameters = w)
print(RDM_model4_format_without_location_abs_dist_k0.5$plot)
#w=0.8-------------------------------------------------------
w=0.8 # The relative weight of distance difference compared to format difference
plot_title=paste0("Format without location information \n (absolute dist.; weight(dist) = ",w,")")
RDM_model4_format_without_location_abs_dist_k3=
  RDM_model_plot(df = df_trials_discrete_24,
                 dis_measure_function = dist_format_without_location_abs_dist_discrete,
                 plot_title = plot_title,
                 parameters = w)
print(RDM_model4_format_without_location_abs_dist_k3$plot)

#Model-6: signed distance information (no format nor position)-------------------------------------
plot_title="Signed distance \n (No format or position)"
RDM_model6_dist_only_signed_dist=
  RDM_model_plot(df = df_trials_discrete_24,
                 dis_measure_function = dist_only_signed_dist_discrete,
                 plot_title = plot_title)
print(RDM_model6_dist_only_signed_dist$plot)

#Model-7: Format with location information(with signed distance information)-------------------------------------
#w is The relative weight of distance difference compared to format difference

#w=0.2-------------------------------------------------------
w=0.2 # The relative weight of distance difference compared to format difference
plot_title=paste0("Format with location information \n (signed dist.; weight(dist) = ",w,")")
RDM_model_format_with_location_signed_dist_k0.2=
  RDM_model_plot(df = df_trials_discrete_24,
                 dis_measure_function = dist_format_with_location_signed_dist_discrete,
                 plot_title = plot_title,parameters = w)
print(RDM_model_format_with_location_signed_dist_k0.2$plot)
#w=0.5-------------------------------------------------------
w=0.5 # The relative weight of distance difference compared to format difference
plot_title=paste0("Format with location information \n (signed dist.; weight(dist) = ",w,")")
RDM_model_format_with_location_signed_dist_k0.5=
  RDM_model_plot(df = df_trials_discrete_24,
                 dis_measure_function = dist_format_with_location_signed_dist_discrete,
                 plot_title = plot_title,parameters = w)
print(RDM_model_format_with_location_signed_dist_k0.5$plot)
#w=0.8-------------------------------------------------------
w=0.8 # The relative weight of distance difference compared to format difference
plot_title=paste0("Format with location information \n (signed dist.; weight(dist) = ",w,")")
RDM_model_format_with_location_signed_dist_k0.8=
  RDM_model_plot(df = df_trials_discrete_24,
                 dis_measure_function = dist_format_with_location_signed_dist_discrete,
                 plot_title = plot_title,parameters = w)
print(RDM_model_format_with_location_signed_dist_k0.8$plot)

#Model-8: Format without location information(with signed distance information)-------------------------------------
#w is The relative weight of distance difference compared to format difference

#w=0.2-------------------------------------------------------
w=0.2 # The relative weight of distance difference compared to format difference
plot_title=paste0("Format without location information \n (signed dist.; weight(dist) = ",w,")")
RDM_model_format_without_location_signed_dist_k0.2=
  RDM_model_plot(df = df_trials_discrete_24,
                 dis_measure_function = dist_format_without_location_signed_dist_discrete,
                 plot_title = plot_title,parameters = w)
print(RDM_model_format_without_location_signed_dist_k0.2$plot)
#w=0.5-------------------------------------------------------
w=0.5 # The relative weight of distance difference compared to format difference
plot_title=paste0("Format without location information \n (signed dist.; weight(dist) = ",w,")")
RDM_model_format_without_location_signed_dist_k0.5=
  RDM_model_plot(df = df_trials_discrete_24,
                 dis_measure_function = dist_format_without_location_signed_dist_discrete,
                 plot_title = plot_title,parameters = w)
print(RDM_model_format_without_location_signed_dist_k0.5$plot)
#w=0.8-------------------------------------------------------
w=0.8 # The relative weight of distance difference compared to format difference
plot_title=paste0("Format without location information \n (signed dist.; weight(dist) = ",w,")")
RDM_model_format_without_location_signed_dist_k0.8=
  RDM_model_plot(df = df_trials_discrete_24,
                 dis_measure_function = dist_format_without_location_signed_dist_discrete,
                 plot_title = plot_title,parameters = w)
print(RDM_model_format_without_location_signed_dist_k0.8$plot)

# Discrete distance RSA ID (unique ID: 1~36)==========================================)======================================================
#Model-1: Format with location information(no distance information)-------------------------------------
plot_title="Format with location information \n (no distance information)"
RDM_model1_format_with_location =
  RDM_model_plot(df = df_trials_discrete_36,
                 dis_measure_function = dist_format_with_location_no_distance,
                 plot_title = plot_title )
print(RDM_model1_format_with_location$plot)
#Model-2: Format without location information(no distance information)-------------------------------------
plot_title="Format without location information \n (no distance information)"
RDM_model2_format_without_location =
  RDM_model_plot(df = df_trials_discrete_36,
                 dis_measure_function = dist_format_without_location_no_distance,
                 plot_title = plot_title )
print(RDM_model2_format_without_location$plot)
#Model-5: absolute distance information (no format or position)-------------------------------------
plot_title="Absolute distance \n (No format or position)"
RDM_model5_dist_only_abs_dist_discrete=
  RDM_model_plot(df = df_trials_discrete_36,
                 dis_measure_function = dist_only_abs_dist_discrete,
                 plot_title = plot_title)
print(RDM_model5_dist_only_abs_dist_discrete$plot)

#Model-3: Format with location information(with absolute distance information)-------------------------------------
#w is The relative weight of distance difference compared to format difference
#w=0.2-------------------------------------------------------
w=0.2 # The relative weight of distance difference compared to format difference
plot_title = paste0("Format with location information \n (absolute dist.; weight(dist) = ",w,")")
RDM_model3_format_with_location_abs_dist_k0.2=
  RDM_model_plot(df = df_trials_discrete_36,
                 dis_measure_function = dist_format_with_location_abs_dist_discrete,
                 plot_title = plot_title,
                 parameters = w )
print(RDM_model3_format_with_location_abs_dist_k0.2$plot)
#w=0.5-------------------------------------------------------
w=0.5 # The relative weight of distance difference compared to format difference
plot_title = paste0("Format with location information \n (absolute dist.; weight(dist) = ",w,")")
RDM_model3_format_with_location_abs_dist_k0.5=
  RDM_model_plot(df = df_trials_discrete_36,
                 dis_measure_function = dist_format_with_location_abs_dist_discrete,
                 plot_title = plot_title,
                 parameters = w )
print(RDM_model3_format_with_location_abs_dist_k0.5$plot)
#w=0.8-------------------------------------------------------
w=0.8 # The relative weight of distance difference compared to format difference
plot_title = paste0("Format with location information \n (absolute dist.; weight(dist) = ",w,")")
RDM_model3_format_with_location_abs_dist_k0.8=
  RDM_model_plot(df = df_trials_discrete_36,
                 dis_measure_function = dist_format_with_location_abs_dist_discrete,
                 plot_title = plot_title,
                 parameters = w )
print(RDM_model3_format_with_location_abs_dist_k0.8$plot)

#Model-4: Format without location information(with absolute distance information)-------------------------------------
#w is The relative weight of distance difference compared to format difference
#w=0.2-------------------------------------------------------
w=0.2 # The relative weight of distance difference compared to format difference
plot_title=paste0("Format without location information \n (absolute dist.; weight(dist) = ",w,")")
RDM_model4_format_without_location_abs_dist_k0.2=
  RDM_model_plot(df = df_trials_discrete_36,
                 dis_measure_function = dist_format_without_location_abs_dist_discrete,
                 plot_title = plot_title,
                 parameters = w)
print(RDM_model4_format_without_location_abs_dist_k0.2$plot)
#w=0.5-------------------------------------------------------
w=0.5 # The relative weight of distance difference compared to format difference
plot_title=paste0("Format without location information \n (absolute dist.; weight(dist) = ",w,")")
RDM_model4_format_without_location_abs_dist_k0.5=
  RDM_model_plot(df = df_trials_discrete_36,
                 dis_measure_function = dist_format_without_location_abs_dist_discrete,
                 plot_title = plot_title,
                 parameters = w)
print(RDM_model4_format_without_location_abs_dist_k0.5$plot)
#w=0.8-------------------------------------------------------
w=0.8 # The relative weight of distance difference compared to format difference
plot_title=paste0("Format without location information \n (absolute dist.; weight(dist) = ",w,")")
RDM_model4_format_without_location_abs_dist_k3=
  RDM_model_plot(df = df_trials_discrete_36,
                 dis_measure_function = dist_format_without_location_abs_dist_discrete,
                 plot_title = plot_title,
                 parameters = w)
print(RDM_model4_format_without_location_abs_dist_k3$plot)

#Model-6: signed distance information (no format nor position)-------------------------------------
plot_title="Signed distance \n (No format or position)"
RDM_model6_dist_only_signed_dist=
  RDM_model_plot(df = df_trials_discrete_36,
                 dis_measure_function = dist_only_signed_dist_discrete,
                 plot_title = plot_title)
print(RDM_model6_dist_only_signed_dist$plot)

#Model-7: Format with location information(with signed distance information)-------------------------------------
#w is The relative weight of distance difference compared to format difference

#w=0.2-------------------------------------------------------
w=0.2 # The relative weight of distance difference compared to format difference
plot_title=paste0("Format with location information \n (signed dist.; weight(dist) = ",w,")")
RDM_model_format_with_location_signed_dist_k0.2=
  RDM_model_plot(df = df_trials_discrete_36,
                 dis_measure_function = dist_format_with_location_signed_dist_discrete,
                 plot_title = plot_title,parameters = w)
print(RDM_model_format_with_location_signed_dist_k0.2$plot)
#w=0.5-------------------------------------------------------
w=0.5 # The relative weight of distance difference compared to format difference
plot_title=paste0("Format with location information \n (signed dist.; weight(dist) = ",w,")")
RDM_model_format_with_location_signed_dist_k0.5=
  RDM_model_plot(df = df_trials_discrete_36,
                 dis_measure_function = dist_format_with_location_signed_dist_discrete,
                 plot_title = plot_title,parameters = w)
print(RDM_model_format_with_location_signed_dist_k0.5$plot)
#w=0.8-------------------------------------------------------
w=0.8 # The relative weight of distance difference compared to format difference
plot_title=paste0("Format with location information \n (signed dist.; weight(dist) = ",w,")")
RDM_model_format_with_location_signed_dist_k0.8=
  RDM_model_plot(df = df_trials_discrete_36,
                 dis_measure_function = dist_format_with_location_signed_dist_discrete,
                 plot_title = plot_title,parameters = w)
print(RDM_model_format_with_location_signed_dist_k0.8$plot)

#Model-8: Format without location information(with signed distance information)-------------------------------------
#w is The relative weight of distance difference compared to format difference

#w=0.2-------------------------------------------------------
w=0.2 # The relative weight of distance difference compared to format difference
plot_title=paste0("Format without location information \n (signed dist.; weight(dist) = ",w,")")
RDM_model_format_without_location_signed_dist_k0.2=
  RDM_model_plot(df = df_trials_discrete_36,
                 dis_measure_function = dist_format_without_location_signed_dist_discrete,
                 plot_title = plot_title,parameters = w)
print(RDM_model_format_without_location_signed_dist_k0.2$plot)
#w=0.5-------------------------------------------------------
w=0.5 # The relative weight of distance difference compared to format difference
plot_title=paste0("Format without location information \n (signed dist.; weight(dist) = ",w,")")
RDM_model_format_without_location_signed_dist_k0.5=
  RDM_model_plot(df = df_trials_discrete_36,
                 dis_measure_function = dist_format_without_location_signed_dist_discrete,
                 plot_title = plot_title,parameters = w)
print(RDM_model_format_without_location_signed_dist_k0.5$plot)
#w=0.8-------------------------------------------------------
w=0.8 # The relative weight of distance difference compared to format difference
plot_title=paste0("Format without location information \n (signed dist.; weight(dist) = ",w,")")
RDM_model_format_without_location_signed_dist_k0.8=
  RDM_model_plot(df = df_trials_discrete_36,
                 dis_measure_function = dist_format_without_location_signed_dist_discrete,
                 plot_title = plot_title,parameters = w)
print(RDM_model_format_without_location_signed_dist_k0.8$plot)

# RSA ID (unique ID: 1~18)==========================================)======================================================
#Model: Null model (not encoding signals)-------------------------------------
plot_title="Null model \n (not encoding any signal)"
RDM =
  RDM_model_plot(df = df_trials_discrete_18,
                 dis_measure_function = dist_null,
                 plot_title = plot_title)
print(RDM$plot)
#Model: FCF: All-the-same model (not encoding signals)-------------------------------------
plot_title="Full cross-format similarity \n (extreme abstraction)"
RDM =
  RDM_model_plot(df = df_trials_discrete_18,
                 dis_measure_function = dist_all_the_same,
                 plot_title = plot_title)
print(RDM$plot)

#Model: NCF: Pair with the same format (no distance information)-------------------------------------
plot_title="No cross-format similarity \n (no distance information)"
RDM =
  RDM_model_plot(df = df_trials_discrete_18,
                 dis_measure_function = dist_same_format_pair_no_distance,
                 plot_title = plot_title)
print(RDM$plot)
#Model: GCF: Format without location information(no distance information)-------------------------------------
plot_title="Graded cross-format similarity \n (no distance information)"
RDM_model2_format_without_location =
  RDM_model_plot(df = df_trials_discrete_18,
                 dis_measure_function = dist_format_without_location_no_distance,
                 plot_title = plot_title )
print(RDM_model2_format_without_location$plot)

# Discrete distance==========================================)======================================================
#Model: absolute distance information (no format or position)-------------------------------------
plot_title="Absolute distance \n Mag. = |v1-v2| \n Dis. = |Mag1-Mag2| \n (No format or position)"
RDM=
  RDM_model_plot(df = df_trials_discrete_18,
                 dis_measure_function = dist_only_abs_dist_discrete,
                 plot_title = plot_title)
print(RDM$plot)
#Model: signed distance information (no format or position)-------------------------------------
plot_title="Signed distance \n Mag. = (v1-v2) \n Dis. = |Mag1-Mag2| \n (No format or position)"
RDM=
  RDM_model_plot(df = df_trials_discrete_18,
                 dis_measure_function = dist_only_signed_dist_discrete,
                 plot_title = plot_title)
print(RDM$plot)

#Model: Difficulty - Accuracy--------------------------
age_group_name=c("adult","2_grade","5_grade","all")
age_group_plot_name=c("Adult","2 grade","5 grade","All")
g_list_ACC=c()
for (age_group_index in 1:length(age_group_name)){
  plot_title=paste0("Difficulty - Accuracy \n (",age_group_plot_name[age_group_index],")")
  parameters=list(age_group=age_group_name[age_group_index])
  RDM=
    RDM_model_plot(df = df_trials_discrete_18,
                   dis_measure_function = dist_diff_acc,
                   parameters = parameters,
                   plot_title = plot_title)
  # print(RDM$plot)
  g_list_ACC[[age_group_index]]=RDM$plot
}
#Model: Difficulty - RT--------------------------
age_group_name=c("adult","2_grade","5_grade","all")
age_group_plot_name=c("Adult","2 grade","5 grade","All")
g_list_RT=c()
for (age_group_index in 1:length(age_group_name)){
  plot_title=paste0("Difficulty - RT \n (",age_group_plot_name[age_group_index],")")
  parameters$age_group=age_group_name[age_group_index]
  RDM=
    RDM_model_plot(df = df_trials_discrete_18,
                   dis_measure_function = dist_diff_RT,
                   parameters = parameters,
                   plot_title = plot_title)
  # print(RDM$plot)
  g_list_RT[[age_group_index]]=RDM$plot
}
# library(gridExtra)
# + theme(legend.position="none")
gridExtra::grid.arrange(grobs=c(g_list_ACC,
                                g_list_RT),
                        ncol=4)
# Hybrid models-------------------)------------------------------------------------
#Model: GCF: Format without location information(with absolute distance information)-------------------------------------
for (w_index in 1:length(LIST_W)){
  w=LIST_W[w_index] # The relative weight of distance difference compared to format difference
  plot_title=paste0("Graded cross-format similarity \n (absolute dist.; weight(dist) = ",w,")")
  RDM=
    RDM_model_plot(df = df_trials_discrete_18,
                   dis_measure_function = dist_format_without_location_abs_dist_discrete,
                   plot_title = plot_title,
                   parameters = w)
  print(RDM$plot)
}
#Model: GCF: Format without location information(with signed distance information)-------------------------------------
#w is The relative weight of distance difference compared to format difference
for (w_index in 1:length(LIST_W)){
  w=LIST_W[w_index] # The relative weight of distance difference compared to format difference
  plot_title=paste0("Graded cross-format similarity \n (signed dist.; weight(dist) = ",w,")")
  RDM=
    RDM_model_plot(df = df_trials_discrete_18,
                   dis_measure_function = dist_format_without_location_signed_dist_discrete,
                   plot_title = plot_title,
                   parameters = w)
  print(RDM$plot)
}
#Model: NCF: Pair with the same format (with absolute distance information)-------------------------------------
#w is The relative weight of distance difference compared to format difference
for (w_index in 1:length(LIST_W)){
  w=LIST_W[w_index] # The relative weight of distance difference compared to format difference
  plot_title=paste0("No cross-format similarity \n (absolute dist.; weight(dist) = ",w,")")
  RDM=
    RDM_model_plot(df = df_trials_discrete_18,
                   dis_measure_function = dist_same_format_pair_abs_distance,
                   plot_title = plot_title,
                   parameters = w)
  print(RDM$plot)
}
#Model: NCF: Pair with the same format (with signed distance information)-------------------------------------
#w is The relative weight of distance difference compared to format difference
for (w_index in 1:length(LIST_W)){
  w=LIST_W[w_index] # The relative weight of distance difference compared to format difference
  plot_title=paste0("No cross-format similarity \n (signed dist.; weight(dist) = ",w,")")
  RDM=
    RDM_model_plot(df = df_trials_discrete_18,
                   dis_measure_function = dist_same_format_pair_signed_distance,
                   plot_title = plot_title,
                   parameters = w)
  print(RDM$plot)
}

#Model: GCF & Difficulty (Accuracy)---------------------------
#w is The relative weight of distance difference compared to format difference
age_group_name=c("adult","2_grade","5_grade","all")
age_group_plot_name=c("Adult","2 grade","5 grade","All")
plot_index=1
g_list=c()
for (age_group_index in 1:length(age_group_name)){
  for (w_index in 1:length(LIST_W)){
    w=LIST_W[w_index] # The relative weight of 2nd distance(i.e.,difficulty) compared to the first difference
    age_group=age_group_name[age_group_index]
    parameters=list(w2=w,age_group=age_group)

    plot_title=paste0("Graded cross-format similarity \n",
                      "Diff (ACC: ",age_group_plot_name[age_group_index],") ; weight = ",w)
    RDM=
      RDM_model_plot(df = df_trials_discrete_18,
                     dis_measure_function = dist_format_without_location_diff_acc,
                     plot_title = plot_title,
                     parameters = parameters)
    # g_list[[plot_index]]=RDM$plot
    # print(RDM$plot)
    ggsave(plot = RDM$plot,filename = paste0("/study3/devfracs/DCM_YunShiuan/RSA/",plot_index,".png"),
           width = 3.5,height = 3)
    plot_index=plot_index+1
  }
}
#Model: NCF & Difficulty (Accuracy)---------------------------
#w is The relative weight of distance difference compared to format difference
age_group_name=c("adult","2_grade","5_grade","all")
age_group_plot_name=c("Adult","2 grade","5 grade","All")
plot_index=1
g_list=c()
for (age_group_index in 1:length(age_group_name)){
  for (w_index in 1:length(LIST_W)){
    w=LIST_W[w_index] # The relative weight of 2nd distance(i.e.,difficulty) compared to the first difference
    age_group=age_group_name[age_group_index]
    parameters=list(w2=w,age_group=age_group)

    plot_title=paste0("No cross-format similarity \n",
                      "Diff (ACC: ",age_group_plot_name[age_group_index],") ; weight = ",w)
    RDM=
      RDM_model_plot(df = df_trials_discrete_18,
                     dis_measure_function = dist_same_format_pair_diff_acc,
                     plot_title = plot_title,
                     parameters = parameters)
    # g_list[[plot_index]]=RDM$plot
    # print(RDM$plot)
    ggsave(plot = RDM$plot,filename = paste0("/study3/devfracs/DCM_YunShiuan/RSA/",plot_index,".png"),
           width = 3.5,height = 3)
    plot_index=plot_index+1
  }
}

# Partial Hybrid models----------------------)---------------------
#Model: GCF: Graded cross-format similarity weighted by non-symbolic absolute dist---------------------------
for (w_index in 1:length(LIST_W)){
  w=LIST_W[w_index] # The relative weight of distance difference compared to format difference
  plot_title=paste0("Graded cross-format similarity \n (non-symbolic abs. dist: w = ",w,")")
  RDM=
    RDM_model_plot(df = df_trials_discrete_18,
                   dis_measure_function = dist_format_without_location_LL_abs_dist_discrete,
                   plot_title = plot_title,
                   parameters = w)
  print(RDM$plot)
}
#Model: GCF: Graded cross-format similarity  weighted by non-symbolic signed dist---------------------------
for (w_index in 1:length(LIST_W)){
  w=LIST_W[w_index] # The relative weight of distance difference compared to format difference
  plot_title=paste0("Graded cross-format similarity \n (non-symbolic signed dist: w = ",w,")")
  RDM=
    RDM_model_plot(df = df_trials_discrete_18,
                   dis_measure_function = dist_format_without_location_LL_signed_dist_discrete,
                   plot_title = plot_title,
                   parameters = w)
  print(RDM$plot)
}
#Model: NCF: No cross-format similarity  weighted by non-symbolic absolute dist---------------------------
for (w_index in 1:length(LIST_W)){
  w=LIST_W[w_index] # The relative weight of distance difference compared to format difference
  plot_title=paste0("No cross-format similarity \n (non-symbolic abs. dist: w = ",w,")")
  RDM=
    RDM_model_plot(df = df_trials_discrete_18,
                   dis_measure_function = dist_same_format_pair_LL_abs_dist_discrete,
                   plot_title = plot_title,
                   parameters = w)
  print(RDM$plot)
}
#Model: NCF: No cross-format similarity  weighted by non-symbolic signed dist---------------------------
for (w_index in 1:length(LIST_W)){
  w=LIST_W[w_index] # The relative weight of distance difference compared to format difference
  plot_title=paste0("No cross-format similarity \n (non-symbolic signed dist: w = ",w,")")
  RDM=
    RDM_model_plot(df = df_trials_discrete_18,
                   dis_measure_function = dist_same_format_pair_LL_signed_dist_discrete,
                   plot_title = plot_title,
                   parameters = w)
  print(RDM$plot)
}


# Log-transformed Discrete distance ============================================================)======================================================
#Model: absolute distance information (no format or position)-------------------------------------
list_parameters=list(parameters=list(log_version="v1"),
                     parameters=list(log_version="v2"),
                     parameters=list(log_version="v3"))
list_plot_names=c("Absolute log-magnitude (ver.1) \n Mag. = |log(v1)-log(v2)| \n Dis. = |Mag1-Mag2| \n (No format or position)",
                  "Absolute log-magnitude (ver.2) \n Mag. = log( |v1-v2| ) \n Dis. = |Mag1-Mag2| \n (No format or position)",
                  "Absolute log-magnitude (ver.3) \n Mag. = log( |log(v1)-log(v2)| ) \n Dis. = |Mag1-Mag2| \n (No format or position)")
g_list_dist=c()
for (plot_index in 1:length(list_parameters)){
  plot_title=list_plot_names[[plot_index]]
  parameters=list_parameters[[plot_index]]
  RDM=
    RDM_model_plot(df = df_trials_discrete_18,
                   dis_measure_function = dist_only_abs_dist_discrete_log,
                   plot_title = plot_title,
                   parameters = parameters)
  # print(RDM$plot)
  g_list_dist[[plot_index]]=RDM$plot
}
#Add in the raw distance
plot_title="Absolute distance \n Mag. = |v1-v2| \n Dis. = |Mag1-Mag2| \n (No format or position)"
RDM_abs_dist=
  RDM_model_plot(df = df_trials_discrete_18,
                 dis_measure_function = dist_only_abs_dist_discrete,
                 plot_title = plot_title)
gridExtra::grid.arrange(grobs=c(list(RDM_abs_dist$plot),g_list_dist),
                        ncol=2)
#Model: signed distance information (no format or position)-------------------------------------

list_parameters=list(parameters=list(log_version="v1"),
                     parameters=list(log_version="v2"))
list_plot_names=c("Signed log-magnitude (ver.1) \n Mag. = log(v1)-log(v2) \n Dis. = |Mag1-Mag2| \n (No format or position)",
                  "Signed log-magnitude (ver.2) \n Mag. = sign*shift(log(|v1-v2|)) \n Dis. = |Mag1-Mag2| \n (No format or position)")
g_list_dist=c()
for (plot_index in 1:length(list_parameters)){
  plot_title=list_plot_names[[plot_index]]
  parameters=list_parameters[[plot_index]]
  RDM=
    RDM_model_plot(df = df_trials_discrete_18,
                   dis_measure_function = dist_only_signed_dist_discrete_log,
                   plot_title = plot_title,
                   parameters = parameters)
  # print(RDM$plot)
  g_list_dist[[plot_index]]=RDM$plot
}
#Add in the raw distance
plot_title="Signed distance \n Mag. = (v1-v2) \n Dis. = |Mag1-Mag2| \n (No format or position)"
RDM_abs_dist=
  RDM_model_plot(df = df_trials_discrete_18,
                 dis_measure_function = dist_only_signed_dist_discrete,
                 plot_title = plot_title)
gridExtra::grid.arrange(grobs=c(list(RDM_abs_dist$plot),g_list_dist),
                        ncol=2)
# Hybrid models-------------------)------------------------------------------------
#Model: GCF: Format without location information(with absolute distance information)------------------
#Model: NCF: Pair with the same format (with absolute distance information)---------------------
list_parameters=list(parameters=list(log_version="v1"),
                     parameters=list(log_version="v2"),
                     parameters=list(log_version="v3"))
list_plot_names_log=c("Absolute log-magnitude (ver.1) \n Mag. = |log(v1) - log(v2)|",
                  "Absolute log-magnitude (ver.2) \n Mag. = log( |v1-v2| ) ",
                  "Absolute log-magnitude (ver.3) \n Mag. = log( |log(v1)-log(v2)| )")
list_plot_names_format=c("Graded Cross-Format (GCF)","No Cross-Format (NCF)")
list_dist_log_hybrid_functions=c(dist_format_without_location_abs_dist_discrete_log,
                                 dist_same_format_pair_abs_distance_log)
#w is The relative weight of distance difference compared to format difference
list_ggplot_mixed_mag_log=c()
plot_index=1
for (format_index in 1: length(list_plot_names_format)){
    plot_title_format=list_plot_names_format[[format_index]]
    dist_function=list_dist_log_hybrid_functions[[format_index]]
    
  for(log_version_index in 1:length(list_parameters)){
    parameters=list_parameters[[log_version_index]]
    plot_title_log_version=list_plot_names_log[[log_version_index]]
    
    for (w_index in 1:length(LIST_W)){
      w=LIST_W[w_index] # The relative weight of distance difference compared to format difference
      parameters$w2=w
      plot_title=paste0(plot_title_format,"\n",plot_title_log_version,"\n weight (mag.) = ",w)
      RDM=
        RDM_model_plot(df = df_trials_discrete_18,
                       dis_measure_function = dist_function,
                       plot_title = plot_title,
                       parameters = parameters)
      list_ggplot_mixed_mag_log[[plot_index]]=RDM$plot
      plot_index=plot_index+1
      # print(RDM$plot)
      print(paste0("format:",format_index,
                   "log-version",log_version_index,
                   "w",w_index))
    }
  }
}
gridExtra::grid.arrange(grobs=list_ggplot_mixed_mag_log[1:12],
                        ncol=4)
gridExtra::grid.arrange(grobs=list_ggplot_mixed_mag_log[13:24],
                        ncol=4)

#Model: GCF: Format without location information(with signed distance information)------------------
#Model: NCF: Pair with the same format (with signed distance information)---------------------
list_parameters=list(parameters=list(log_version="v1"),
                     parameters=list(log_version="v2"))
list_plot_names_log=c("Signed log-magnitude (ver.1) \n Mag. = log(v1)-log(v2) \n Dis. = |Mag1-Mag2| \n (No format or position)",
                  "Signed log-magnitude (ver.2) \n Mag. = sign*shift(log(|v1-v2|)) \n Dis. = |Mag1-Mag2| \n (No format or position)")
list_plot_names_format=c("Graded Cross-Format (GCF)","No Cross-Format (NCF)")
list_dist_log_hybrid_functions=c(dist_format_without_location_signed_dist_discrete_log,
                                 dist_same_format_pair_LL_signed_dist_discrete_log)
#w is The relative weight of distance difference compared to format difference
list_ggplot_mixed_mag_log=c()
plot_index=1
for (format_index in 1: length(list_plot_names_format)){
  plot_title_format=list_plot_names_format[[format_index]]
  dist_function=list_dist_log_hybrid_functions[[format_index]]
  
  for(log_version_index in 1:length(list_parameters)){
    parameters=list_parameters[[log_version_index]]
    plot_title_log_version=list_plot_names_log[[log_version_index]]
    
    for (w_index in 1:length(LIST_W)){
      w=LIST_W[w_index] # The relative weight of distance difference compared to format difference
      parameters$w2=w
      plot_title=paste0(plot_title_format,"\n",plot_title_log_version,"\n weight (mag.) = ",w)
      RDM=
        RDM_model_plot(df = df_trials_discrete_18,
                       dis_measure_function = dist_function,
                       plot_title = plot_title,
                       parameters = parameters)
      list_ggplot_mixed_mag_log[[plot_index]]=RDM$plot
      plot_index=plot_index+1
      # print(RDM$plot)
      print(paste0("format:",format_index,
                   "log-version",log_version_index,
                   "w",w_index))
    }
  }
}
gridExtra::grid.arrange(grobs=list_ggplot_mixed_mag_log[1:(length(list_parameters)*length(LIST_W))],
                        ncol=3)
gridExtra::grid.arrange(grobs=list_ggplot_mixed_mag_log[((length(list_parameters)*length(LIST_W))+1):
                                                          length(list_ggplot_mixed_mag_log)],
                        ncol=3)

# Partial Hybrid models----------------------)---------------------
#Model: GCF: Graded cross-format similarity weighted by non-symbolic absolute dist---------------------------
#Model: NCF: No cross-format similarity  weighted by non-symbolic absolute dist---------------------------
list_parameters=list(parameters=list(log_version="v1"),
                     parameters=list(log_version="v2"),
                     parameters=list(log_version="v3"))
list_plot_names_log=c("Absolute log-magnitude (ver.1) \n Mag. = |log(v1) - log(v2)|",
                      "Absolute log-magnitude (ver.2) \n Mag. = log( |v1-v2| ) ",
                      "Absolute log-magnitude (ver.3) \n Mag. = log( |log(v1)-log(v2)| )")
list_plot_names_format=c("Graded Cross-Format (GCF)","No Cross-Format (NCF)")
list_dist_log_hybrid_functions=c(dist_format_without_location_LL_abs_dist_discrete_log,
                                 dist_same_format_pair_abs_distance_log)
#w is The relative weight of distance difference compared to format difference
list_ggplot_mixed_mag_log=c()
plot_index=1
for (format_index in 1: length(list_plot_names_format)){
  plot_title_format=list_plot_names_format[[format_index]]
  dist_function=list_dist_log_hybrid_functions[[format_index]]
  
  for(log_version_index in 1:length(list_parameters)){
    parameters=list_parameters[[log_version_index]]
    plot_title_log_version=list_plot_names_log[[log_version_index]]
    
    for (w_index in 1:length(LIST_W)){
      w=LIST_W[w_index] # The relative weight of distance difference compared to format difference
      parameters$w2=w
      plot_title=paste0(plot_title_format,"\n",plot_title_log_version,"\n weighted in the LL batch (mag.) = ",w)
      RDM=
        RDM_model_plot(df = df_trials_discrete_18,
                       dis_measure_function = dist_function,
                       plot_title = plot_title,
                       parameters = parameters)
      list_ggplot_mixed_mag_log[[plot_index]]=RDM$plot
      plot_index=plot_index+1
      # print(RDM$plot)
      print(paste0("format:",format_index,
                   "log-version",log_version_index,
                   "w",w_index))
    }
  }
}
gridExtra::grid.arrange(grobs=list_ggplot_mixed_mag_log[1:9],
                        ncol=3)
gridExtra::grid.arrange(grobs=list_ggplot_mixed_mag_log[10:18],
                        ncol=3)
#Caution:----------------------------------)------------------------------------------
#The inevitable noise within the CN resulted from collapsing trials with different location -----------------------------------------------
# The inevitable noise within the CN resulted from collapsing trials with different location  (FL & LF => CN).
# 18 IDs: 3 format x 3 distance levels x 2 direction.
# The location information could not be not encoded when each ID has 2 repetition.
#Model: Format with location information(no distance information)-------------------------------------
plot_title="Format with location information \n (no distance information)"
RDM_model1_format_with_location =
  RDM_model_plot(df = df_trials_discrete_18,
                 dis_measure_function = dist_format_with_location_no_distance,
                 plot_title = plot_title )
print(RDM_model1_format_with_location$plot)
#Model: Format with location information(with absolute distance information)-------------------------------------
#w is The relative weight of distance difference compared to format difference
#w=0.2-------------------------------------------------------
w=0.2 # The relative weight of distance difference compared to format difference
plot_title = paste0("Format with location information \n (absolute dist.; weight(dist) = ",w,")")
RDM_model3_format_with_location_abs_dist_k0.2=
  RDM_model_plot(df = df_trials_discrete_18,
                 dis_measure_function = dist_format_with_location_abs_dist_discrete,
                 plot_title = plot_title,
                 parameters = w )
print(RDM_model3_format_with_location_abs_dist_k0.2$plot)
#w=0.5-------------------------------------------------------
w=0.5 # The relative weight of distance difference compared to format difference
plot_title = paste0("Format with location information \n (absolute dist.; weight(dist) = ",w,")")
RDM_model3_format_with_location_abs_dist_k0.5=
  RDM_model_plot(df = df_trials_discrete_18,
                 dis_measure_function = dist_format_with_location_abs_dist_discrete,
                 plot_title = plot_title,
                 parameters = w )
print(RDM_model3_format_with_location_abs_dist_k0.5$plot)
#w=0.8-------------------------------------------------------
w=0.8 # The relative weight of distance difference compared to format difference
plot_title = paste0("Format with location information \n (absolute dist.; weight(dist) = ",w,")")
RDM_model3_format_with_location_abs_dist_k0.8=
  RDM_model_plot(df = df_trials_discrete_18,
                 dis_measure_function = dist_format_with_location_abs_dist_discrete,
                 plot_title = plot_title,
                 parameters = w )
print(RDM_model3_format_with_location_abs_dist_k0.8$plot)
#Model: Format with location information(with signed distance information)-------------------------------------
#w is The relative weight of distance difference compared to format difference
#w=0.2-------------------------------------------------------
w=0.2 # The relative weight of distance difference compared to format difference
plot_title=paste0("Format with location information \n (signed dist.; weight(dist) = ",w,")")
RDM_model_format_with_location_signed_dist_k0.2=
  RDM_model_plot(df = df_trials_discrete_18,
                 dis_measure_function = dist_format_with_location_signed_dist_discrete,
                 plot_title = plot_title,parameters = w)
print(RDM_model_format_with_location_signed_dist_k0.2$plot)
#w=0.5-------------------------------------------------------
w=0.5 # The relative weight of distance difference compared to format difference
plot_title=paste0("Format with location information \n (signed dist.; weight(dist) = ",w,")")
RDM_model_format_with_location_signed_dist_k0.5=
  RDM_model_plot(df = df_trials_discrete_18,
                 dis_measure_function = dist_format_with_location_signed_dist_discrete,
                 plot_title = plot_title,parameters = w)
print(RDM_model_format_with_location_signed_dist_k0.5$plot)
#w=0.8-------------------------------------------------------
w=0.8 # The relative weight of distance difference compared to format difference
plot_title=paste0("Format with location information \n (signed dist.; weight(dist) = ",w,")")
RDM_model_format_with_location_signed_dist_k0.8=
  RDM_model_plot(df = df_trials_discrete_18,
                 dis_measure_function = dist_format_with_location_signed_dist_discrete,
                 plot_title = plot_title,parameters = w)
print(RDM_model_format_with_location_signed_dist_k0.8$plot)
