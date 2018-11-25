%% Part 3 - ii: RSA : Prepare data for RSA usage (before construct RDM)
% trial 4 in the README.txt.
% Part3: 2nd order analysis.
% ii: use t maps (36 trials with discrete distance) 
% instead of beta maps (144  unique trials with continuous distance)
% [Run info has already been averaged out. Therefore, the number of runs no longer matters.]
addpath(genpath('D:\Yun-Shiuan_LAMBDA\rsatoolbox\'));

%% Set up
% Constants
% path
path_1st_result='D:\Yun-Shiuan_LAMBDA\RSA\Part2_1st_oder_analysis_RSA_ID_discrete_unique_36';
path_2nd_result='D:\Yun-Shiuan_LAMBDA\RSA\Part3_2nd_order_analysis_RSA_ID_discrete_unique_36';
version_num='4';
% file
file_RDM_result='D:\Yun-Shiuan_LAMBDA\RSA\Part1_neural_and_conceptual_RDMs_RSA_ID_discrete_repeated_36\RDMs\RDMs_and_options_neural_and_conceptual_trial_4_06-Mar-2018';
file_1st_user_option=fullfile(path_1st_result,'user_options_visualize_RDM_4_06-Mar-2018.mat');
% Parameters
list_VOI_name={"R_V1","L_V1",...
               "R_IPS","L_IPS"}';
amount_each_format=6;
% Load in needed data
load(file_RDM_result,...
     'neural_RDMs','neural_RDMs_avg_subjects','conceptual_RDMs',...
     'user_options_RDM_construction')
load(file_1st_user_option,'user_options_visualize_RDM');
% Options for2nd order analysis
clear user_options_2nd_order
% Edit the color bar color due to Matlab version issue
%(R2015b change default color maps. Should now be colormap(jet(256)))
% user_options_visualize_RDM.colourScheme=colormap(jet(256)); [Already edify the source codes]
user_options_2nd_order=user_options_visualize_RDM;

user_options_2nd_order.analysisName = ['RSA_trial_',version_num,'_',date];
user_options_2nd_order.rootPath=path_2nd_result;                   
user_options_2nd_order.RDMcorrelationType='Kendall_taua';
user_options_2nd_order.distanceMeasure='Kendall_taua';
reduced_conceptual_RDMs=conceptual_RDMs(1,[3,5,7,10,13,15,17,22,24,26,29,32,34,36,39,40]);


%% (Step 1) Visualizatopm: 2nd order simmilarity matrix and 2nd order MDS
% Display a 2nd order simmilarity matrix between neural RDMs and conceptual
% RDMs
pairwiseCorrelateRDMs({neural_RDMs_avg_subjects, conceptual_RDMs},...
                       user_options_2nd_order);

% MDS: Plot all RDMs on a MDS plot to visualise pairwise
% distances.-----------------------------------------------------------
% Reduce the conceptual RDMs to be plotted.
% Without reducing: Points in the configuration have co-located.  Try a different
% starting point, or use a different criterion.
MDSRDMs({neural_RDMs_avg_subjects, reduced_conceptual_RDMs}, user_options_2nd_order,...
    struct('titleString', 'MDS of neural and conceptual'));
% If one want to save manually instead of using the automatic saving function of the
% MDSRDMs()
% saveas(gcf,fullfile(path_2nd_result,[user_options_2nd_order.analysisName,'_SecondOrderMDS_.pdf']))


%% (Step 2)stats_p_r = compareRefRDM2candRDMs(refRDM, candRDMs[, userOptions]) - statistical inference:
clear models
for conceptual_index=1:numel(reduced_conceptual_RDMs)
    models{conceptual_index}=reduced_conceptual_RDMs(conceptual_index);
end
for neural_index=1:numel(neural_RDMs_avg_subjects)
    models{conceptual_index+neural_index}=neural_RDMs_avg_subjects(neural_index);
end
% test the relatedness and compare the candidate RDMs
user_options_2nd_order.RDMcorrelationType='Kendall_taua';
user_options_2nd_order.RDMrelatednessTest = 'subjectRFXsignedRank';
user_options_2nd_order.RDMrelatednessThreshold = 0.05;
user_options_2nd_order.figureIndex = [1:2];
user_options_2nd_order.RDMrelatednessMultipleTesting = 'FDR';
user_options_2nd_order.candRDMdifferencesTest = 'subjectRFXsignedRank';
user_options_2nd_order.candRDMdifferencesThreshold = 0.05;
user_options_2nd_order.candRDMdifferencesMultipleTesting = 'none';
%Remain to be
%debugged: Undefined function 'mod' for input arguments of type 'matlab.ui.Figure'.----------------------------------------------------------------
% R_V1
stats_p_r=compareRefRDM2candRDMs(neural_RDMs(1,:), models, user_options_2nd_order);
% L_V1
stats_p_r=compareRefRDM2candRDMs(neural_RDMs(2,:), models, user_options_2nd_order);
% R IPS
stats_p_r=compareRefRDM2candRDMs(neural_RDMs(3,:), models, user_options_2nd_order);
% L IPS
stats_p_r=compareRefRDM2candRDMs(neural_RDMs(4,:), models, user_options_2nd_order);