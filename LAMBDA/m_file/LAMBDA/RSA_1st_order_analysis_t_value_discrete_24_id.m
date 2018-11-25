%% Part 2 - ii: RSA : 1st order analysis 
% trial 1~3 in the README.txt.
% Part2: 1st order analysis (visualize RDM, MDS, clustering)
% ii: use t maps (24 unique trials with discrete distance) 
% instead of beta maps (144  unique trials with continuous distance)
% [Run info has already been averaged out. Therefore, the number of runs no longer matters.]
%% Set up
% Constants
% path
path_RDM_result='D:\Yun-Shiuan_LAMBDA\RSA\Part1_neural_and_conceptual_RDMs';
file_RDM_result=fullfile(path_RDM_result,'RDMs','RDMs_and_options_neural_and_conceptual_trial_3_02-Mar-2018.mat');
path_1st_result='D:\Yun-Shiuan_LAMBDA\RSA\Part2_1st_oder_analysis';
% Parameters
list_VOI_name={"R_V1","L_V1",...
               "R_IPS","L_IPS"}';
amount_each_format=6;
condition_labels=cellstr(string(1:24));
version_num='3';
% Load in needed data
load(file_RDM_result,...
     'neural_RDMs','neural_RDMs_avg_subjects','conceptual_RDMs',...
     'user_options_RDM_construction')
% Options for visualization
clear user_options_visualize_RDM
% Edit the color bar color due to Matlab version issue
%(R2015b change default color maps. Should now be colormap(jet(256)))
% user_options_visualize_RDM.colourScheme=colormap(jet(256)); [Already edify the source codes]
user_options_visualize_RDM=user_options_RDM_construction;
user_options_visualize_RDM.saveFiguresPDF=true;
user_options_visualize_RDM.saveFigurePDF=true;
user_options_visualize_RDM.tightInset=true;
user_options_visualize_RDM.analysisName = ['RSA_trial_',version_num,'_',date];
user_options_visualize_RDM.conditionColours=[repmat([1 0 0],amount_each_format,1);...
                                             repmat([0 1 0],amount_each_format,1);...
                                             repmat([1 0.7 0],amount_each_format,1);...
                                             repmat([0 0 1],amount_each_format,1)];
user_options_visualize_RDM.conditionLabels=cellstr(string(1:24));
                       
%% (Step 1) Visulize the Neural and Conceptual RDMs: Each subject and average across subjects
% Save the RDM matrix plot into the RDMs folder instead of the 1st order
% one
user_options_visualize_RDM.rootPath=path_RDM_result;

% Plot Neural RDMs--------------------------------
% Avg subjects
figureRDMs(concatenateRDMs(neural_RDMs_avg_subjects), user_options_visualize_RDM,...
    struct('fileName', 'neural_RDMs_avg_subjects', 'figureNumber', 1));
% for each VOI
for v=1:size(list_VOI_name,1)
    figureRDMs(concatenateRDMs(neural_RDMs(v,:)), user_options_visualize_RDM,...
        struct('fileName', ['neural_RDMs_',char(list_VOI_name{v})], 'figureNumber', 1));
end
% Plot Conceptual RDMs--------------------------------
figureRDMs(concatenateRDMs(conceptual_RDMs), user_options_visualize_RDM,...
    struct('fileName', 'conceptual_RDMs', 'figureNumber', 1));

%% (Step 2) Determine dendrograms for the clustering of the conditions
%%%%%%%%%%%%%WIERD! The condition color and condition numbers do not match!
user_options_visualize_RDM.rootPath=path_1st_result; %Save the results starting from here to the 1st order part
dendrogramConditions(neural_RDMs_avg_subjects, user_options_visualize_RDM,...
    struct(...
    'useAlternativeConditionLabels', true,...
    'alternativeConditionLabels', {condition_labels}));

%% (Step 3) Determine dendrograms for the clustering of the conditions
MDSConditions(neural_RDMs_avg_subjects, user_options_visualize_RDM,...
    struct('alternativeConditionLabels', {condition_labels}));

file_save=...
    fullfile(user_options_visualize_RDM.rootPath,['user_options_visualize_RDM_',version_num,'_',date,'.mat']);

save(file_save,'user_options_visualize_RDM')
