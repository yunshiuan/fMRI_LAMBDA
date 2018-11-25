%% Part 2 - iv-b: RSA : 1st order analysis 
% trial 6 in the README.txt.
% Part2: 1st order analysis (visualize RDM, MDS, clustering)
%
% iv: RDA_ID_discrete 18 ID:
% - RSA ID : 1~18 (3 formats x 3 distance levels x 2 directions)
% - The motivation is to increase the statistical efficiency (Conditions in the Trial 4 is too noisy),
% while still keep that same-notation trials have the same statistical efficiency
% as the cross-notation trials (as Trial 4 did).
% [LL: 12 trials = 3 distance levels x 2 direction x 2 repetition]
% [FF: 12 trials = 3 distance levels x 2 direction x 2 repetition]
% [cross-notation: 12 trials = 3 distance levels x 2 direction x 2 [repetition]
% (Having more repetition in a condition entails a smaller degree of random noise.)

% The motivation for having 18 trials is to ensure beta estimate of each condition is subjected to same amount of noise.
% (Having more repetition in a condition entails a smaller degree of random noise.)

% [Run info has already been averaged out. Therefore, the number of runs no longer matters.]
% b: Split neural RDM by age: 2nd grader and 5th grader. 
% Run 1st and 2nd order RSA on each age group separately.
%% Set up
% Constants
% path
path_RDM_result='D:\Yun-Shiuan_LAMBDA\RSA\trial_6_id_18_split_age\Part1_neural_and_conceptual_RDMs';
file_RDM_result=fullfile(path_RDM_result,'RDMs','RDMs_and_options_neural_and_conceptual_trial_6_09-Mar-2018.mat');
path_1st_result='D:\Yun-Shiuan_LAMBDA\RSA\trial_6_id_18_split_age\Part2_1st_oder_analysis';
%files
file_valid_run='D:\Yun-Shiuan_LAMBDA\Run_inclusion_info\inclusive_runs_indexes.csv';
file_demographic='D:\Yun-Shiuan_LAMBDA\demographic\tidy_demographic.csv';
% Parameters
list_VOI_name={"R_V1","L_V1",...
               "R_IPS","L_IPS"}';
list_group_name={'all_sub','2_grade','5_grade'};
number_each_format=6;
num_trials=18; % 12 trials per format but only 6 are unique
condition_labels=cellstr(string(1:num_trials));
version_num='6';
% Load in needed data
load(file_RDM_result,...
     'neural_RDMs',...
     'neural_RDMs_avg_2_grade','neural_RDMs_avg_5_grade','neural_RDMs_avg_all',...
     'conceptual_RDMs',...
     'user_options_RDM_construction')
 list_neural_RDMs={neural_RDMs_avg_all,neural_RDMs_avg_2_grade,neural_RDMs_avg_5_grade};
% Read in run inclusion index info (note that only subject matters for
% t-map based analysis)
run_inclusion_index=cellfun(@(x) regexprep(x,'"',''),...
    read_mixed_csv(file_valid_run,','),'un',0);
run_inclusion_index=table(run_inclusion_index(2:end,2),run_inclusion_index(2:end,3),...
    'VariableNames',{'sub_id','run_num'});
% Derive subjects with valid runs
subject_list=unique(run_inclusion_index.sub_id);

% Read in demographic data
table_demographic=read_mixed_csv_to_table(file_demographic);
table_demographic=table_demographic(:,{'sub_id','grade'});
table_demographic.sub_id=lower(table_demographic.sub_id);
table_demographic=join(table(subject_list,'VariableNames',{'sub_id'}),...
    table_demographic);

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
user_options_visualize_RDM.conditionColours=[repmat([1 0 0],number_each_format,1);...
                                             repmat([0 1 0],number_each_format,1);...
                                             repmat([0 0 1],number_each_format,1)];
user_options_visualize_RDM.conditionLabels=cellstr(string(1:num_trials));
user_options_visualize_RDM.rootPath=path_RDM_result;

%% (Step 1) Visulize the Neural and Conceptual RDMs: Each subject and average across subjects
% Save the RDM matrix plot into the RDMs folder instead of the 1st order
% one

% Plot Neural RDMs--------------------------------
% Avg subjects-------------------
% all subjects
figureRDMs(concatenateRDMs(neural_RDMs_avg_all), user_options_visualize_RDM,...
    struct('fileName', 'neural_RDMs_avg_all'));
% 2nd grade
figureRDMs(concatenateRDMs(neural_RDMs_avg_2_grade), user_options_visualize_RDM,...
    struct('fileName', 'neural_RDMs_avg_2_grade'));
% 5th grade
figureRDMs(concatenateRDMs(neural_RDMs_avg_5_grade), user_options_visualize_RDM,...
    struct('fileName', 'neural_RDMs_avg_5_grade' ));
% For each VOI-------------------------------------
% all subjects----------------
for v=1:size(list_VOI_name,1)
    figureRDMs(concatenateRDMs(neural_RDMs(v,:)), user_options_visualize_RDM,...
        struct('fileName', ['neural_RDMs_',char(list_VOI_name{v}),'_all']));
end
% 2nd grade----------------
for v=1:size(list_VOI_name,1)
    figureRDMs(concatenateRDMs(neural_RDMs(v,strcmp(table_demographic.grade,'2'))), user_options_visualize_RDM,...
        struct('fileName', ['neural_RDMs_',char(list_VOI_name{v}),'_2_grade']));
end
% 5th grade----------------
for v=1:size(list_VOI_name,1)
    figureRDMs(concatenateRDMs(neural_RDMs(v,strcmp(table_demographic.grade,'5'))), user_options_visualize_RDM,...
        struct('fileName', ['neural_RDMs_',char(list_VOI_name{v}),'_5_grade']));
end

% Plot Conceptual RDMs--------------------------------
figureRDMs(concatenateRDMs(conceptual_RDMs), user_options_visualize_RDM,...
    struct('fileName', 'conceptual_RDMs', 'figureNumber', 1));

%% (Step 2) Determine dendrograms for the clustering of the conditions
%%%%%%%%%%%%%WIERD! The condition color and condition numbers do not match!
user_options_visualize_RDM.rootPath=path_1st_result; %Save the results starting from here to the 1st order part
% (all/2nd grade/5th grade) x  VOIs----------------
for group_index=1:numel(list_group_name)
    dendrogramConditions(list_neural_RDMs{group_index}, user_options_visualize_RDM,...
        struct(...
        'useAlternativeConditionLabels', true,...
        'alternativeConditionLabels', {condition_labels}));
    % The pdf file name is et as default, so I need to modify them into
    % desired names.
    for v=1:numel(list_VOI_name)
        old_name = [user_options_visualize_RDM.analysisName,'_Dendrogram_',...
            regexprep(char(list_VOI_name{v}),'_','-'),'.pdf'];
        path_old_name = fullfile(user_options_visualize_RDM.rootPath,old_name);
        new_name= [user_options_visualize_RDM.analysisName,'_Dendrogram_',char(list_VOI_name{v}),...
            '_',list_group_name{group_index},'.pdf'];
        path_new_name = fullfile(user_options_visualize_RDM.rootPath,new_name);
        movefile(path_old_name,path_new_name);
    end
end

%% (Step 3) MDS of conditions
% (all/2nd grade/5th grade/child/adult)  x  VOIs----------------
for group_index=1:numel(list_group_name)
    MDSConditions(list_neural_RDMs{group_index}, user_options_visualize_RDM,...
        struct('alternativeConditionLabels', {condition_labels}));
    % The pdf file name is et as default, so I need to modify them into
    % desired names.
    for v=1:numel(list_VOI_name)
        old_name = [user_options_visualize_RDM.analysisName,'_FirstOrderMDS_',...
            regexprep(char(list_VOI_name{v}),'_','-'),'.pdf'];
        path_old_name = fullfile(user_options_visualize_RDM.rootPath,old_name);
        new_name= [user_options_visualize_RDM.analysisName,'_FirstOrderMDS_',char(list_VOI_name{v}),...
            '_',list_group_name{group_index},'.pdf'];
        path_new_name = fullfile(user_options_visualize_RDM.rootPath,new_name);
        movefile(path_old_name,path_new_name);
    end
end

file_save=...
    fullfile(user_options_visualize_RDM.rootPath,['user_options_visualize_RDM_',version_num,'_',date,'.mat']);

save(file_save,'user_options_visualize_RDM')
