%% Part 3 - iv-c: RSA : 2nd order analysis - MDS
% trial 5 in the README.txt.
% Part3: 2nd order analysis - MDS
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

% c: Include adult neural data (both collapse split by age group)
% Split into 2nd grader, 5th grader, and the adult

% [Run info has already been averaged out. Therefore, the number of runs no longer matters.]
PATH_ROOT='D:\Yun-Shiuan_LAMBDA';
PATH_TOOL_CODE='D:\GoogleDrive\Lambda_code\m_file\LAMBDA\tool_code';
PATH_HELPER_FUNCTION='D:\GoogleDrive\Lambda_code\m_file\LAMBDA\helper_function';
addpath(genpath(fullfile(PATH_ROOT,'rsatoolbox')));
addpath(PATH_TOOL_CODE);% enable read_mixed_csv()
addpath(genpath(PATH_HELPER_FUNCTION));% enable my rsa helper functions


%% Set up
% Constants----------------------------------------------
% path-------------------------
PATH_1ST_RESULT=fullfile(PATH_ROOT,'RSA','trial_13','Part2_1st_oder_analysis');
PATH_2ND_RESULT=fullfile(PATH_ROOT,'RSA','trial_13','Part3_2nd_order_analysis','base_models');
% file--------------------------
FILE_RDM_RESULT=fullfile(PATH_ROOT,'RSA','trial_13','Part1_neural_and_conceptual_RDMs',...
                        'RDMs','RDMs_and_options_neural_and_conceptual_trial_13_18-Apr-2018.mat');
FILE_1ST_USER_OPTION=fullfile(PATH_1ST_RESULT,'user_options_visualize_RDM_13_18-Apr-2018.mat');
FILE_VALID_RUN_CHILD=fullfile(PATH_ROOT,'Run_inclusion_info','inclusive_runs_indexes_new_April_5.csv');
FILE_VALID_RUN_ADULT=fullfile(PATH_ROOT,'Adult','Run_inclusion_info','inclusive_runs_indexes.csv');
FILE_DEMOGRAPHIC_CHILD=fullfile(PATH_ROOT,'demographic','tidy_demographic.csv');
FILE_DEMOGRAPHIC_ADULT=fullfile(PATH_ROOT,'Adult','demographic','tidy_demographic.csv');
% FILE_SETUP='D:\GoogleDrive\Lambda_code\m_file\LAMBDA\RSA_2nd_order_analysis_t_value_discrete_18_include_adult_set_up.m';
% Parameters----------------------
LIST_VOI_NAME={"R_V1","L_V1",...
               "R_IPS","L_IPS","R_DLPFC"}';
VERSION_NUM='13';
VERSION_DATE= ['RSA_trial_',VERSION_NUM,'_',date];
LIST_AGE_GROUP_NAME={'2_grade','5_grade','adult','all'}; %No need to collapse children or all subjects
% For accuracy models (note that 'all' leads to nothing, it is included
% only for making the dimension consistent
LIST_AGE_ABRV={'S','F','A','all'}; 
LIST_MDS_CRITERION={'stress','metricstress'}; % The options for MDS
BOOLEAN_ACC_INCLUDE=false;%Whether accuracy models should be included (could be adjusted)
PATTERN_BASE_CONCEPTUAL_RDMS=['^g$|^n$|',...%Format
    '^a$|^s$|',...%Distance
    '^agv1$|^agv2$|^agv3$|^sgv1$|^sgv2$|',...%log distance
    '^null$'];%null
if BOOLEAN_ACC_INCLUDE
    PATH_2ND_RESULT=fullfile(PATH_2ND_RESULT,'with_ACC');
else
    PATH_2ND_RESULT=fullfile(PATH_2ND_RESULT,'without_ACC');
end
PATH_2ND_RESULT=fullfile(PATH_2ND_RESULT,'MDS');
mkdir(PATH_2ND_RESULT);
% run(FILE_SETUP); %Remember to adjust the parameters in the setup script


parameters_for_set_up=struct('file_user_options',FILE_1ST_USER_OPTION,...
    'file_rdm_result',FILE_RDM_RESULT,...
    'file_valid_run_child',FILE_VALID_RUN_CHILD,...
    'file_valid_run_adult',FILE_VALID_RUN_ADULT,...
    'file_demographic_child',FILE_DEMOGRAPHIC_CHILD,...
    'file_demographic_adult',FILE_DEMOGRAPHIC_ADULT,...
    'path_2nd_result',PATH_2ND_RESULT,...
    'version_date',VERSION_DATE,...
    'pattern_reduced_conceptual_RDMs',PATTERN_BASE_CONCEPTUAL_RDMS,...
    'list_age_abrv',{LIST_AGE_ABRV},...
    'list_age_group_name',{LIST_AGE_GROUP_NAME},...
    'acc_include',BOOLEAN_ACC_INCLUDE);

[preprocessed_parameters]=RSA_setUpSecondOrderEnvironment(parameters_for_set_up);

list_reduced_conceptual_RDMs=preprocessed_parameters.list_reduced_conceptual_RDMs;
user_options_2nd_order=preprocessed_parameters.user_options_2nd_order;
list_neural_avg_RDMs=preprocessed_parameters.list_neural_avg_RDMs;
list_group_index=preprocessed_parameters.list_group_index;
neural_RDMs=preprocessed_parameters.neural_RDMs;


user_options_2nd_order.RDMcorrelationType='Kendall_taua';
user_options_2nd_order.distanceMeasure='Kendall_taua';


%% 2nd order RDMs and Visualization: 2nd order simmilarity matrix and 2nd order MDS
% Display a 2nd order simmilarity matrix between neural RDMs and conceptual
% RDMs             
for group_index=1:numel(LIST_AGE_GROUP_NAME)
    pairwiseCorrelateRDMs({list_neural_avg_RDMs{group_index},...
                           list_reduced_conceptual_RDMs{group_index}},...
                           user_options_2nd_order);
    % The pdf file name is et as default, so I need to modify them into
    % desired names.
    old_name = [user_options_2nd_order.analysisName,'_secondOrderSM.pdf'];
    path_old_name = fullfile(user_options_2nd_order.rootPath,old_name);
    new_name= [user_options_2nd_order.analysisName,'_secondOrderSM_',...
        LIST_AGE_GROUP_NAME{group_index},'.pdf'];
    path_new_name = fullfile(user_options_2nd_order.rootPath,new_name);
    movefile(path_old_name,path_new_name);
end

% MDS: Plot all models as an MDS plot to visualise pairwise
% distances.-----------------------------------------------------------
% Reduce the conceptual RDMs to be plotted.
% Without reducing: Points in the configuration have co-located.  Try a different
% starting point, or use a different criterion.
user_options_2nd_order.dotSize=8;
user_options_2nd_order.saveFiguresPDF=1;
user_options_2nd_order.saveFigurePDF=1;

for MDS_criterion=1:numel(LIST_MDS_CRITERION) %Test different MDS criteria
    name_MDS_type=LIST_MDS_CRITERION{MDS_criterion};
    user_options_2nd_order.criterion=name_MDS_type;
    path_MDS_type=fullfile(user_options_2nd_order.rootPath,name_MDS_type);
    mkdir(path_MDS_type)
    
    for group_index=1:numel(LIST_AGE_GROUP_NAME)
        MDS_results=MDSRDMsReturnResults({list_neural_avg_RDMs{group_index},...
                    list_reduced_conceptual_RDMs{group_index}}, user_options_2nd_order,...
                    struct('titleString', 'MDS of neural and conceptual'));
        name_group=LIST_AGE_GROUP_NAME{group_index};
        
        %Save the MDS results in the structure
        valid_name_group=matlab.lang.makeValidName(name_group);
        valid_name_MDS_type=matlab.lang.makeValidName(name_MDS_type);
        
        list_MDS_results.(valid_name_MDS_type).(valid_name_group)=MDS_results;
        
        % The pdf file name is set as default, so I need to modify them into
        % desired names.
        old_name = [user_options_2nd_order.analysisName,'_SecondOrderMDS_.pdf'];
        path_old_name = fullfile(user_options_2nd_order.rootPath,old_name);
        new_name= [user_options_2nd_order.analysisName,'_secondOrderMDS_',...
            name_group,'.pdf'];
        path_new_name = fullfile(user_options_2nd_order.rootPath,name_MDS_type,new_name);
        movefile(path_old_name,path_new_name);
    end
end
% If one want to save manually instead of using the automatic saving function of the
% MDSRDMs()
% saveas(gcf,fullfile(path_2nd_result,[user_options_2nd_order.analysisName,'_SecondOrderMDS_.pdf']))

%Save the MDS results for visualization in R
file_MDS_result_save=...
        fullfile(user_options_2nd_order.rootPath,['MDS_Results_',VERSION_NUM,'_',date,'.mat']);
save(file_MDS_result_save,'list_MDS_results');

file_MDS_user_option_save=...
        fullfile(user_options_2nd_order.rootPath,['MDS_User_Option_',VERSION_NUM,'_',date,'.mat']);
save(file_MDS_user_option_save,'user_options_2nd_order');






