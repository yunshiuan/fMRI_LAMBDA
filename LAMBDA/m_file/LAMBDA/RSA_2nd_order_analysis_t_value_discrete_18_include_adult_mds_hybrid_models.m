%% Part 5 - iv-c: RSA : 2nd order analysis - Hybrid best base models - MDS 
% trial 5 in the README.txt.
% Part5: 2nd order analysis - Hibrid best format model with best magnitude
% model for each age group ROI - MDS 
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Update the subject list for children: use "inclusive_runs_indexes_new_Jun_10.csv" (2018/6/10)
%Update the subject list for children: use "inclusive_runs_indexes_new_Jun_11.csv" (2018/6/10)

PATH_ROOT='D:\Yun-Shiuan_LAMBDA';
PATH_TOOL_CODE='D:\GoogleDrive\Lambda_code\m_file\LAMBDA\tool_code';
PATH_HELPER_FUNCTION='D:\GoogleDrive\Lambda_code\m_file\LAMBDA\helper_function';
addpath(genpath(fullfile(PATH_ROOT,'rsatoolbox')));
addpath(PATH_TOOL_CODE);% enable read_mixed_csv()
addpath(genpath(PATH_HELPER_FUNCTION));% enable my rsa helper functions


%% Set up
% Constants----------------------------------------------
% Parameters----------------------

VERSION_NUM='22';
VERSION_DATE= ['RSA_trial_',VERSION_NUM,'_',date];
% LIST_VOI_NAME={'R_V1','L_V1',...
%                'R_IPS','L_IPS','R_DLPFC'}';
% LIST_BEST_BASE_MAG={'sgv2','sgv2',....
%                     'agv2','agv2','agv2'}';
% LIST_BEST_BASE_FORMAT={'g','g',....
%                        'g','g','g'}';
LIST_VOI_NAME={'R-V1','L-V1',...
               'R-IPS','L-IPS'}';
LIST_BEST_BASE_MAG={'sgv2','sgv2',....
                    'agv2','agv2'}';
LIST_BEST_BASE_FORMAT={'g','g',....
                       'g','g'}';
LIST_BASE_MODEL_TYPE={'mag_base','format_base'}';
LIST_AGE_GROUP_NAME={'2_grade','5_grade','adult','all'}'; 
LIST_AGE_GROUP_NAME_VALID=matlab.lang.makeValidName(LIST_AGE_GROUP_NAME);
LIST_AGE_ABRV={'S','F','A'}; % For accuracy models
LIST_MDS_CRITERION={'stress','metricstress'}; % The options for MDS
LIST_ACC_VERSION={'without_ACC'}'; %Ignore the Accuracy since it does not count as format or magnitude
BOOLEAN_ACC_INCLUDE=false;%Whether accuracy models should be included (could be adjusted)


VAR_NAMES={'VOI_name','age_group','base_type','best_model'}';
NUM_ROWS=numel(LIST_VOI_NAME)*numel(LIST_AGE_GROUP_NAME)*numel(LIST_BASE_MODEL_TYPE);
TABLE_BEST_BASE_MODELS=table(...
    rep_to_length(LIST_VOI_NAME,NUM_ROWS),...
    rep_to_length(rep_each_times(LIST_AGE_GROUP_NAME,numel(LIST_VOI_NAME),1),NUM_ROWS),...
    rep_to_length(rep_each_times(...
        LIST_BASE_MODEL_TYPE,...
        numel(LIST_VOI_NAME)*numel(LIST_AGE_GROUP_NAME),1),...
        NUM_ROWS),...
    [rep_each_times(...
        LIST_BEST_BASE_MAG,...
        1,numel(LIST_AGE_GROUP_NAME));...
    rep_each_times(...
        LIST_BEST_BASE_FORMAT,...
        1,numel(LIST_AGE_GROUP_NAME))],...
     'VariableNames',VAR_NAMES);
PATTERN_REDUCED_HYBRID_CONCEPTUAL_RDMS=['^gagv2(1|4|7)0|','^gsgv2(1|4|7)0|',...
    '^g$|','^agv2$|','^sgv2$'];

% path-------------------------
PATH_2ND_RESULT_HYBRID=fullfile(PATH_ROOT,'RSA','trial_22','Part3_2nd_order_analysis','hybrid_models');
PATH_2ND_RESULT_BASE=fullfile(PATH_ROOT,'RSA','trial_13','Part3_2nd_order_analysis','base_models');

if BOOLEAN_ACC_INCLUDE
    PATH_2ND_RESULT_HYBRID=fullfile(PATH_2ND_RESULT_HYBRID,'with_ACC');
else
    PATH_2ND_RESULT_HYBRID=fullfile(PATH_2ND_RESULT_HYBRID,'without_ACC');
end
PATH_2ND_RESULT_HYBRID=fullfile(PATH_2ND_RESULT_HYBRID,'MDS');
mkdir(PATH_2ND_RESULT_HYBRID);

% file--------------------------
FILE_RDM_RESULT=fullfile(PATH_ROOT,'RSA','trial_22','Part1_neural_and_conceptual_RDMs',...
    'RDMs','RDMs_and_options_neural_and_conceptual_trial_22_10-Jun-2018.mat');
FILE_2ND_MDS_USER_OPTION=fullfile(PATH_2ND_RESULT_BASE,LIST_ACC_VERSION,'MDS',...
    'MDS_User_Option_13_27-Apr-2018.mat');
% FILE_2ND_RELATEDNESS_TEST_USER_OPTION=fullfile(PATH_2ND_RESULT_BASE,LIST_ACC_VERSION,'relatedness_test',...
%     'Relatedness_Test_User_Option_13_18-Apr-2018.mat');
FILE_VALID_RUN_CHILD=fullfile(PATH_ROOT,'Run_inclusion_info','inclusive_runs_indexes_new_June_11.csv');
FILE_VALID_RUN_ADULT=fullfile(PATH_ROOT,'Adult','Run_inclusion_info','inclusive_runs_indexes.csv');
FILE_DEMOGRAPHIC_CHILD=fullfile(PATH_ROOT,'demographic','tidy_demographic.csv');
FILE_DEMOGRAPHIC_ADULT=fullfile(PATH_ROOT,'Adult','demographic','tidy_demographic.csv');


%Set up necessary parameters
parameters_for_set_up=struct('file_user_options',FILE_2ND_MDS_USER_OPTION,...
    'file_rdm_result',FILE_RDM_RESULT,...
    'file_valid_run_child',FILE_VALID_RUN_CHILD,...
    'file_valid_run_adult',FILE_VALID_RUN_ADULT,...
    'file_demographic_child',FILE_DEMOGRAPHIC_CHILD,...
    'file_demographic_adult',FILE_DEMOGRAPHIC_ADULT,...
    'path_2nd_result',PATH_2ND_RESULT_HYBRID,...
    'version_date',VERSION_DATE,...
    'pattern_reduced_conceptual_RDMs',PATTERN_REDUCED_HYBRID_CONCEPTUAL_RDMS,...
    'list_age_abrv',{LIST_AGE_ABRV},...
    'list_age_group_name',{LIST_AGE_GROUP_NAME},...
    'acc_include',BOOLEAN_ACC_INCLUDE,...
    'list_interested_brain_region',{LIST_VOI_NAME});

[preprocessed_parameters]=RSA_setUpSecondOrderEnvironment(parameters_for_set_up);

list_reduced_conceptual_RDMs=preprocessed_parameters.list_reduced_conceptual_RDMs;
user_options_2nd_order=preprocessed_parameters.user_options_2nd_order;
list_neural_avg_RDMs=preprocessed_parameters.list_neural_avg_RDMs;
list_group_index=preprocessed_parameters.list_group_index;
neural_RDMs=preprocessed_parameters.neural_RDMs;


user_options_2nd_order.RDMcorrelationType='Kendall_taua';
user_options_2nd_order.distanceMeasure='Kendall_taua';


%% 2nd order MDS for each VOI x age group
% Display a 2nd order simmilarity matrix between neural RDMs and conceptual
% RDMs             
% for group_index=1:numel(LIST_AGE_GROUP_NAME)
%     
%     pairwiseCorrelateRDMs({list_neural_avg_RDMs{group_index},...
%         list_reduced_conceptual_RDMs{group_index}},...
%         user_options_2nd_order);
%     % The pdf file name is et as default, so I need to modify them into
%     % desired names.
%     old_name = [user_options_2nd_order.analysisName,'_secondOrderSM.pdf'];
%     path_old_name = fullfile(user_options_2nd_order.rootPath,old_name);
%     new_name= [user_options_2nd_order.analysisName,'_secondOrderSM_',...
%         LIST_AGE_GROUP_NAME{group_index},'.pdf'];
%     path_new_name = fullfile(user_options_2nd_order.rootPath,new_name);
%     movefile(path_old_name,path_new_name);
% end

% MDS: Plot all models as an MDS plot to visualise pairwise
% distances.-----------------------------------------------------------
% Reduce the conceptual RDMs to be plotted.
% Without reducing: Points in the configuration have co-located.  Try a different
% starting point, or use a different criterion.
user_options_2nd_order.dotSize=8;
user_options_2nd_order.saveFiguresPDF=0;
user_options_2nd_order.saveFigurePDF=0;

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
        if logical(user_options_2nd_order.saveFigurePDF)
            old_name = [user_options_2nd_order.analysisName,'_SecondOrderMDS_.pdf'];
            path_old_name = fullfile(user_options_2nd_order.rootPath,old_name);
            new_name= [user_options_2nd_order.analysisName,'_secondOrderMDS_',...
                name_group,'.pdf'];
            path_new_name = fullfile(user_options_2nd_order.rootPath,name_MDS_type,new_name);
            movefile(path_old_name,path_new_name);
        end
    end
end

%Save the MDS results for visualization in R
file_MDS_result_save=...
        fullfile(user_options_2nd_order.rootPath,['MDS_Results_',VERSION_NUM,'_',date,'.mat']);
save(file_MDS_result_save,'list_MDS_results');

file_MDS_user_option_save=...
        fullfile(user_options_2nd_order.rootPath,['MDS_User_Option_',VERSION_NUM,'_',date,'.mat']);
save(file_MDS_user_option_save,'user_options_2nd_order');



