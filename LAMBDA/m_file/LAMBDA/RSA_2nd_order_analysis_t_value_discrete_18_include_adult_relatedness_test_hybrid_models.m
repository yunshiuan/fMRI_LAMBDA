%% Part 6 - iv-c: RSA : 2nd order analysis - Hybrid best base models - Relatedness Test
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Redefine the scale and decrease the increment to 0.01.
% Yet for the relatedness plot, still only include 9 mixed models per combination
% (i.e., w=0.10, 0.20,...). This is because I do not want to include too 
% many candidate RDMs in the plot. In contrast, for tau-profile analysis,
% all 99 mixed models per combination are included. (2018/6/18)

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
LIST_VOI_NAME={'R-V1','L-V1',...
               'R-IPS','L-IPS'}';
% LIST_BEST_BASE_MAG={'sgv2','sgv2',....
%                     'agv2','agv2'}';
% LIST_BEST_BASE_FORMAT={'g','g',....
%                        'g','g'}';
% % LIST_VOI_NAME={'R_V1','L_V1',...
%                'R_IPS','L_IPS','R_DLPFC'}';
% LIST_BEST_BASE_MAG={'sgv2','sgv2',....
%                     'agv2','agv2','agv2'}';
% LIST_BEST_BASE_FORMAT={'g','g',....
%                        'g','g','g'}';
LIST_BASE_MODEL_TYPE={'mag_base','format_base'}';
LIST_AGE_GROUP_NAME={'2_grade','5_grade','adult','all'}'; 
LIST_AGE_GROUP_NAME_VALID=matlab.lang.makeValidName(LIST_AGE_GROUP_NAME);
LIST_AGE_ABRV={'S','F','A'}; % For accuracy models
LIST_MDS_CRITERION={'stress','metricstress'}; % The options for MDS
LIST_ACC_VERSION={'without_ACC'}'; %Ignore the Accuracy since it does not count as format or magnitude
BOOLEAN_ACC_INCLUDE=false;%Whether accuracy models should be included (could be adjusted)

VAR_NAMES={'VOI_name','age_group','base_type','best_model'}';
NUM_ROWS=numel(LIST_VOI_NAME)*numel(LIST_AGE_GROUP_NAME)*numel(LIST_BASE_MODEL_TYPE);
% TABLE_BEST_BASE_MODELS=table(...
%     rep_to_length(LIST_VOI_NAME,NUM_ROWS),...
%     rep_to_length(rep_each_times(LIST_AGE_GROUP_NAME,numel(LIST_VOI_NAME),1),NUM_ROWS),...
%     rep_to_length(rep_each_times(...
%         LIST_BASE_MODEL_TYPE,...
%         numel(LIST_VOI_NAME)*numel(LIST_AGE_GROUP_NAME),1),...
%         NUM_ROWS),...
%     [rep_each_times(...
%         LIST_BEST_BASE_MAG,...
%         1,numel(LIST_AGE_GROUP_NAME));...
%     rep_each_times(...
%         LIST_BEST_BASE_FORMAT,...
%         1,numel(LIST_AGE_GROUP_NAME))],...
%      'VariableNames',VAR_NAMES);
PATTERN_REDUCED_HYBRID_CONCEPTUAL_RDMS=['^gagv2\d0|','^gsgv2\d0|','^g$|',...
    '^agv2$|','^sgv2$'];

% path-------------------------
%output
PATH_2ND_RESULT_HYBRID=fullfile(PATH_ROOT,'RSA','trial_22','Part3_2nd_order_analysis','hybrid_models');
PATH_2ND_RESULT_BASE=fullfile(PATH_ROOT,'RSA','trial_22','Part3_2nd_order_analysis','base_models');

if BOOLEAN_ACC_INCLUDE
    PATH_2ND_RESULT_HYBRID=fullfile(PATH_2ND_RESULT_HYBRID,'with_ACC');
else
    PATH_2ND_RESULT_HYBRID=fullfile(PATH_2ND_RESULT_HYBRID,'without_ACC');
end
PATH_2ND_RESULT_HYBRID=fullfile(PATH_2ND_RESULT_HYBRID,'relatedness_test');
mkdir(PATH_2ND_RESULT_HYBRID);

% file--------------------------
%input
FILE_RDM_RESULT=fullfile(PATH_ROOT,'RSA','trial_22','Part1_neural_and_conceptual_RDMs',...
    'RDMs','RDMs_and_options_neural_and_conceptual_trial_22_10-Jun-2018.mat');
% FILE_2ND_MDS_USER_OPTION=fullfile(PATH_2ND_RESULT_BASE,LIST_ACC_VERSION,'MDS',...
%     'MDS_User_Option_13_27-Apr-2018.mat');
FILE_2ND_RELATEDNESS_TEST_USER_OPTION=fullfile(PATH_2ND_RESULT_BASE,LIST_ACC_VERSION,'relatedness_test',...
    'Relatedness_Test_User_Option_22_10-Jun-2018.mat');
FILE_VALID_RUN_CHILD=fullfile(PATH_ROOT,'Run_inclusion_info','inclusive_runs_indexes_new_June_11.csv');
FILE_VALID_RUN_ADULT=fullfile(PATH_ROOT,'Adult','Run_inclusion_info','inclusive_runs_indexes.csv');
FILE_DEMOGRAPHIC_CHILD=fullfile(PATH_ROOT,'demographic','tidy_demographic.csv');
FILE_DEMOGRAPHIC_ADULT=fullfile(PATH_ROOT,'Adult','demographic','tidy_demographic.csv');


%Set up necessary parameters
parameters_for_set_up=struct('file_user_options',FILE_2ND_RELATEDNESS_TEST_USER_OPTION,...
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

%Set up options for 2nd test
user_options_2nd_order.RDMcorrelationType='Kendall_taua';
user_options_2nd_order.distanceMeasure='Kendall_taua';
user_options_2nd_order.RDMrelatednessTest = 'subjectRFXsignedRank';
user_options_2nd_order.RDMrelatednessThreshold = 0.05;
user_options_2nd_order.figureIndex = [1:2];
user_options_2nd_order.RDMrelatednessMultipleTesting = 'FDR';
user_options_2nd_order.candRDMdifferencesTest = 'subjectRFXsignedRank'; %none
user_options_2nd_order.candRDMdifferencesThreshold = 0.05;
user_options_2nd_order.candRDMdifferencesMultipleTesting = 'none';

user_options_2nd_order.RDMcorrelationType='Kendall_taua';
user_options_2nd_order.distanceMeasure='Kendall_taua';

%% Testing hybrid models: 
% stats_p_r = compareRefRDM2candRDMs(refRDM, candRDMs[, userOptions]) - statistical inference:
 % The scope of the base models search scope (e.g., magnitude base models)
 
% LIST_TEST_SCOPE_BASE_MODELS={'all_base','mag_base','format_base'};
    
% test the relatedness and compare the candidate RDMs
name_figure1_filename_base='BarGraph';
name_figure2_filename_base='CompareCandidateRDMs';

%Number of iteration
num_cand_conceptual=numel(list_reduced_conceptual_RDMs);
num_cand_neural=numel(LIST_VOI_NAME);

%To collect the candidates RDMs for each age gruop
candidate_RDMs=cell(numel(num_cand_conceptual+num_cand_neural),1);

%To collect the relatedness test results
% collect_stats_p_r: a structure storing the results

user_options_2nd_order.saveFiguresPDF=0;
user_options_2nd_order.saveFigurePDF=0;

% for test_scope_index=1:numel(LIST_TEST_SCOPE_BASE_MODELS)
%
%     name_scope=LIST_TEST_SCOPE_BASE_MODELS{test_scope_index};

%The path where result pdf wii print out to
user_options_2nd_order.resultsPath=user_options_2nd_order.rootPath;
mkdir(user_options_2nd_order.resultsPath);

for group_index=1:numel(LIST_AGE_GROUP_NAME)
    
    reduced_conceptual_RDMs_this_group=list_reduced_conceptual_RDMs{group_index};
    
    base_all_indexes=...
        true(size(reduced_conceptual_RDMs_this_group));
    % No need to consider scope for the hybrid models, since not
    % distinguising the format and mag models.
%     base_mag_indexes=...
%         ismember({reduced_conceptual_RDMs_this_group.name},NAMES_MAG_BASE);
%     base_format_indexes=...
%         ismember({reduced_conceptual_RDMs_this_group.name},PATTERN_FORMAT_BASE);
%     list_test_scope_base_models_indexes=...
%         {base_all_indexes,base_mag_indexes,base_format_indexes};
    
%     scope_indexes=list_test_scope_base_models_indexes{test_scope_index};
    
    
    name_group=LIST_AGE_GROUP_NAME{group_index};
    
    % Retrieve conceptual models
    % when include_ACC=true, include difficulty models (age group-specific)
    reduced_conceptual_RDMs=list_reduced_conceptual_RDMs{group_index};
%     
%     reduced_conceptual_RDMs=reduced_conceptual_RDMs(scope_indexes);
    
    clear candidate_RDMs
    
    % Candidate Models: Conceptual RDMS (averaged across subjects)
    for conceptual_index=1:numel(reduced_conceptual_RDMs)
        candidate_RDMs{conceptual_index}=reduced_conceptual_RDMs(conceptual_index);
    end
    
    % Candidate Models: Neural RDMS (averaged across subjects)
    for neural_index=1:numel(LIST_VOI_NAME)
        candidate_RDMs{conceptual_index+neural_index}=list_neural_avg_RDMs{group_index}(neural_index);
    end
    
    for v=1:numel(LIST_VOI_NAME)
        % The reference model (the neural RDMs from every subjects)
        reference_neural_RDMs=neural_RDMs(v,list_group_index{group_index});
        
        %Figures name
        name_VOI=LIST_VOI_NAME{v};
        
        name_figure1_filename=[name_figure1_filename_base,'_',...
            char(name_group),'_',char(name_VOI)];
        name_figure2_filename=[name_figure2_filename_base,'_',...
            char(name_group),'_',char(name_VOI)];
        
        user_options_2nd_order.figure1filename=name_figure1_filename;
        user_options_2nd_order.figure2filename=name_figure2_filename;
        
        % Test!
        stats_p_r=compareRefRDM2candRDMs(reference_neural_RDMs, ...
            candidate_RDMs, user_options_2nd_order);
        % Save the results
%         valid_name_scope=matlab.lang.makeValidName(char(name_scope));
        valid_name_group=matlab.lang.makeValidName(char(name_group));
        valid_name_VOI=matlab.lang.makeValidName(char(name_VOI));
        
        collect_stats_p_r.(valid_name_group).(valid_name_VOI)=...
            stats_p_r;
        
    end
end
% end

%Save the relatedness test result
file_test_result_save=...
    fullfile(user_options_2nd_order.rootPath,['Relatedness_Test_Results_',VERSION_NUM,'_',date,'.mat']);
save(file_test_result_save,'collect_stats_p_r');


try
    user_options_2nd_order=...
        rmfield(user_options_2nd_order,{'figure1filename'});
    user_options_2nd_order=...
        rmfield(user_options_2nd_order,{'figure2filename'});
catch
    try
        user_options_2nd_order=...
            rmfield(user_options_2nd_order,{'figure2filename'});
    catch
    end
end

file_option_2nd_order=...
    fullfile(user_options_2nd_order.rootPath,['Relatedness_Test_User_Option_',VERSION_NUM,'_',date,'.mat']);
save(file_option_2nd_order,'user_options_2nd_order')
