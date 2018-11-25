%% RSA : 2nd order analysis - Hybrid models relatedness test for the tau-profile -  Bootstrapping

%Bootstrapping the tau profile
% (1)2nd order tau profile (tau ~ weigh):  
%     (1)[MATLAB] Bootstrap to construct a set of tau profile for each age group x brain region.
%             -full sample: 18 x 2nd/ 20 x 5th /24 x Adult
%             -bootstrap procedure: 
%                     -sampling with replacement (resample 1000 times)
%                             2nd: sample(size=18,replacement=T) , repeate 1000 times
%                             5th: sample(size=20,replacement=T) , repeate 1000 times
%                             Adult: sample(size=24,replacement=T) , repeate 1000 times
%                     -bootstrapped sample: 1000 taus for each age group
%     (2)[R]Compute the delta_area (i.e.,right - left) for each profile per age group x brain region.
%     (2)[R]Derive the distribution of the delta_area per age group x brain region.

%Note: Only focus on GxA (but not GxS).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
NUM_BOOT_STRAP=1000;
LIST_VOI_NAME={'R-V1','L-V1',...
               'R-IPS','L-IPS'}';
LIST_BASE_MODEL_TYPE={'mag_base','format_base'}';

LIST_AGE_GROUP_NAME={'2_grade','5_grade','adult'}'; 
% LIST_AGE_GROUP_NAME={'2_grade','5_grade','adult','all'}'; 

LIST_AGE_GROUP_NAME_VALID=matlab.lang.makeValidName(LIST_AGE_GROUP_NAME);
LIST_AGE_ABRV={'S','F','A'}; % For accuracy models

LIST_MDS_CRITERION={'stress'}; % The options for MDS
% LIST_MDS_CRITERION={'stress','metricstress'}; % The options for MDS

LIST_ACC_VERSION={'without_ACC'}'; %Ignore the Accuracy since it does not count as format or magnitude
BOOLEAN_ACC_INCLUDE=false;%Whether accuracy models should be included (could be adjusted)

VAR_NAMES={'VOI_name','age_group','base_type','best_model'}';
NUM_ROWS=numel(LIST_VOI_NAME)*numel(LIST_AGE_GROUP_NAME)*numel(LIST_BASE_MODEL_TYPE);

PATTERN_REDUCED_HYBRID_CONCEPTUAL_RDMS=['^gagv2\d0|','^g$|',...
    '^agv2$'];
% PATTERN_REDUCED_HYBRID_CONCEPTUAL_RDMS=['^gagv2\d0|','^gsgv2\d0|','^g$|',...
%     '^agv2$|','^sgv2$'];

% path-------------------------
%output
PATH_2ND_RESULT_HYBRID_BOOTSTRAP=fullfile(PATH_ROOT,...
    'RSA','trial_22','Part3_2nd_order_analysis','hybrid_models','bootstrap');
PATH_2ND_RESULT_BASE=fullfile(PATH_ROOT,...
    'RSA','trial_13','Part3_2nd_order_analysis','base_models');

if BOOLEAN_ACC_INCLUDE
    PATH_2ND_RESULT_HYBRID_BOOTSTRAP=fullfile(PATH_2ND_RESULT_HYBRID_BOOTSTRAP,'with_ACC');
else
    PATH_2ND_RESULT_HYBRID_BOOTSTRAP=fullfile(PATH_2ND_RESULT_HYBRID_BOOTSTRAP,'without_ACC');
end
PATH_2ND_RESULT_HYBRID_BOOTSTRAP=fullfile(PATH_2ND_RESULT_HYBRID_BOOTSTRAP,'relatedness_test');
mkdir(PATH_2ND_RESULT_HYBRID_BOOTSTRAP);

% file--------------------------
%input
FILE_RDM_RESULT=fullfile(PATH_ROOT,'RSA','trial_22','Part1_neural_and_conceptual_RDMs',...
    'RDMs','RDMs_and_options_neural_and_conceptual_trial_22_10-Jun-2018.mat');
% FILE_2ND_MDS_USER_OPTION=fullfile(PATH_2ND_RESULT_BASE,LIST_ACC_VERSION,'MDS',...
%     'MDS_User_Option_13_27-Apr-2018.mat');
FILE_2ND_RELATEDNESS_TEST_USER_OPTION=fullfile(PATH_2ND_RESULT_BASE,LIST_ACC_VERSION,'relatedness_test',...
    'Relatedness_Test_User_Option_13_18-Apr-2018.mat');
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
    'path_2nd_result',PATH_2ND_RESULT_HYBRID_BOOTSTRAP,...
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
% user_options_2nd_order.candRDMdifferencesTest='none';
% for test_scope_index=1:numel(LIST_TEST_SCOPE_BASE_MODELS)
%
%     name_scope=LIST_TEST_SCOPE_BASE_MODELS{test_scope_index};

%The path where result pdf wii print out to
user_options_2nd_order.resultsPath=user_options_2nd_order.rootPath;
mkdir(user_options_2nd_order.resultsPath);

for group_index=1:numel(LIST_AGE_GROUP_NAME)
    name_group=LIST_AGE_GROUP_NAME{group_index};
       
        
    for boot_strap_index=1:5%NUM_BOOT_STRAP        
        
        %Bootstrap the subject index for each iteration
        subject_index_find=find(list_group_index{group_index});
        n_subject_this_age_group=numel(subject_index_find);
        n_subject_all=numel(list_group_index{1});
        
        %Seed for bootstrapping (to ensure reproducibility)
        seed = RandStream.create('mrg32k3a','NumStreams',NUM_BOOT_STRAP,'Seed',boot_strap_index); 
        
        subject_index_bootstrap=...
            datasample(seed,subject_index_find,n_subject_this_age_group,'Replace',true);       
        
        clear candidate_RDMs
        
        % Candidate Models: Conceptual RDMS (averaged across subjects)
        
        %Conceptual model for this age group (will matter if ACC is included)
        reduced_conceptual_RDMs=list_reduced_conceptual_RDMs{group_index};
        
        for conceptual_index=1:numel(reduced_conceptual_RDMs)
            candidate_RDMs{conceptual_index}=reduced_conceptual_RDMs(conceptual_index);
        end
        
        % Candidate Models: Neural RDMS (averaged across subjects)
        for neural_index=1:numel(LIST_VOI_NAME)
            candidate_RDMs{conceptual_index+neural_index}=list_neural_avg_RDMs{group_index}(neural_index);
        end
        
        for v=1:numel(LIST_VOI_NAME)
            % The reference model (the neural RDMs from every subjects)
            reference_neural_RDMs=neural_RDMs(v,subject_index_bootstrap);
            
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
            
            boot_strap_num=['b' num2str(boot_strap_index)];
            collect_stats_p_r.(valid_name_group).(valid_name_VOI).(boot_strap_num)=...
                stats_p_r;
             sprintf(strcat(valid_name_group,'\n',valid_name_VOI,'\n',boot_strap_num))
        end
    end
end
% end

%Save the relatedness test result
file_test_result_save=...
    fullfile(user_options_2nd_order.rootPath,['test_Relatedness_Test_Results_',VERSION_NUM,'_',date,'.mat']);
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
