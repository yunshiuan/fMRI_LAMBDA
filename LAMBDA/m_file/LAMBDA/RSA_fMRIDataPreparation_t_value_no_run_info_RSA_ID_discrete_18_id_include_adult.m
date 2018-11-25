%% [Used] Part 1 - iv-3-c: RSA : Prepare data for RSA usage (before construct RDM) [No run info]
% trial 7 in the README.txt.
% Part1: To generate both neural and conceptual model RDMs for 1st and 2nd oder
% analysis afterwards.

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

% 3: [No runs info]
% Since betas from different runs have already been collapsed into single
% t maps, there is only one t map for each condition per subject VOI.

% c: Include adult neural data (both collapse split by age group)
% Split into 2nd grader, 5th grader, and the adult

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%IPS adjusted by emperical results (2018/06/10)
%Update subject list for children (i.e., inclusive_runs_indexes_new_June_10.csv) (2018/6/10)
%Update subject list for children (i.e., inclusive_runs_indexes_new_June_11.csv) (2018/6/10)

%% All together (regard distance as discrete - RSA ID : 1~18, repeated across all subject runs)
PATH_ROOT='D:\Yun-Shiuan_LAMBDA';
PATH_TOOL_CODE='D:\GoogleDrive\Lambda_code\m_file\LAMBDA\tool_code';
addpath(genpath(fullfile(PATH_ROOT,'rsatoolbox')));
addpath(PATH_TOOL_CODE);% enable read_mixed_csv()

%% Some constants
% Parameters
LIST_VOI_NAME={"R_V1","L_V1",...
    "R_IPS","L_IPS","R_DLPFC"}';
% The mapping information used for renaming the conceptual models
LIST_CONCEPTUAL_MODEL_OLD_TERMS=... 
    cellstr({
    'RDMdist',...
    'allthesame','formatwithoutlocation','sameformatpair',...%Pure format models
    'signeddistance','absdistance',...%Linear Distance
    'signeddist','absdist','discrete','nodistance','only','LL',...    
    '2grade','5grade','adult','diff','acc',...%Difficulty
    'log' %Log transform
    })';
LIST_CONCEPTUAL_MODEL_NEW_TERMS=...
    cellstr({'','f','g','n',...%Pure format models
    's','a',...%Linear Distance
    's','a','','','','l',...
    'S','F','A','','ac',...%Difficulty,
    'g'%Difficulty
    })';
TABLE_CONCEPTUAL_MODEL_OLD2NEW_MAP=...
    table(LIST_CONCEPTUAL_MODEL_OLD_TERMS,LIST_CONCEPTUAL_MODEL_NEW_TERMS);

NUMBER_VOI=numel(LIST_VOI_NAME);
RSA_IDS_AMOUNT=18;
VERSION_NUM='22'; %The version number

%files
FILE_VALID_RUN_CHILD=fullfile(PATH_ROOT,'Run_inclusion_info','inclusive_runs_indexes_new_June_11.csv');
FILE_VALID_RUN_ADULT=fullfile(PATH_ROOT,'Adult','Run_inclusion_info','inclusive_runs_indexes.csv');
FILE_DEMOGRAPHIC_CHILD=fullfile(PATH_ROOT,'demographic','tidy_demographic.csv');
FILE_DEMOGRAPHIC_ADULT=fullfile(PATH_ROOT,'Adult','demographic','tidy_demographic.csv');
%paths
PATH_THIS_TRIAL=fullfile(PATH_ROOT,'RSA',['trial_',num2str(VERSION_NUM)]);
PATH_VOI_CHILD=fullfile(PATH_ROOT,'RSA','trial_21','Part0-b_VOIs_t_value_extracted','child'); % Reuse
PATH_VOI_ADULT=fullfile(PATH_ROOT,'RSA','trial_21','Part0-b_VOIs_t_value_extracted','adult'); % Reuse
PATH_RDMS=fullfile(PATH_THIS_TRIAL,'Part1_neural_and_conceptual_RDMs');
PATH_CONCEPTUAL_MODEL=fullfile(PATH_ROOT,'RSA','trial_19',...
    'Part1_neural_and_conceptual_RDMs','conceptual_model_RDMs');
% Read in run inclusion index info 
%Child
run_inclusion_child_index=read_mixed_csv_to_table(FILE_VALID_RUN_CHILD);
subject_list_child=unique(run_inclusion_child_index.sub_id);% Derive subjects with valid runs
%Adult
run_inclusion_adult_index=read_mixed_csv_to_table(FILE_VALID_RUN_ADULT);
subject_list_adult=unique(run_inclusion_adult_index.sub_id);% Derive subjects with valid runs
% subject_list_adult=subject_list_adult(~ismember(subject_list_adult,{'XFC305'}));
% %305 have missing run (this has been solved)

% Read in demographic data (for age group split)------------------------
%Child
table_demographic_child=read_mixed_csv_to_table(FILE_DEMOGRAPHIC_CHILD);
table_demographic_child=table_demographic_child(:,{'sub_id','grade'});
table_demographic_child.sub_id=lower(table_demographic_child.sub_id);
table_demographic_child=join(table(subject_list_child,'VariableNames',{'sub_id'}),...
    table_demographic_child);
%Adult
table_demographic_adult=read_mixed_csv_to_table(FILE_DEMOGRAPHIC_ADULT);
table_demographic_adult=table_demographic_adult(:,{'sub_id','grade'});
table_demographic_adult=join(table(subject_list_adult,'VariableNames',{'sub_id'}),...
    table_demographic_adult);
%Join the Two demographic Table (Child+Adult)
table_demographic_all=vertcat(table_demographic_child,table_demographic_adult);
index_2_grade=strcmp(table_demographic_all.grade,'2');
index_5_grade=strcmp(table_demographic_all.grade,'5');
index_child=logical(index_2_grade+index_5_grade);%All children including 2 and 5 graders
index_adult=strcmp(table_demographic_all.grade,'Adult');
subject_list=table_demographic_all.sub_id;

%% (Step1) Pre-preparation (before setting the parameters for the toolbox)--------------------------------------------------
%% Create Beta file for each condition x subject within each ROI
% Note: this part is done by GLM to extract t maps for each unique trial
% See Part0-a and Part0-b for details
%% (Step2) constructRDMs(responsePatterns, betaCorrespondence, userOptions) - Construct the neural RDM in the tool-box specific format from fMRI data
%% (Step2-a) Creat "responsePatterns" for  constructRDMs - a structure such that responsePatterns.(mask).(subject) is a [nMaskedVoxels nConditions nSessions]
% Fil in "responsePatterns.(VOI_name{v}).(subject_list{id})=[nMaskedVoxels nConditions nSessions]"
% e.g., 'Beta_red01_true.mat',
% which each stores 343 voxel values within the ROI for the 1st condition for Sub1(i.e.,Sub 1--Beta_red01_true)
% This corresponds to "responsePatterns_true.trueSimRoI.subject1(:,1)" [343 voxels x 64 conditions]
% Note that this is true: "isequal(responsePatterns_true.trueSimRoI.subject1(:,1) , reshape(betaImage,343,1))"
% Subject loop
for id=1:numel(subject_list)    
    % Load in the VOI names of the subject run and arrange rows by RSA_ID
    if ismember(subject_list{id},subject_list_child)
        PATH_VOI=PATH_VOI_CHILD;
    elseif ismember(subject_list{id},subject_list_adult)
        PATH_VOI=PATH_VOI_ADULT;
    end
    
    path_VOI_sub=fullfile(PATH_VOI,subject_list{id});
    files_VOI=dir(fullfile(path_VOI_sub,'t_map*')); %Use dir instead of ls so that it works both on linux and other OS
    files_VOI={files_VOI.name}';
    
    RSA_ID=regexp(files_VOI,'(?<=RSA_)\d+(?=_(L|R))','match');
    RSA_ID=str2double([RSA_ID{:}]');
    
    VOI_name=regexp(files_VOI,'(?<=RSA_\d+_).*(?=\.mat)','match');
    VOI_name=[VOI_name{:}]';
    
    table_VOI_files=table(files_VOI,RSA_ID,VOI_name);
    % arrange rows by RSA_ID
    table_VOI_files=sortrows(table_VOI_files,'RSA_ID');
    % VOI loop
    for v=1:NUMBER_VOI
        % t values of the specific VOI and subject 
        VOI_files_of_the_VOI=table_VOI_files{strcmp(table_VOI_files.VOI_name,LIST_VOI_NAME{v}),'files_VOI'};
        full_name_VOI_files_of_the_VOI=fullfile(PATH_VOI,subject_list{id},VOI_files_of_the_VOI);
        
        %Initialize the empty array for "[nMaskedVoxels nConditions nSessions]"
        responsePatterns.(char(LIST_VOI_NAME{v})).(subject_list{id})=[];
        
        % For collecting a single 2d nMaskedVoxels x nConditions layer
        new_voxel_by_cond_layer=[];
        
        % RSA ID loop
        for RSA_id=1:RSA_IDS_AMOUNT
            %Load in the beta image
            beta_image=load(full_name_VOI_files_of_the_VOI{RSA_id});
            
            % responsePatterns.(list_VOI_name{v}).(subject_list{id})=[nMaskedVoxels nConditions nSessions];
            %                                                   [231 x 24 x 6]
            new_voxel_by_cond_layer=...
                [new_voxel_by_cond_layer, beta_image.t_map_in_ROI']; %#ok<AGROW>
        end
        % Concatenate the 2d nMaskedVoxels x nConditions layer in the 3rd dimension ( if runs exist)
        responsePatterns.(char(LIST_VOI_NAME{v})).(subject_list{id})=...
            cat(3,responsePatterns.(char(LIST_VOI_NAME{v})).(subject_list{id}),new_voxel_by_cond_layer);
        
    end
    strcat('Done: sub_id: ', subject_list{id}) 
end

%% (Step2-b) Creat "betaCorrespondence" - the size of  [nConditions nSessions]
% Note that the only information extracted in the constructRDMs()is the size of betaCorrespondece
% The size of betaCorrespondence is 
% 	nSessions = size(betas, 1);
% 	nConditions = size(betas, 2);
% Hence, assign betaCorrespondece as any matrix with the correct size is
% fine.
betaCorrespondence = zeros(1,RSA_IDS_AMOUNT);

%% (Step2-c) Create "userOptions"
% userOptions.analysisName: A string which is prepended to the saved files.
% userOptions.rootPath: A string describing the root path where files will be saved (inside created directories).
% userOptions.maskNames: A cell array containing strings identifying the mask names. It defaults to the field names of the first subject of responsePatterns.
% userOptions.subjectNames: A cell array containing strings identifying the subject names. Defaults to the field names in (the first mask of) responsePatterns.
% userOptions.distance: A string indicating the distance measure with which to calculate the RDMs. Defaults to �Correlation�, but can be set to any Matlab distance measure.
% userOptions.RoIColor: A triple indicating the [R G B] value of the colour which should be used to indicated RoI RDMs on various diagrams. Defaults to black ([0 0 0]).
clear userOptions;
userOptions.analysisName = ['RSA_trial_',VERSION_NUM,'_',date];
userOptions.rootPath = PATH_RDMS; %the path of the root of the to-be-saved files
mkdir(userOptions.rootPath);

%% constructRDMs(responsePatterns, betaCorrespondence, userOptions) - Construct the tool-box specific format from fMRI data
neural_RDMs = constructRDMs(responsePatterns, betaCorrespondence, userOptions);  %Output: RDM_test [nMasks, nSubjects, nSessions]
%(Averaging the sessions has no effect on the RDMs, since t map collapsed run info.
% This is just used to get read of the session labels in the RDMs' names.)
neural_RDMs = averageRDMs_subjectSession(neural_RDMs, 'session'); %Output: RDM [nMasks, 1] 
neural_RDMs_avg_2_grade=averageRDMs_subjectSession(neural_RDMs(:,index_2_grade), 'subject'); %Output: RDM [nMasks, 1]
neural_RDMs_avg_5_grade=averageRDMs_subjectSession(neural_RDMs(:,index_5_grade), 'subject'); %Output: RDM [nMasks, 1]
neural_RDMs_avg_child=averageRDMs_subjectSession(neural_RDMs(:,index_child), 'subject'); %Output: RDM [nMasks, 1]
neural_RDMs_avg_adult=averageRDMs_subjectSession(neural_RDMs(:,index_adult), 'subject'); %Output: RDM [nMasks, 1]

neural_RDMs_avg_all=averageRDMs_subjectSession(neural_RDMs, 'subject'); %Output: RDM [nMasks, 1]

%% (Step4) constructModelRDMs(rawModels, userOptions)- Construct conceptual RDMs
% rawModels:  a structure in which rawModels.(modelName) is the model RDM [24x24 double]
% userOptions.analysisName: A string which is prepended to the saved files.
% userOptions.rootPath: A string describing the root path where files will be saved (inside created directories).
% userOptions.ModelColor: A triple indicating the [R G B] value of the colour which should be used to indicated model RDMs on various diagrams. Defaults to black ([0 0 0]).

%% (Step4-a) Create the conceptual models: "rawModels"
% List all RDMs I want to match with the neural RDM
% from csv file names
list_conceptual_RDMs_csv=dir(fullfile(PATH_CONCEPTUAL_MODEL,'RDM*'));
list_conceptual_RDMs_csv={list_conceptual_RDMs_csv.name}';
list_conceptual_RDMs_csv=list_conceptual_RDMs_csv(~cellfun(@isempty,regexp(list_conceptual_RDMs_csv,'^RDM.*','match')));
list_conceptual_RDMs_csv=fullfile(PATH_CONCEPTUAL_MODEL,list_conceptual_RDMs_csv);

% RDM names
list_conceptual_RDMs_name=regexp(list_conceptual_RDMs_csv,'(?<=RDM_dist_).*(?=.csv)','match');
list_conceptual_RDMs_name=[list_conceptual_RDMs_name{:}]';
list_conceptual_RDMs_name=regexprep(list_conceptual_RDMs_name,'0.','');

clear raw_conceptual_models
for RDM_index=1:numel(list_conceptual_RDMs_csv)
        RDM=csvread(list_conceptual_RDMs_csv{RDM_index});
        raw_conceptual_models.(char(list_conceptual_RDMs_name{RDM_index})) = RDM;
end

%% (Step3-b) Create the "userOptions"
% Note that userOptions.analysisName & userOptions.rootPath have already
% been assigned in the previous step.
userOptions.ModelColor= [0 1 0];
%% constructModelRDMs(rawModels, userOptions)- Construct conceptual RDMs
conceptual_RDMs=constructModelRDMs(raw_conceptual_models, userOptions);
user_options_RDM_construction=userOptions;
%% Save RDMs (neural RDMs, conceptual RDMs, run-averaged conceptual RDMs
% Rename(shorten) the subject id for the sake of visualization
for VOI_index=1:size(neural_RDMs,1)
    for id=1:size(neural_RDMs,2)
    neural_RDMs(VOI_index,id).name=...
        regexprep(neural_RDMs(VOI_index,id).name,'(df10)|(\|)|(\s)|(XFC)','');
    end
end

% Rename(shorten) the conceptual RDMs for the sake of visualization
old_name=conceptual_RDMs;
new_name=regexprep({old_name.name}','\s','');

for term_index=1:size(TABLE_CONCEPTUAL_MODEL_OLD2NEW_MAP,1)
    new_name = ...
        regexprep(new_name,...
        TABLE_CONCEPTUAL_MODEL_OLD2NEW_MAP.LIST_CONCEPTUAL_MODEL_OLD_TERMS(term_index),...
        TABLE_CONCEPTUAL_MODEL_OLD2NEW_MAP.LIST_CONCEPTUAL_MODEL_NEW_TERMS(term_index));
end
    % check: [{old_name.name}',new_name]
for model_index=1:size(conceptual_RDMs,2)    
    conceptual_RDMs(model_index).name=...
        new_name{model_index};
end

path_RDM_output=fullfile(userOptions.rootPath,'RDMs',['RDMs_and_options_neural_and_conceptual_trial_',VERSION_NUM,'_',date,'.mat']);
save(path_RDM_output,...
     'neural_RDMs',...
     'neural_RDMs_avg_2_grade','neural_RDMs_avg_5_grade',...
     'neural_RDMs_avg_child','neural_RDMs_avg_adult',...
     'neural_RDMs_avg_all','conceptual_RDMs',...
     'user_options_RDM_construction')