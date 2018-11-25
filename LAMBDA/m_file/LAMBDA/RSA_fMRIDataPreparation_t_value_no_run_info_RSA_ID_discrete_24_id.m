%% [Used] Part 1 - ii-3: RSA : Prepare data for RSA usage (before construct RDM) [No run info]
% Part1: To generate both neural and conceptual model RDMs for 1st and 2nd oder
% analysis afterwards.
% ii: use t maps (24 unique trials using discrete distance) 
% instead of beta maps (144  unique trials with continuous distance)
% 3: [No runs info]
% Since betas from different runs have already been collapsed into single
% t maps, there is only one t map for each condition per subject VOI.

%% All together (regard distance as discrete - RSA ID : 1~24, repeated across all subject runs)
addpath(genpath('D:\Yun-Shiuan_LAMBDA\rsatoolbox\'));
addpath('D:\GoogleDrive\Lambda_code\m_file\LAMBDA\tool_code');% enable read_mixed_csv()

%% Some constants
% Parameters
list_VOI_name={"R_V1","L_V1",...
               "R_IPS","L_IPS"}';
VOI_amount=size(list_VOI_name,1);
RSA_IDs_amount=24;
version_num='4'; %The version number
%files
file_valid_run='D:\Yun-Shiuan_LAMBDA\Run_inclusion_info\inclusive_runs_indexes.csv';
%paths
path_VOI='D:\Yun-Shiuan_LAMBDA\RSA\Part0-b_VOIs_t_value_extracted_RSA_ID_discrete_repeated_24';
path_RDMs='D:\Yun-Shiuan_LAMBDA\RSA\Part1_neural_and_conceptual_RDMs_RSA_ID_discrete_repeated_24';
path_conceptual_model=fullfile(path_RDMs,'conceptual_model_RDMs');
% Read in run inclusion index info (note that only subject matters for
% t-map based analysis)
run_inclusion_index=cellfun(@(x) regexprep(x,'"',''),...
    read_mixed_csv(file_valid_run,','),'un',0);
run_inclusion_index=table(run_inclusion_index(2:end,2),run_inclusion_index(2:end,3),...
    'VariableNames',{'sub_id','run_num'});

% Derive subjects with valid runs
subject_list=unique(run_inclusion_index.sub_id);


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
for id=1:size(subject_list,1)    
    % Load in the VOI names of the subject run and arrange rows by RSA_ID
    path_VOI_sub=fullfile(path_VOI,subject_list{id});
    files_VOI=cellstr(ls(path_VOI_sub));
    files_VOI=files_VOI(~ismember(files_VOI,{'.','..'}));
    
    RSA_ID=regexp(files_VOI,'(?<=RSA_)\d+(?=_(L|R))','match');
    RSA_ID=str2double([RSA_ID{:}]');
    
    VOI_name=regexp(files_VOI,'(?<=RSA_\d+_).*(?=\.mat)','match');
    VOI_name=[VOI_name{:}]';
    
    table_VOI_files=table(files_VOI,RSA_ID,VOI_name);
    % arrange rows by RSA_ID
    table_VOI_files=sortrows(table_VOI_files,'RSA_ID');
    % VOI loop
    for v=1:VOI_amount
        % Beta values of the specific VOI and subject run
        VOI_files_of_the_VOI=table_VOI_files{strcmp(table_VOI_files.VOI_name,list_VOI_name{v}),'files_VOI'};
        full_name_VOI_files_of_the_VOI=fullfile(path_VOI,subject_list{id},VOI_files_of_the_VOI);
        
        %Initialize the empty array for "[nMaskedVoxels nConditions nSessions]"
        responsePatterns.(char(list_VOI_name{v})).(subject_list{id})=[];
        
        % For collecting a single 2d nMaskedVoxels x nConditions layer
        new_voxel_by_cond_layer=[];
        
        % RSA ID loop
        for RSA_id=1:RSA_IDs_amount
            %Load in the beta image
            beta_image=load(full_name_VOI_files_of_the_VOI{RSA_id});
            
            % responsePatterns.(list_VOI_name{v}).(subject_list{id})=[nMaskedVoxels nConditions nSessions];
            %                                                   [231 x 24 x 6]
            new_voxel_by_cond_layer=...
                [new_voxel_by_cond_layer, beta_image.t_map_in_ROI']; %#ok<AGROW>
        end
        % Concatenate the 2d nMaskedVoxels x nConditions layer in the 3rd dimension ( if runs exist)
        responsePatterns.(char(list_VOI_name{v})).(subject_list{id})=...
            cat(3,responsePatterns.(char(list_VOI_name{v})).(subject_list{id}),new_voxel_by_cond_layer);
        
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
betaCorrespondence = zeros(1,RSA_IDs_amount);

%% (Step2-c) Create "userOptions"
% userOptions.analysisName: A string which is prepended to the saved files.
% userOptions.rootPath: A string describing the root path where files will be saved (inside created directories).
% userOptions.maskNames: A cell array containing strings identifying the mask names. It defaults to the field names of the first subject of responsePatterns.
% userOptions.subjectNames: A cell array containing strings identifying the subject names. Defaults to the field names in (the first mask of) responsePatterns.
% userOptions.distance: A string indicating the distance measure with which to calculate the RDMs. Defaults to “Correlation”, but can be set to any Matlab distance measure.
% userOptions.RoIColor: A triple indicating the [R G B] value of the colour which should be used to indicated RoI RDMs on various diagrams. Defaults to black ([0 0 0]).
clear userOptions;
userOptions.analysisName = ['RSA_trial_',version_num,'_',date];
userOptions.rootPath = path_RDMs; %the path of the root of the to-be-saved files


%% constructRDMs(responsePatterns, betaCorrespondence, userOptions) - Construct the tool-box specific format from fMRI data
neural_RDMs = constructRDMs(responsePatterns, betaCorrespondence, userOptions);  %Output: RDM_test [nMasks, nSubjects, nSessions]
%(Averaging the sessions has no effect on the RDMs, since t map collapsed run info.
% This is just used to get read of the session labels in the RDMs' names.)
neural_RDMs = averageRDMs_subjectSession(neural_RDMs, 'session'); %Output: RDM [nMasks, 1] 
neural_RDMs_avg_subjects=averageRDMs_subjectSession(neural_RDMs, 'subject'); %Output: RDM [nMasks, 1]
%% (Step4) constructModelRDMs(rawModels, userOptions)- Construct conceptual RDMs
% rawModels:  a structure in which rawModels.(modelName) is the model RDM [24x24 double]
% userOptions.analysisName: A string which is prepended to the saved files.
% userOptions.rootPath: A string describing the root path where files will be saved (inside created directories).
% userOptions.ModelColor: A triple indicating the [R G B] value of the colour which should be used to indicated model RDMs on various diagrams. Defaults to black ([0 0 0]).

%% (Step4-a) Create the "rawModels"
% List all RDMs I want to match with the neural RDM
% csv file names
list_conceptual_RDMs_csv=cellstr(ls(path_conceptual_model));
list_conceptual_RDMs_csv=list_conceptual_RDMs_csv(~cellfun(@isempty,regexp(list_conceptual_RDMs_csv,'^RDM.*','match')));
list_conceptual_RDMs_csv=fullfile(path_conceptual_model,list_conceptual_RDMs_csv);

% RDM names
list_concpetual_RDMs_name=regexp(list_conceptual_RDMs_csv,'RDM_dist.*(?=.csv)','match');
list_concpetual_RDMs_name=[list_concpetual_RDMs_name{:}]';
list_concpetual_RDMs_name=regexprep(list_concpetual_RDMs_name,'0.','');

clear raw_conceptual_models
for RDM_index=1:size(list_conceptual_RDMs_csv,1)
    RDM=csvread(list_conceptual_RDMs_csv{RDM_index});
    raw_conceptual_models.(char(list_concpetual_RDMs_name{RDM_index})) = RDM;
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
    for(id=1:size(neural_RDMs,2))
    neural_RDMs(VOI_index,id).name=...
        regexprep(neural_RDMs(VOI_index,id).name,'(df10)|(\|)|(\s)','');
    end
end
% Rename(shorten) the conceptual RDMs for the sake of visualization
old_name=conceptual_RDMs;
new_name=regexprep({old_name.name}','\s','');
% format abbreviation: fnl/fl/nf(format without location/ format with location/ no format)
format_without_location_index=~cellfun(@isempty,regexp(new_name,'formatwithoutlocation','match'))*1;
format_with_location_index=~cellfun(@isempty,regexp(new_name,'formatwithlocation','match'))*2;
no_format_index=cellfun(@isempty,regexp(new_name,'format','match'))*3;

format_abv_index=format_with_location_index+format_without_location_index+no_format_index;
format_abv=cellstr(["fnl","fl","nf"]);
format_abv=format_abv(format_abv_index);
% distance abbreviation: abd/sd/_(absolute distance/ signed distanct/ no distance)
abs_dist_index=~cellfun(@isempty,regexp(new_name,'absdist','match'))*1;
signed_dist_index=~cellfun(@isempty,regexp(new_name,'signeddist','match'))*2;
no_dist_index=~cellfun(@isempty,regexp(new_name,'nodist','match'))*3;

dist_abv_index=abs_dist_index+signed_dist_index+no_dist_index;
dist_abv=cellstr(["abd","sd",""]);
dist_abv=dist_abv(dist_abv_index);
% weigthts
weight_abv=regexp({old_name.name}','\d$','match');
empties=cellfun('isempty',weight_abv);
weight_abv(empties)={""};
weight_abv=[weight_abv{:}];
% Paste0 all abrreviated sub-strings together
strcat(format_abv,dist_abv,weight_abv)
% check: [new_name,strcat(format_abv,dist_abv,weight_abv)']
new_name=strcat(format_abv,dist_abv,weight_abv);
for model_index=1:size(conceptual_RDMs,2)    
    conceptual_RDMs(model_index).name=...
        new_name{model_index};
end

path_RDM_output=fullfile(userOptions.rootPath,'RDMs',['RDMs_and_options_neural_and_conceptual_trial_',version_num,'_',date,'.mat']);
save(path_RDM_output,...
     'neural_RDMs','neural_RDMs_avg_subjects','conceptual_RDMs',...
     'user_options_RDM_construction')