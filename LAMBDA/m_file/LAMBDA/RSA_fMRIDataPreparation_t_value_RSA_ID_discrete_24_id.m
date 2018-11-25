%% [Used] Part 1 - ii: RSA : Prepare data for RSA usage (before construct RDM)
%% All together (regard distance as discrete - RSA ID : 1~36, repeated across all subject runs)
addpath(genpath('D:\Yun-Shiuan_LAMBDA\rsatoolbox\'));
addpath('D:\GoogleDrive\Lambda_code\m_file\LAMBDA\tool_code');% enable read_mixed_csv()

%% Some constants
% Parameters
list_VOI_name={"R_V1","L_V1",...
    "R_IPS","L_IPS"}';
RSA_IDs_amount=24;
VOI_amount=size(list_VOI_name,1);
%files
file_valid_run='D:\Yun-Shiuan_LAMBDA\Run_inclusion_info\inclusive_runs_indexes.csv';
%paths
path_VOI='D:\Yun-Shiuan_LAMBDA\RSA\Part0-b_VOIs_t_value_extracted_RSA_ID_discrete';
path_first_order_RDM='D:\Yun-Shiuan_LAMBDA\RSA\Part1_first_order_RDM';
% Read in run inclusion index info
run_inclusion_index=cellfun(@(x) regexprep(x,'"',''),...
    read_mixed_csv(file_valid_run,','),'un',0);
run_inclusion_index=table(run_inclusion_index(2:end,2),run_inclusion_index(2:end,3),...
    'VariableNames',{'sub_id','run_num'});
% Derive subjects with valid runs
subject_list=unique(run_inclusion_index.sub_id);

% % Count the run amounts of each subject
% amount_unique_runs=cell(size(subject_list,1),1);
% for id=1:size(subject_list,1)
%     amount_unique_runs{id,1}=length(run_inclusion_index{strcmp(run_inclusion_index.sub_id,subject_list{id}),...
%                                                        'run_num'});
% end
% unique_runs=table(subject_list,cell2mat(amount_unique_runs),...
%                   'VariableNames',{'sub_id','unique_runs'});
%
% run_inclusion_index = join(run_inclusion_index,unique_runs);
%   Value    Count   Percent
%       1        0      0.00%
%       2        0      0.00%
%       3        3      7.50%
%       4        5     12.50%
%       5        9     22.50%
%       6       23     57.50%

%% Pre-preparation (before setting the parameters for the toolbox)--------------------------------------------------
%% Create Beta file for each condition x subject within each ROI
% Note: this part is done by GLM to extract t maps for each unique trial
% See Part0-a and Part0-b for details


%% [Not used] Preparation (Setting the parameters for the toolbox)--------------------------------------------------------------
% for id=1:40%size(subject_list,1)
%     % Read in the inclusion run indexes and use them to index runs in the following for loop
%     run_rows=strcmp(run_inclusion_index.sub_id,subject_list{id});
%     run_valid=table2cell(run_inclusion_index(run_rows,'run_num'));
%     run_valid=regexp(run_valid,'(?<=r)\d+','match');
%     run_valid=str2double([run_valid{:}]');
%
%     % Create RDM.mat per subject run
%     for run=1:length(run_valid)
%         r=run_valid(run);% The index refer to the real  run number of the task (so it might not be consecutive due to exclusive runs)
%
%         % Load in the VOI names of the subject run and arrange rows by RSA_ID
%         path_VOI_sub_run=fullfile(path_VOI,subject_list{id},['run' num2str(r)]);
%         files_VOI=cellstr(ls(path_VOI_sub_run));
%         files_VOI=files_VOI(~ismember(files_VOI,{'.','..'}));
%
%         RSA_ID=regexp(files_VOI,'(?<=RSA_)\d+(?=_(L|R))','match');
%         RSA_ID=str2double([RSA_ID{:}]');
%
%         VOI_name=regexp(files_VOI,'(?<=RSA_\d+_).*(?=\.mat)','match');
%         VOI_name=[VOI_name{:}]';
%
%         table_VOI_files=table(files_VOI,RSA_ID,VOI_name);
%         % arrange rows by RSA_ID
%         table_VOI_files=sortrows(table_VOI_files,'RSA_ID');
%
%         for v=1:size(VOI_name,1)
%             % Beta values of the specific VOI and subject run
%             VOI_files_of_the_VOI=table_VOI_files{strcmp(table_VOI_files.VOI_name,VOI_name{v}),'files_VOI'};
%             RSA_ID_of_the_VOI=table_VOI_files{strcmp(table_VOI_files.VOI_name,VOI_name{v}),'RSA_ID'};
%
%
%             %% fMRIDataPreparation(betaCorrespondence, userOptions)
%             % take fMRI data and load it into the correct format for the rest of the toolbox to use.
%             % userOptions
%             % - userOptions.analysisName : A string which is prepended to the saved files.
%             % - userOptions.rootPath : the path of the root of the to-be-saved files
%             % - userOptions.subjectNames : cell array containing subject names
%             % - userOptions.betaPath : A string which contains the absolute path to the location of the beta images
%             %                          (can contain the following wildcards which would
%             %                          be replaced - subjectName(by userOptions.subjectNames) &
%             %                          betaIdentifier (by file names provided by betaCorrespondence)
%             % - userOptions.conditionLabels : A Cell array containing the names of the conditions in this experiment
%             userOptions.analysisName = 'RSA_test';
%             userOptions.rootPath = path_first_order_RDM; %the path of the root of the to-be-saved files
%             userOptions.subjectNames = subject_list';
%             userOptions.betaPath=fullfile(path_VOI,'[[subjectName]]','[[betaIdentifier]]');
%             userOptions.conditionLabels=strtrim(cellstr(num2str(RSA_ID_of_the_VOI))');
%
%             % betaCorrespondence: an array of beta filenames (nSessionsx n Conditions. e.g.,6 runs x 24 conditions sturct array per person.)
%             % That is, betaCorrespondence_true(nCondition).identifier = the string for the nth Condition file name
%             % e.g.,betaCorrespondence_true(1).identifier = 'Beta_red01_true.mat',
%             % which each stores 343 voxel values within the ROI for the 1st condition for Sub1(i.e.,Sub 1--Beta_red01_true)
%             % This corresponds to "responsePatterns_true.trueSimRoI.subject1(:,1)" [343 voxels x 64 conditions]
%             % Note that this is ture: "isequal(responsePatterns_true.trueSimRoI.subject1(:,1) , reshape(betaImage,343,1))"
%             betaCorrespondence=struct('identifier',{});
%             for id=1:size(subject_list,1)
%                 % Read in the inclusion run indexes and use them to index runs in the following for loop
%                 run_rows=strcmp(run_inclusion_index.sub_id,subject_list{id});
%                 run_valid=table2cell(run_inclusion_index(run_rows,'run_num'));
%                 run_valid=regexp(run_valid,'(?<=r)\d+','match');
%                 run_valid=str2double([run_valid{:}]');
%                 for RSA_id=1:RSA_IDs_amount
%                     for run=1:length(run_valid)
%                         r=run_valid(run);% The index refer to the real  run number of the task (so it might not be consecutive due to exclusive runs)
%                         betaCorrespondence(run,RSA_id).identifier = fullfile(['run' num2str(r)],VOI_files_of_the_VOI{RSA_id});
%                     end
%                 end
%             end
%
%
%     % Output
%     % userOptions.analysisName_ImageData.mat
%     % - conatinas the raw beta images in a structure such that
%     %   fullBrainVols.(subject) is a [nVoxels nConditions nSessions] matrix.
%     fullBrainVols = fMRIDataPreparation(betaCorrespondence, userOptions);
%
%     % responsePatterns stores all the Beta values for every Subject x ROI mask x Conditions x Sessions
%     % responsePatterns.(mask).(Subject)[nMaskedVoxels nConditions nSessions]:
%     % e.g.,Subject1: [343(voxles in the ROI),64(conditions),1(run)]
%     % Note:
%     % responsePatterns_true.trueSimRoI.subject1(:,1) corresponds to 343 voxel
%     % values within the ROI for the 1st condition for Sub1(i.e.,Sub 1--Beta_red01_true)
%
%
%         end
%     end

%% Creat "responsePatterns" for  constructRDMs - a structure such that responsePatterns.(mask).(subject) is a [nMaskedVoxels nConditions nSessions]
% Fil in "responsePatterns.(VOI_name{v}).(subject_list{id})=[nMaskedVoxels nConditions nSessions]"
% e.g., 'Beta_red01_true.mat',
% which each stores 343 voxel values within the ROI for the 1st condition for Sub1(i.e.,Sub 1--Beta_red01_true)
% This corresponds to "responsePatterns_true.trueSimRoI.subject1(:,1)" [343 voxels x 64 conditions]
% Note that this is ture: "isequal(responsePatterns_true.trueSimRoI.subject1(:,1) , reshape(betaImage,343,1))"
% Subject loop
for id=1:size(subject_list,1)
    % Read in the inclusion run indexes and use them to index runs in the following for loop
    run_rows=strcmp(run_inclusion_index.sub_id,subject_list{id});
    run_valid=table2cell(run_inclusion_index(run_rows,'run_num'));
    run_valid=regexp(run_valid,'(?<=r)\d+','match');
    run_valid=str2double([run_valid{:}]');
    
    % Run loop
    for run=1:length(run_valid)
        r=run_valid(run);% The index refer to the real  run number of the task (so it might not be consecutive due to exclusive runs)
        
        % Load in the VOI names of the subject run and arrange rows by RSA_ID
        path_VOI_sub_run=fullfile(path_VOI,subject_list{id},['run' num2str(r)]);
        files_VOI=cellstr(ls(path_VOI_sub_run));
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
            full_name_VOI_files_of_the_VOI=fullfile(path_VOI,subject_list{id},['run',num2str(r)],VOI_files_of_the_VOI);
            
            %Initialize the empty array for "[nMaskedVoxels nConditions nSessions]"
            if run==1
                responsePatterns.(char(list_VOI_name{v})).(subject_list{id})=[];
            end
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
            % Concatenate the 2d nMaskedVoxels x nConditions layer in the 3rd dimension (runs)
            responsePatterns.(char(list_VOI_name{v})).(subject_list{id})=...
                cat(3,responsePatterns.(char(list_VOI_name{v})).(subject_list{id}),new_voxel_by_cond_layer);
            
        end
        strcat('Done: sun_id: ', subject_list{id},' ; run: ',num2str(r))
    end
end

%% Creat "betaCorrespondence " - a structure such that responsePatterns.(mask).(subject) is a [nMaskedVoxels nConditions nSessions]
betaCorrespondence=struct('identifier',{});
for id=1:size(subject_list,1)
     % Read in the inclusion run indexes and use them to index runs in the following for loop
    run_rows=strcmp(run_inclusion_index.sub_id,subject_list{id});
    run_valid=table2cell(run_inclusion_index(run_rows,'run_num'));
    run_valid=regexp(run_valid,'(?<=r)\d+','match');
    run_valid=str2double([run_valid{:}]');
    for RSA_id=1:RSA_IDs_amount
        for run=1:length(run_valid)
            betaCorrespondence(run,RSA_id).identifier = VOI_files_of_the_VOI;
        end
    end
end

% Name the ROI for both streams of data
RoIName = 'SimROI';
responsePatterns.(['true' RoIName]) = fullBrainVols_true;

% the name of the root of the to-be-saved files
userOptions.analysisName
% the path of the root of the to-be-saved files
userOptions.rootPath

% cell array containing subject names
userOptions.subjectNames

% The distance measure for RDM calculation (default to correlation)
userOptions.distance

%[R G B] value of the colour which should be used to indicated RoI RDMs on various diagrams
userOptions.RoIColor
%% Construct the tool-box specific format from fMRI data
fullBrainVols_true = fMRIDataPreparation(betaCorrespondence_true, userOptions_true);
