%% Part3: DCM - Generate all sorts of models with B matrix nuances (Leave-one-out at one time) [For Bayesian model selection afterwards] 
addpath('D:\GoogleDrive\Lambda_code\m_file\LAMBDA\helper_function');% enable DCM_specify_model()
addpath('D:\GoogleDrive\Lambda_code\m_file\LAMBDA\tool_code');% enable read_mixed_csv()
%% Some constants (research specific)
% To get the subject list with valid runs which has valid first level contrast, i.e., not being excluded)
path_valid_id='D:\Yun-Shiuan_LAMBDA\First_level_parametric_design_matrix_and_result\with_reslice'; 
% Retrive the run list not being excluded due to excess of head motion or
% poor behavioral performance
path_inclusion_index='D:\Yun-Shiuan_LAMBDA\Run_inclusion_info\inclusive_runs_indexes.csv'; 
% Root path for DCM specified and estimated result files
path_DCM_workspace_root='D:\Yun-Shiuan_LAMBDA\DCM\Part2_DCM_Specify_and_Estimate_with_reslice_normalized\full_duration\resliced_peaks\Houde_resliced_peak_4mm_C_matrix_only_V1\B_leave_one_out';
% Root path for 1st level GLM for DCM specification
path_GLM_root='D:\Yun-Shiuan_LAMBDA\DCM\Part0_First_GLM_for_DCM_single_run_with_reslice\full_duration\';
% Root path for VOIs
path_VOIs_root='D:\Yun-Shiuan_LAMBDA\DCM\Part1_VOIs_fixed_ROI_for_DCM_with_reslice\Houde_resliced_peak_4mm\';
% Constant parameters for my DCM specification
interested_inputs=[1 2 3];
TE=0.022;
A_matrix=ones(8,8);
C_matrix=[0 0 0; %L_IPS
          0 0 0; %L_M1
          0 0 0; %L_SMA
          1 1 1; %L_V1
          0 0 0; %R_IPS
          0 0 0; %R_SMA
          1 1 1; %R_V1
          0 0 0;]; %R_DLPFC
      
% A list of B interested matrices
interested_B=[% Same-hemisphere vertical forward flows
              [1,4];[3,1];[2,3];...%L_V1->L_IPS->L_SMA->L_M1
              [5,7];[6,5];...%R_V1->R_IPS->R_SMA
              
              % Horizontal flows
              [7,4];[4,7];... %R_V1 <-> L_V1
              [5,1];[1,5];... %R_IPS <-> L_IPS
              [6,3];[3,6];...%R_SMA <-> L_SMA
              
              [1,8];[2,8];[3,8];[4,8];[5,8];[6,8];[7,8]; % DLPFC to all other regions
              
              % Cross-hemispher vertical forward flows
              [1,7];[5,4]; % R_V1->L_IPS, L_V1->R_IPS
              [3,5];[6,1]; % R_IPS->L_SMA, L_IPS->R_SMA
              [2,6] % R_SMA->L_M1
              
              % Same-hemisphere vertical backward flows
              [4,1];[7,5]; % L_IPS->L_V1, R_IPS->R_V1
              % Cross-hemisphere vertical backward flows
              [7,1];[4,5]; % L_IPS->R_V1, R_IPS->L_V1   

              ];
interested_B_names={% Same-hemisphere vertical forward flows
                    'L_V1-L_IPS';'L_IPS-L_SMA';'L_SMA-L_M1';...
                    'R_V1-R_IPS';'R_IPS-R_SMA';...
                    'L_V1-R_V1';'R_V1-L_V1';...
                    'L_IPS-R_IPS';'R_IPS-L_IPS';...
                    'L_SMA-R_SMA';'R_SMA-L_SMA';...
                    % Horizontal flows
                    'R_DLPFC-L_IPS';'R_DLPFC-L_M1';'R_DLPFC-L_SMA';'R_DLPFC-L_V1';...
                    'R_DLPFC-R_IPS';'R_DLPFC-R_SMA';'R_DLPFC-R_V1';...
                    % Cross-hemispher vertical forward flows
                    'R_V1-L_IPS';'L_V1-R_IPS';...
                    'R_IPS-L_SMA';'L_IPS-R_SMA';...
                    'R_SMA-L_M1';...
                    % Same-hemisphere vertical backward flows
                    'L_IPS-L_V1'; 'R_IPS-R_V1';
                    % Cross-hemisphere vertical backward flows
                    'L_IPS-R_V1'; 'R_IPS_L_V1';
                    };
B_matrices=cell(size(interested_B,1),1); % Collect these B matrices

for i=1:size(interested_B,1)
    B=ones(8,8,3);
    B(interested_B(i,1),interested_B(i,2),:)=0;
    B_matrices{i}=B;
end
%% Get valid IDs
subject_list=cellstr(ls(path_valid_id)); % subject list (which has valid first level contrast, i.e., not being excluded)
subject_list=subject_list(~cellfun(@isempty,regexp(subject_list,'(?<=df)\d+','match')));% include only subject name

%% Read in run inclusion index info
run_inclusion_index=cellfun(@(x) regexprep(x,'"',''),...
    read_mixed_csv(path_inclusion_index,','),'un',0);
run_inclusion_index=table(run_inclusion_index(2:end,2),run_inclusion_index(2:end,3),...
    'VariableNames',{'sub_id','run_num'});

%% DCM specification (per subject run)
for B=1:size(B_matrices,1)
    % Change the B matrix at each model-wise iteration
    B_matrix=B_matrices{B};
    % For this specific model, the path for DCM specified and estimated result files
    path_DCM_workspace_model=fullfile(path_DCM_workspace_root,interested_B_names{B});
    for id=1:size(subject_list,1)
        % Read in the inclusion run indexes and use them to index runs in the following for loop
        run_rows=strcmp(run_inclusion_index.sub_id,subject_list{id});
        run_valid=table2cell(run_inclusion_index(run_rows,'run_num'));
        run_valid=regexp(run_valid,'(?<=r)\d+','match');
        run_valid=str2double([run_valid{:}]');
        output_dir=fullfile(path_DCM_workspace_model,subject_list{id});
        mkdir(output_dir)

        % Create DCM.mat per subject run    
        for run=1:length(run_valid)
            clear DCM
            % The r index refer to the real run number of the task (so it might not be consecutive due to exclusive runs)
            % Could be used when loading the VOI_runx.mat, which named by the real run number
            r=run_valid(run);                 
            file_output=fullfile(output_dir,strcat('DCM_model_specified_run',num2str(r),'.mat'));
            % Skip the iteration if already done.
            if exist(file_output,'file')==2
                warning(['Already done - Model: ',char(num2str(B)),'; Subject: ',char(subject_list{id}),...
                         '; Run: ',char(num2str(r)),'.'])
            else
            %% Parameters for DCM specification (changes dynamically for each run)
            % Load in the GLM.mat per subject run
            GLM_mat=fullfile(path_GLM_root,subject_list{id},strcat('run',num2str(r)),'SPM.mat');
            
            % Load in VOIs per subject run
            VOI_dir=fullfile(path_VOIs_root,subject_list{id},'VOI_data_normalized');
            VOIs=cellstr(ls(strcat(VOI_dir,'\*run',num2str(r),'*mat'))); % List all VOI.mat files
            VOIs=fullfile(VOI_dir,VOIs);
            
            % Specify DCM model
            DCM=DCM_specify_model(GLM_mat,interested_inputs,VOIs,TE,A_matrix,B_matrix,C_matrix);
            strcat(['Done - Model: ',char(num2str(B)),'; Subject: ',char(subject_list{id}),...
                         '; Run: ',char(num2str(r)),'.'])
            save(file_output,'DCM');
            end
        end
    end
end