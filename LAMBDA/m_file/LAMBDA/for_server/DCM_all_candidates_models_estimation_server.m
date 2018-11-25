%% Part4: DCM - Estimate all sorts of models with B matrix nuances (Leave-one-out at one time) [For Bayesian model selection afterwards] 
addpath('/study3/devfracs/DCM_YunShiuan/m_script/tool_code');% enable read_mixed_csv()
%% Some constants (research specific)
%Path
% Root path for DCM specified and estimated result files
PATH_DCM_WORKSPACE_ROOT='/study3/devfracs/DCM_YunShiuan/DCM/Part2_DCM_Specify_and_Estimate_with_reslice_normalized/full_duration/resliced_peaks/Houde_resliced_peak_4mm_C_matrix_only_V1/B_leave_one_out';
% path_DCM_workspace_root='D:\Yun-Shiuan_LAMBDA\DCM\Part2_DCM_Specify_and_Estimate_with_reslice_normalized_no_SMA\full_duration\resliced_peaks\Houde_resliced_peak_4mm_C_matrix_only_V1';

% File
% Retrive the run list not being excluded due to excess of head motion or
% poor behavioral performance
FILE_VALID_RUN_CHILD='/study3/devfracs/DCM_YunShiuan/Run_inclusion_infoinclusive_runs_indexes.csv'; 

% A list of B interested matrix names
INTERESTED_B_NAMES={% Same-hemisphere vertical forward flows
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
% For only focusing on the fully specified model                
% interested_B_names={'B_full_specified'};                
%% Get valid IDs & Read in run inclusion index info
run_inclusion_child_index=read_mixed_csv_to_table(FILE_VALID_RUN_CHILD);
subject_list_child=unique(run_inclusion_child_index.sub_id);% Derive subjects with valid runs

%% DCM estimation
for B=27%1:size(interested_B_names,1)
    % For this specific model, the path for DCM specified and estimated result files
    path_DCM_workspace_model=fullfile(PATH_DCM_WORKSPACE_ROOT,INTERESTED_B_NAMES{B});

    for id=1:size(subject_list,1)
        % Read in the inclusion run indexes and use them to index runs in the following for loop
        run_rows=strcmp(run_inclusion_index.sub_id,subject_list{id});
        run_valid=table2cell(run_inclusion_index(run_rows,'run_num'));
        run_valid=regexp(run_valid,'(?<=r)\d+','match');
        run_valid=str2double([run_valid{:}]');
        output_dir=fullfile(path_DCM_workspace_model,subject_list{id});
        cd(output_dir);
        % Create DCM.mat per subject run    
        for run=1:length(run_valid)
            r=run_valid(run);% The index refer to the real run number of the task (so it might not be consecutive due to exclusive runs)
            clear matlabbatch
            % Load in the GLM.matper subject run
            DCM_mat=fullfile(path_DCM_workspace_model,...
                subject_list{id},strcat('DCM_model_specified_run',num2str(r),'.mat'));

            % The name of the DCM result file
            estimated_mat=fullfile(path_DCM_workspace_model,...
                 subject_list{id},strcat('DCM_model_result_estimated_run',num2str(r),'.mat'));
            % Prevent runing if DCM analysis on the run was already done
            if exist(estimated_mat, 'file')==2
                strcat('Already done: Model:',num2str(B),';ID:',subject_list{id},';run',num2str(r))
            else
            % Run the DCM (Should first convert to matlabbatch)
            matlabbatch{1}.spm.dcm.estimate.dcms.subj.dcmmat = {DCM_mat};
            matlabbatch{1}.spm.dcm.estimate.output.separate = struct([]);
            matlabbatch{1}.spm.dcm.estimate.est_type = 1;
            matlabbatch{1}.spm.dcm.estimate.fmri.analysis = 'time';
            spm_jobman('run',matlabbatch);

            % Rename the result DCM and move to the right place
            % (otherwise the estimated results will override the previous DCM.mat)
            old_name=ls('DCM-*Mar*');
            new_name=estimated_mat;
            movefile(old_name,new_name);

            strcat('Finish: Model',num2str(B),'-',INTERESTED_B_NAMES{B},'--',...
                subject_list{id},';run',num2str(r))
            end
        end
    end
end
