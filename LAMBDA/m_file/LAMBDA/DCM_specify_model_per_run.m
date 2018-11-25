%% Part2-a: DCM - DCM specification and estimation (per subject run) (for distance effect DCM)
addpath(char("D:\GoogleDrive\Lambda_code\m_file\LAMBDA\tool_code"));% enable read_mixed_csv()
%spm_jobman('initcfg') -- before running the batch

path='D:\Yun-Shiuan_LAMBDA\First_level_parametric_design_matrix_and_result\with_reslice'; % To get the subject list with valid runs
cd(path);
subject_list=cellstr(ls('*')); % subject list (which has valid first level contrast, i.e., not being excluded)
subject_list=subject_list(~cellfun(@isempty,regexp(subject_list,'(?<=df)\d+','match')));% include only subject name

%% Read in run inclusion index info
run_inclusion_index=cellfun(@(x) regexprep(x,'"',''),...
    read_mixed_csv("D:\Yun-Shiuan_LAMBDA\Run_inclusion_info\inclusive_runs_indexes.csv",','),'un',0);
run_inclusion_index=table(run_inclusion_index(2:end,2),run_inclusion_index(2:end,3),...
    'VariableNames',{'sub_id','run_num'});

%% DCM specification (per subject run)
part2_dir='D:\Yun-Shiuan_LAMBDA\DCM\Part2_DCM_Specify_and_Estimate_with_reslice_normalized\full_duration\resliced_peaks\Houde_resliced_peak_4mm_C_matrix_only_V1';
for id=1:size(subject_list,1)
    % Read in the inclusion run indexes and use them to index runs in the following for loop
    run_rows=strcmp(run_inclusion_index.sub_id,subject_list{id});
    run_valid=table2cell(run_inclusion_index(run_rows,'run_num'));
    run_valid=regexp(run_valid,'(?<=r)\d+','match');
    run_valid=str2double([run_valid{:}]');
    output_dir=fullfile(part2_dir,subject_list{id});
    mkdir(output_dir)
      
    % Create DCM.mat per subject run    
    for run=1:length(run_valid)
        clear DCM
        r=run_valid(run);% The index refer to the real run number of the task (so it might not be consecutive due to exclusive runs)
        % Load in the GLM.matper subject run
        GLM_mat=fullfile('D:\Yun-Shiuan_LAMBDA\DCM\Part0_First_GLM_for_DCM_single_run_with_reslice\full_duration\',...
            subject_list{id},strcat('run',num2str(r)),'SPM.mat');
        load(GLM_mat);
        % The index refer to the real run number of the task (so it might not be consecutive due to exclusive runs)
        % Could be used when loading the VOI_runx.mat, which nemaed by the real
        % run number
       
        % Load regions of interest
        %--------------------------------------------------------------------------
        VOI_dir=fullfile('D:\Yun-Shiuan_LAMBDA\DCM\Part1_VOIs_fixed_ROI_for_DCM_with_reslice\Houde_resliced_peak_4mm\',...
            subject_list{id},'VOI_data_normalized');
        cd(VOI_dir);
        VOIs=cellstr(ls(strcat('*run',num2str(r),'*mat'))); % List all VOI.mat files
        
        for v=1:size(VOIs,1)
            load(fullfile(VOI_dir,VOIs{v}),'xY');
            DCM.xY(v) = xY;
        end
        
        DCM.n = length(DCM.xY);      % number of regions
        DCM.v = length(DCM.xY(1).u); % number of time points
        
        % Time series
        %--------------------------------------------------------------------------
        DCM.Y.dt  = SPM.xY.RT;% TR of the EPI run (2s)
        DCM.Y.X0  = DCM.xY(1).X0; % time-series of the covariates of the GLM (nVolumes x nConfounds))
        for i = 1:DCM.n % Loop through the 8 VOIs
            DCM.Y.y(:,i)  = DCM.xY(i).u; %time course of ROI i (105 x 1)
            DCM.Y.name{i} = DCM.xY(i).name; %name of ROI i
        end
        DCM.Y.Q    = spm_Ce(ones(1,DCM.n)*DCM.v); %Error covariance constraints (cell 1 x nROIs)
        
        % Experimental inputs (per run)
        %--------------------------------------------------------------------------
        DCM.U.dt   =  SPM.Sess.U(1).dt; % TR(s)/nslices = 2/38
        DCM.U.name = [SPM.Sess.U(1:3).name];
        DCM.U.u    = [SPM.Sess.U(1).u(33:end,1) ... % Onset time (encoded in an unknown unit)
                      SPM.Sess.U(2).u(33:end,1) ...
                      SPM.Sess.U(3).u(33:end,1)];
        
        % DCM parameters and options
        %--------------------------------------------------------------------------
        DCM.delays =repmat(SPM.xY.RT/2,DCM.n,1);
        DCM.TE=0.022;
        DCM.options.nonlinear  = 0;
        DCM.options.two_state  = 0;
        DCM.options.stochastic=0;
        DCM.options.centre  = 0;
        DCM.options.induced  = 0;
        
        % Connectivity matrices for model
        %--------------------------------------------------------------------------
        DCM.a = ones(8,8); %Fully specified intrinsic connecitons (between 8 VOIs)
        DCM.b = ones(8,8,3); %Fully specified modulating inputs
%         DCM.c = ones(8,3); %Fully specified driving inputs
        DCM.c = [0 0 0; %L_IPS
                 0 0 0; %L_M1
                 0 0 0; %L_SMA
                 1 1 1; %L_V1
                 0 0 0; %R_IPS
                 0 0 0; %R_SMA
                 1 1 1; %R_V1
                 0 0 0;] %R_DLPFC
        DCM.d = zeros(8,8,0);
        save(fullfile(output_dir,strcat('DCM_model_specified_run',num2str(r),'.mat')),...
            'DCM');
    end
end

%% DCM estimation
part2_dir='D:\Yun-Shiuan_LAMBDA\DCM\Part2_DCM_Specify_and_Estimate_with_reslice_normalized\full_duration\resliced_peaks\Houde_resliced_peak_4mm_C_matrix_only_V1';

for id=31:40%size(subject_list,1)
    % Read in the inclusion run indexes and use them to index runs in the following for loop
    run_rows=strcmp(run_inclusion_index.sub_id,subject_list{id});
    run_valid=table2cell(run_inclusion_index(run_rows,'run_num'));
    run_valid=regexp(run_valid,'(?<=r)\d+','match');
    run_valid=str2double([run_valid{:}]');
    output_dir=fullfile(part2_dir,subject_list{id});
    cd(output_dir);
    % Create DCM.mat per subject run    
    for run=1:length(run_valid)
        r=run_valid(run);% The index refer to the real run number of the task (so it might not be consecutive due to exclusive runs)
        clear matlabbatch
        % Load in the GLM.matper subject run
        DCM_mat=fullfile(part2_dir,...
            subject_list{id},strcat('DCM_model_specified_run',num2str(r),'.mat'));
           
        % The name of the DCM result file
        estimated_mat=fullfile(part2_dir,...
         subject_list{id},strcat('DCM_model_result_estimated_run',num2str(r),'.mat'));
        % Prevent runing if DCM analysis on the run was already done
        if exist(estimated_mat, 'file')==2
            strcat('Already done:',subject_list{id},';run',num2str(r))
        else
        % Run the DCM (Should first convert to matlabbatch)
        matlabbatch{1}.spm.dcm.estimate.dcms.subj.dcmmat = {DCM_mat};
        matlabbatch{1}.spm.dcm.estimate.output.separate = struct([]);
        matlabbatch{1}.spm.dcm.estimate.est_type = 1;
        matlabbatch{1}.spm.dcm.estimate.fmri.analysis = 'time';
        spm_jobman('run',matlabbatch);
        
        % Rename the result DCM and move to the right place
        % (otherwise the estimated results will override the previous DCM.mat)
        old_name=ls('DCM-*Feb*');
        new_name=estimated_mat;
        movefile(old_name,new_name);
        
        strcat('Finish:',subject_list{id},';run',num2str(r))
        end
    end
end

% abs((Ep_raw.A-Ep_normalized.A)./Ep_normalized.A)
% (Ep_raw.A(1,1)-Ep_normalized.A(1,1))./Ep_normalized.A(1,1)
