%% Part 0-a-i: RSA : Estimate Beta for each trial
% RDA_ID_continuous:
% Utalize GLM which sets each trial as a unique predictor (36x6=218 unique trials in total per person)
% By regarding distance as a continuous variable, there's no repetition of trials.
% Disadvantages: 
% 1) The toolbox does not allow conceptual model to differ accross session
% 2) The statistical efficiency decrease without trial repetition
% Solution: 
% Regard distance as descrete variable(N/M/F) so that 
% 1)Each run is composed of a identical set of trials (4 formats x 2 positions x 3 distance =  24 unique ID per each subject run)
% See "RSA_GLM_estimate_t_value_RSA_ID_discrete.m"

%spm_jobman('initcfg') -- before running the batch
addpath('D:\GoogleDrive\Lambda_code\m_file\LAMBDA\tool_code');% enable read_mixed_csv()

%% Some constants
%files
file_valid_run='D:\Yun-Shiuan_LAMBDA\Run_inclusion_info\inclusive_runs_indexes.csv';
GLM_template='D:\Yun-Shiuan_LAMBDA\RSA\Part0-a_GLM_estimate_beta_value\template.mat';
%paths
RSA_GLM_path='D:\Yun-Shiuan_LAMBDA\RSA\Part0-a_GLM_estimate';
preprocess_data_path='D:\Yun-Shiuan_LAMBDA\preprocessed_data_with_reslice_trial2';
beh_data_path='D:\Yun-Shiuan_LAMBDA\EprimeData_raw\DevFracs';
%Other constants
onset_var={'fix_onset', 'fraccomp_sti_onset','fraccomp_resp_onset'};
time_var={'fix_onset','fix_duration','fraccomp_sti_onset','fraccomp_resp_onset','fraccomp_resp_RT','fraccomp_sti_duration'};
volumes=strtrim(cellstr(num2str([6:110]', ',%d')));%Create the trimed indexex(6~110) for indexing EIP .nii file
cov_list={'Dummy','Resp','Fixation'};
%% Read in run inclusion index info
run_inclusion_index=cellfun(@(x) regexprep(x,'"',''),...
    read_mixed_csv(file_valid_run,','),'un',0);
run_inclusion_index=table(run_inclusion_index(2:end,2),run_inclusion_index(2:end,3),...
    'VariableNames',{'sub_id','run_num'});
%% Derive subjects with valid runs
subject_list=unique(run_inclusion_index.sub_id);

%% Specify the GLM to extract the beta series
error_first_matrix=cell(size(subject_list,1),1);
for i=1:size(subject_list,1)
    clear matlabbatch
    load(GLM_template)% Load in the matlabbatch
    
    % Specify the output path
    output_dir=fullfile(RSA_GLM_path,subject_list{i});
    matlabbatch{1}.spm.stats.fmri_spec.dir = {output_dir};
    
    % Read the inclusion run indexs and use them to index runs in the following while loop
    run_rows=strcmp(run_inclusion_index.sub_id,subject_list{i});
    run_valid=table2cell(run_inclusion_index(run_rows,'run_num'));
    run_valid=regexp(run_valid,'(?<=r)\d+','match');
    run_valid=str2double([run_valid{:}]');
    if(isnan(run_valid)==1)
        error_first_matrix{i,1}=strcat('Sub',subject_list{i},' no valid runs');
    else
        % SPM require 'sess' sumber to be consecutive(i.e.,1,2,3,...), this index will differ from 'r', which refer to the true run number
        mkdir(output_dir);
        
        %Loop through each valid run
        batch_index=0;%Because SPM require 'sess' sumber to be consecutive(i.e.,1,2,3,...), this index will differ from 'r', which refer to the true run number
        for run=1:length(run_valid)
            r=run_valid(run);% The index refer to the real run number of the task (so it might not be consecutive due to excluded runs)
            batch_index=batch_index+1;% see comment above
            
            % Insert the correct images for this run 
            %(EPI (normalized but UNSMOOTHED images -- retain the detailed pattern))
            scans=fullfile(preprocess_data_path,subject_list{i},...
                strcat('warrun_',num2str(r),'.nii',volumes));
            matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).scans = scans;
            
            % Read in run info (Onset and Values)
            csv=fullfile(strcat(beh_data_path,regexprep(subject_list{i},'df','')),...
                char(strcat(subject_list{i},"_Run",num2str(r),"_tidy.csv")));
            run_info=cellfun(@(x) regexprep(x,'"',''),...
                read_mixed_csv(csv,','),'un',0);
            run_info=cell2table(run_info(2:end,2:end), 'VariableNames',run_info(1,2:end));
            % Mutate the following variable type to numeric
            run_info.fix_onset=str2double(run_info.fix_onset);
            run_info.fix_duration=str2double(run_info.fix_duration);
            run_info.fraccomp_sti_onset=str2double(run_info.fraccomp_sti_onset);
            run_info.fraccomp_sti_dis=str2double(run_info.fraccomp_sti_dis);
            run_info.fraccomp_resp_onset=str2double(run_info.fraccomp_resp_onset);
            run_info.fraccomp_resp_RT=str2double(run_info.fraccomp_resp_RT);
            run_info.fraccomp_resp_acc=str2double(run_info.fraccomp_resp_acc);
            run_info.fraccomp_sti_value_left=str2double(run_info.fraccomp_sti_value_left);
            run_info.fraccomp_sti_value_right=str2double(run_info.fraccomp_sti_value_right);
            run_info.fraccomp_sti_duration=str2double(run_info.fraccomp_sti_duration);
            
            % Create a column for absolute distance
            run_info.fraccomp_sti_dis_abs=abs(run_info.fraccomp_sti_dis);
            % Adjust the onset time: set the time=0 point to t=10s point(first fixation appears), instead of the starting time of Eprime
            t0=min(run_info.fix_onset);
            for v=1:length(onset_var)
                run_info{:,onset_var{v}}=run_info{:,onset_var{v}}-t0;
            end
            
            for v=1:length(time_var) %Covert from ms to second
                run_info{:,time_var{v}}=run_info{:,time_var{v}}/1000;
            end
            % Get the row indexes of non-dummy and hummy trials
            row_dummy=strcmp(run_info.trial_order,'NA');
            row_non_dummy=~row_dummy;
            % Get all design matix variables: e.g., onset times and distance value
            trial_id_RSA=run_info{row_non_dummy,'trial_id_RSA'};
            trial_onset=run_info{row_non_dummy,'fraccomp_sti_onset'};
            trial_duration=run_info{row_non_dummy,'fraccomp_sti_duration'};
            % Deal with the missing onset time in df1002 run 1.
            % Use interpolation to fill it up.
            if (strcmp(subject_list{i},'df1002') && run==1)
                dummy_onset=199.678;
            else
                dummy_onset=run_info{row_dummy,'fraccomp_sti_onset'};
            end
            
            current_miss=strcmp(run_info.fraccomp_resp,'');% all miss response(no matter dummy or not)
            resp_onset=run_info{row_non_dummy & ~current_miss,'fraccomp_resp_onset'}; % Only include those which is not null or dummy
            resp_error=double(run_info{row_non_dummy & ~current_miss,'fraccomp_resp_acc'}==0);%Null should be excluded otherewise in the first placeit will be misclassified as 'error'
            
            fix_onset=run_info{:,'fix_onset'};%Including the fixation in the dummy trials
            fix_previous_miss=[0;current_miss(1:end-1)];%Miss of the previous trial affect the following fixation activity
            
            % Insertion into the
            % 'matlabbatch'--------------------------------------------
            %% RSA trial ID (Note to set its duration to 4s to maximize the activition captured)
            % Note:
            % Since SPM by default names contrast/parameter from 1 consecutively,
            % one coulds leave the output file numbers the way they are. By this way,
            % in order to take the RSA trial ID into account,
            % one should deal with it while creating the RDMs.
            % (But still name the condition by RSA trail ID so that it might be helpful in the future analysis)
            
            
            for trial_id_RSA_index=1:length(trial_id_RSA)
                % Trail-by-Trail preictors
                matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(trial_id_RSA_index).name = trial_id_RSA{trial_id_RSA_index};
                matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(trial_id_RSA_index).onset = trial_onset(trial_id_RSA_index);
                matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(trial_id_RSA_index).duration = trial_duration(trial_id_RSA_index);
                
                matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(trial_id_RSA_index).tmod = 0;
                matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(trial_id_RSA_index).pmod = struct('name', {}, 'param', {}, 'poly', {});
                matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(trial_id_RSA_index).orth = 0;
            end
            
            % Task-related Covariates (e.g., dummy, responses, fixation)
            for cov=1:size(cov_list,2)
                matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(trial_id_RSA_index+cov).duration=0;
                matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(trial_id_RSA_index+cov).orth = 0;
                matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(trial_id_RSA_index+cov).tmod = 0;
                matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(trial_id_RSA_index+cov).pmod= struct('name', {}, 'param', {}, 'poly', {}) ;
            end
            
            % Dummy --------------------
            matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(trial_id_RSA_index+1).name = 'Dummy';
            matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(trial_id_RSA_index+1).onset = dummy_onset;
            
            % Resp --------------------
            matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(trial_id_RSA_index+2).name = 'Resp';
            matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(trial_id_RSA_index+2).onset = resp_onset;
            
            % Remove the 'incorrect response' parametric modulator if there were no such a thing in this run
            if (sum(resp_error)>0)
                matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(trial_id_RSA_index+2).pmod(1).name='Error';
                matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(trial_id_RSA_index+2).pmod.param=resp_error;
                matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(trial_id_RSA_index+2).pmod.poly = 1;
            else
                matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(trial_id_RSA_index+2).pmod=struct([]) ;
            end
            
            % Fixation -----------------
            matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(trial_id_RSA_index+3).name = 'Fixation';
            matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(trial_id_RSA_index+3).onset = fix_onset;
            % Remove the 'Previous Miass response' parametric modulator if there were no such a thing in this run
            if (sum(fix_previous_miss)>0)
                matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(trial_id_RSA_index+3).pmod(1).name='Previous_Miss';
                matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(trial_id_RSA_index+3).pmod.param=fix_previous_miss;
                matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(trial_id_RSA_index+3).pmod.poly = 1;
            else
                matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(trial_id_RSA_index+3).pmod=struct([]);
            end
            motion_parameters_file=strcat('D:\Yun-Shiuan_LAMBDA\preprocessed_data_without_reslice\',...
                subject_list{i},'\rp_arun_',...
                num2str(r),'.txt');
            matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).multi_reg = {motion_parameters_file};
        end
        % Remove redundant sessions if less than 6 valid runs
        remove_index=length(run_valid)+1; %remove the index which equals to amount+1
        while(batch_index < 6)
            batch_index=batch_index+1;
            matlabbatch{1}.spm.stats.fmri_spec.sess(remove_index)=[];
        end
        % Save the matlabbatch
        save(strcat(output_dir,'\first_level_batch'),'matlabbatch');
%         % run the batch
%         spm_jobman('run',matlabbatch);
    end
end

%% Run the GLM to extract the beta series
error_first_matrix=cell(size(subject_list,1),1);
for i=1:size(subject_list,1)
    try
        clear matlabbatch
         % Specify the SPM file
        output_dir=fullfile(RSA_GLM_path,subject_list{i});
        SPM_file=fullfile(output_dir,'first_level_batch');

        % Load in the matlabbatch
        load(SPM_file)

        % run the batch  
        spm_jobman('run',matlabbatch);
    catch ME
        error_first_matrix{i,1}=strcat('id:',num2str(i),':',ME.message);
    end
end
