%% Part0-a: Create new GLMs(for distance effect DCM) - first level design matrix (EPI without reslice)
%% Outline
% Part0-a-i: 6 runs in a single GLM.mat- for VOIs extraction (allow GLM.mat containing multiple or single run)
% --------------------------------------------------------------------------
% Include multiple runs for Part0-a-i is because one could easily use the
% outcomes for second level GLM analysis (categorical appraoch.)
% 1) for "adjusting data" (mean-center the time-series extracted)
% 2) for specifying where the raw data during VOI extraction,

% In order to adjust the data, the GLMs should include the F contrast
% covering all interested conditions (i.e, the three categorical distance
% conditions.)

% Part0-a-ii: 1 single run in a single GLm.mat- for DCM specification (only allow GLM.mat containing single run)
% --------------------------------------------------------------------------
% Include onle one single run is because the DCM program only read in the
% first occured run upon specification.
% 1) GLM for defining the timing of driving/modulating inputs



% Note1: No need to include any covariates. But still need to estimate the GLM to
% derive the F contrast.
% These GLMs are only used for specifying where the raw data are,
% defining the driving/modulating inputs, and the F contrast.

% Note2: Note2: Including covariates has no effect on DCM results if
% they aren’t defined as driving/modulating inputs, since they won’t be included in the DCM at all.

% Note3: However, I still include covariates in the GLMs, since the GLMs are
% going to be estimated anyways. I could take advantage of them(the Part0-a-i GLM.mat)
% and regard these GLMs results as categorical version results.

%% Read in run inclusion index info
addpath(char("D:\GoogleDrive\Lambda_code\m_file\LAMBDA\tool_code"));% enable read_mixed_csv()
run_inclusion_index=cellfun(@(x) regexprep(x,'"',''),...
    read_mixed_csv("D:\Yun-Shiuan_LAMBDA\Run_inclusion_info\inclusive_runs_indexes.csv",','),'un',0);
run_inclusion_index=table(run_inclusion_index(2:end,2),run_inclusion_index(2:end,3),...
    'VariableNames',{'sub_id','run_num'});

%% Read in and Insert task-related information(e.g., onset, response, condition)
onset_var={'fix_onset', 'fraccomp_sti_onset','fraccomp_resp_onset'};
time_var={'fix_onset','fix_duration','fraccomp_sti_onset','fraccomp_resp_onset','fraccomp_resp_RT'};
volumes=strtrim(cellstr(num2str([6:110]', ',%d')));%Create the trimed indexex(6~110) for indexing EIP .nii file

%% Part0-a-i: This GLM is for VOI extraction (which each GLM includes 6 runs)
% --------------------------------------------------------------------------
%spm_jobman('initcfg') -- before running the batch
path='D:\Yun-Shiuan_LAMBDA\First_level_parametric_design_matrix_and_result';
path_first='D:\Yun-Shiuan_LAMBDA\DCM\Part0_First_GLM_for_DCM_multiple_run';
cd(path);
subject_list=cellstr(ls('*')); % subject list
subject_list=subject_list(~cellfun(@isempty,regexp(subject_list,'(?<=df)\d+','match')));% include only subject name
subject_list=subject_list(~ismember(subject_list,{'df1001','df1054'})); %Exclude df1001(no T1) and df1054(incomplete)
%(note: df1014 & de1018 has abnormal amount of runs)
%% Create First Level Batch for all participants
error_first_matrix={};
for i=1:size(subject_list,1)
    clear matlabbatch
    load('D:\Yun-Shiuan_LAMBDA\DCM\First_GLM_for_DCM_single_run\template.mat');
    % Specify outputs of SPM.mat and all first level results
    output_dir=strcat(path_first,'\',subject_list{i});
    matlabbatch{1}.spm.stats.fmri_spec.dir = {output_dir};
    
    % Read the inclusion run indexs and use them to index runs in the following for loop
    run_rows=strcmp(run_inclusion_index.sub_id,subject_list{i});
    run_valid=table2cell(run_inclusion_index(run_rows,'run_num'));
    run_valid=regexp(run_valid,'(?<=r)\d+','match');
    run_valid=str2double([run_valid{:}]');
    if(isnan(run_valid)==1)
        error_first_matrix{i,1}=strcat('Sub',subject_list{i},' no valid runs');
    else
        mkdir(output_dir);
        %Loop through each valid run
        batch_index=0;%Because SPM require 'sess' sumber to be consecutive(i.e.,1,2,3,...), this index will differ from 'r', which refer to the true run number
        
        for run=1:length(run_valid)
            r=run_valid(run);% The index refer to the real run number of the task (so it might not be consecutive due to exclusive runs)
            batch_index=batch_index+1;% see comment above
            
            % Insert the correct images for this run
            scans=strcat('D:\Yun-Shiuan_LAMBDA\preprocessed_data_without_reslice\',...
                subject_list{i},'\',...
                'swarun_',num2str(r),'.nii',volumes);
            matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).scans = scans;
            
            % Read in run info (Onset and Values)
            csv=strcat("D:\Yun-Shiuan_LAMBDA\EprimeData_raw\DevFracs",...
                regexprep(subject_list{i},'df',''),"\",subject_list{i},"_Run",num2str(r),"_tidy.csv");
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
            run_info.fraccomp_sti_dis_abs=abs(run_info.fraccomp_sti_dis);% Create a column for absolute distance
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
            % Get all design matix variables: e.g., onset times
            Near_onset=run_info{strcmp(run_info.fraccomp_sti_dis_type,'Near') & row_non_dummy,'fraccomp_sti_onset'};
            Medium_onset=run_info{strcmp(run_info.fraccomp_sti_dis_type,'Medium') & row_non_dummy,'fraccomp_sti_onset'};
            Far_onset=run_info{strcmp(run_info.fraccomp_sti_dis_type,'Far') & row_non_dummy,'fraccomp_sti_onset'};
            
            % Deal with the missing dummy onset time in df1002 run 1.
            % Use interpolation to fill it up.
            if (strcmp(subject_list{i},'df1002') && run==1)
                dummy_onset=199.678;
            else
                dummy_onset=run_info{row_dummy,'fraccomp_sti_onset'};
            end
            
            current_miss=strcmp(run_info.fraccomp_resp,'');% all miss response(no matter dummy or not)
            resp_onset=run_info{row_non_dummy & ~current_miss,'fraccomp_resp_onset'}; % Only include those which is not null or dummy
            resp_error=double(run_info{row_non_dummy & ~current_miss,'fraccomp_resp_acc'}==0);%Null should be excluded otherewise in the first placeit will be misclassified as 'error'
            
            %Consider to remove the fixation event so that it could be regarded
            %as the
            %baseline-----------------------------------------------------------------------------------------------------
            fix_onset=run_info{:,'fix_onset'};%Including the fixation in the dummy trials
            fix_previous_miss=[0;current_miss(1:end-1)];%Miss of the previous trial affect the following fixation activity
            
            % Insertion into the
            % 'matlabbatch'--------------------------------------------
            matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(1).name = 'Near';
            matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(1).onset = Near_onset;
            
            matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(2).name = 'Medium';
            matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(2).onset = Medium_onset;
            
            matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(3).name = 'Far';
            matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(3).onset = Far_onset;
            
            
            matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(4).name = 'Dummy';
            matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(4).onset = dummy_onset;
            
            matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(5).name = 'Resp';
            matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(5).onset = resp_onset;
            % Remove the 'incorrect response' parametric modulator if there were no such a thing in this run
            if (sum(resp_error)>0)
                matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(5).pmod.name='Error';
                matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(5).pmod.param=resp_error;
            else
                matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(5).pmod=struct([]) ;
            end
            
            matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(6).name = 'Fixation';
            matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(6).onset = fix_onset;
            % Remove the 'Previous Miass response' parametric modulator if there were no such a thing in this run
            if (sum(fix_previous_miss)>0)
                matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(6).pmod.name='Previous_Miss';
                matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(6).pmod.param=fix_previous_miss;
            else
                matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(6).pmod=struct([]);
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
    end
end
%% Run All First Level Batches (Part0-a-i)
% To find subjects that have valid runs(n=40)
cd('D:\Yun-Shiuan_LAMBDA\First_level_parametric_design_matrix_and_result');
subject_list=cellstr(ls('*')); % subject list (which has valid first level batch, i.e., not being excluded)
subject_list=subject_list(~cellfun(@isempty,regexp(subject_list,'(?<=df)\d+','match')));% include only subject name

error_first_level_trail=cell(length(subject_list),1);
for i=1:length(subject_list)
    try
        load(strcat('D:\Yun-Shiuan_LAMBDA\DCM\First_GLM_for_DCM\',...
            subject_list{i},'\first_level_batch.mat'));
        spm_jobman('run',matlabbatch);
        strcat('Finish:',subject_list{i})
    catch ME
        error_first_level_trail{i,1}=strcat('Sub',subject_list{i},': ',ME.message)
    end
end
error_first_level_trail2=error_first_level_trail;
save('error_first_level_record_trail2.mat','error_first_level_trail')

%% Part0-a-ii: This GLM is for DCM specification (which each GLM includes 1 single run)
% --------------------------------------------------------------------------
path='D:\Yun-Shiuan_LAMBDA\First_level_parametric_design_matrix_and_result';
path_first='D:\Yun-Shiuan_LAMBDA\DCM\Part0_First_GLM_for_DCM_single_run';
cd(path);
subject_list=cellstr(ls('*')); % subject list
subject_list=subject_list(~cellfun(@isempty,regexp(subject_list,'(?<=df)\d+','match')));% include only subject name
subject_list=subject_list(~ismember(subject_list,{'df1001','df1054'})); %Exclude df1001(no T1) and df1054(incomplete)
%(note: df1014 & de1018 has abnormal amount of runs)

error_first_matrix={};
for i=1:size(subject_list,1)
    % Read the inclusion run indexs and use them to index runs in the following for loop
    run_rows=strcmp(run_inclusion_index.sub_id,subject_list{i});
    run_valid=table2cell(run_inclusion_index(run_rows,'run_num'));
    run_valid=regexp(run_valid,'(?<=r)\d+','match');
    run_valid=str2double([run_valid{:}]');
    
    %Loop through each valid run
    for run=1:length(run_valid)
        r=run_valid(run);% The index refer to the real run number of the task (so it might not be consecutive due to exclusive runs)
        clear matlabbatch
        load('D:\Yun-Shiuan_LAMBDA\DCM\Part0_First_GLM_for_DCM_single_run\template.mat');
        % Specify outputs of SPM.mat and all first level results
        output_dir=fullfile(path_first,subject_list{i},strcat('run',num2str(r)));
        matlabbatch{1}.spm.stats.fmri_spec.dir = {output_dir};
        
        if(isnan(run_valid)==1)
            error_first_matrix{i,1}=strcat('Sub',subject_list{i},' no valid runs');
        else
            mkdir(output_dir);
            % Insert the correct images for this run
            scans=strcat('D:\Yun-Shiuan_LAMBDA\preprocessed_data_without_reslice\',...
                subject_list{i},'\',...
                'swarun_',num2str(r),'.nii',volumes);
            matlabbatch{1}.spm.stats.fmri_spec.sess.scans = scans;
            
            % Read in run info (Onset and Values)
            csv=strcat("D:\Yun-Shiuan_LAMBDA\EprimeData_raw\DevFracs",...
                regexprep(subject_list{i},'df',''),"\",subject_list{i},"_Run",num2str(r),"_tidy.csv");
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
            run_info.fraccomp_sti_dis_abs=abs(run_info.fraccomp_sti_dis);% Create a column for absolute distance
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
            % Get all design matix variables: e.g., onset times
            Near_onset=run_info{strcmp(run_info.fraccomp_sti_dis_type,'Near') & row_non_dummy,'fraccomp_sti_onset'};
            Medium_onset=run_info{strcmp(run_info.fraccomp_sti_dis_type,'Medium') & row_non_dummy,'fraccomp_sti_onset'};
            Far_onset=run_info{strcmp(run_info.fraccomp_sti_dis_type,'Far') & row_non_dummy,'fraccomp_sti_onset'};
            
            % Deal with the missing dummy onset time in df1002 run 1.
            % Use interpolation to fill it up.
            if (strcmp(subject_list{i},'df1002') && run==1)
                dummy_onset=199.678;
            else
                dummy_onset=run_info{row_dummy,'fraccomp_sti_onset'};
            end
            
            current_miss=strcmp(run_info.fraccomp_resp,'');% all miss response(no matter dummy or not)
            resp_onset=run_info{row_non_dummy & ~current_miss,'fraccomp_resp_onset'}; % Only include those which is not null or dummy
            resp_error=double(run_info{row_non_dummy & ~current_miss,'fraccomp_resp_acc'}==0);%Null should be excluded otherewise in the first placeit will be misclassified as 'error'
            
            %Consider to remove the fixation event so that it could be regarded
            %as the
            %baseline-----------------------------------------------------------------------------------------------------
            fix_onset=run_info{:,'fix_onset'};%Including the fixation in the dummy trials
            fix_previous_miss=[0;current_miss(1:end-1)];%Miss of the previous trial affect the following fixation activity
            
            % Insertion into the
            % 'matlabbatch'--------------------------------------------
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).name = 'Near';
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).onset = Near_onset;
            
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).name = 'Medium';
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).onset = Medium_onset;
            
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).name = 'Far';
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).onset = Far_onset;
            
            
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).name = 'Dummy';
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).onset = dummy_onset;
            
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).name = 'Resp';
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).onset = resp_onset;
            % Remove the 'incorrect response' parametric modulator if there were no such a thing in this run
            if (sum(resp_error)>0)
                matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).pmod.name='Error';
                matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).pmod.param=resp_error;
            else
                matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).pmod=struct([]) ;
            end
            
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).name = 'Fixation';
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).onset = fix_onset;
            % Remove the 'Previous Miass response' parametric modulator if there were no such a thing in this run
            if (sum(fix_previous_miss)>0)
                matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).pmod.name='Previous_Miss';
                matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).pmod.param=fix_previous_miss;
            else
                matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).pmod=struct([]);
            end
            motion_parameters_file=strcat('D:\Yun-Shiuan_LAMBDA\preprocessed_data_without_reslice\',...
                subject_list{i},'\rp_arun_',...
                num2str(r),'.txt');
            matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {motion_parameters_file};
            
        end
        
        % Save the matlabbatch
        save(strcat(output_dir,'\first_level_batch'),'matlabbatch');
    end
end
%% Run All First Level Batches (Part0-a-ii)
%spm_jobman('initcfg') -- before running the batch
% To find subjects that have valid runs(n=40)
cd('D:\Yun-Shiuan_LAMBDA\First_level_parametric_design_matrix_and_result');
subject_list=cellstr(ls('*')); % subject list (which has valid first level batch, i.e., not being excluded)
subject_list=subject_list(~cellfun(@isempty,regexp(subject_list,'(?<=df)\d+','match')));% include only subject name

error_first_level_trail=cell(length(subject_list),1);
for i=1:length(subject_list)
    % Read the inclusion run indexs and use them to index runs in the following for loop
    run_rows=strcmp(run_inclusion_index.sub_id,subject_list{i});
    run_valid=table2cell(run_inclusion_index(run_rows,'run_num'));
    run_valid=regexp(run_valid,'(?<=r)\d+','match');
    run_valid=str2double([run_valid{:}]');
    for run=1:length(run_valid)
        r=run_valid(run);% The index refer to the real run number of the task (so it might not be consecutive due to exclusive runs)
        try
            load(strcat('D:\Yun-Shiuan_LAMBDA\DCM\Part0_First_GLM_for_DCM_single_run\',...
                subject_list{i},'\',strcat('run',num2str(r)),'\first_level_batch.mat'));
            spm_jobman('run',matlabbatch);
            strcat('Finish:',subject_list{i},';run',num2str(r))
        catch ME
            error_first_level_trail{i,1}=strcat('Sub',subject_list{i},': ',ME.message)
        end
    end
end
error_first_level_trail2=error_first_level_trail;
save('error_first_level_record_trail2.mat','error_first_level_trail')
