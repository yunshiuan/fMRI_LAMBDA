%% Part 0-a-iv: RSA : Estimate Beta for each trial - t value
% Part0-a: Specifying GLMs to extract t avlues of each trial
% iii: use t maps (36 trials with discrete distance) 
% instead of beta maps (144  unique trials with continuous distance, or 24 unique trials.)

% iv: RDA_ID_discrete 18 ID:
% - RSA ID : 1~18 (3 formats x 3 distance levels x 2 directions)
% (2) The motivation is to increase the statistical efficiency (Conditions in the Trial 4 is too noisy),
% while still keep that same-notation trials have the same statistical efficiency
% as the cross-notation trials (as Trial 4 did).
% [LL: 12 trials = 3 distance levels x 2 direction x 2 repetition]
% [FF: 12 trials = 3 distance levels x 2 direction x 2 repetition]
% [cross-notation: 12 trials = 3 distance levels x 2 direction x 2 [repetition]
% (Having more repetition in a condition entails a smaller degree of random noise.)

% The motivation for having 36 trials is to ensure beta estimate of each condition is subjected to same amount of noise.
% (Having more repetition in a condition entails a smaller degree of random noise.)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Update the subject list for children: use
%"inclusive_runs_indexes_new_Jun_10.csv" (2018/6/10)
spm_jobman('initcfg');


PATH_ROOT='D:\Yun-Shiuan_LAMBDA';
PATH_TOOL_CODE='D:\GoogleDrive\Lambda_code\m_file\LAMBDA\tool_code';
addpath(genpath(fullfile(PATH_ROOT,'rsatoolbox')));
addpath(PATH_TOOL_CODE);% enable read_mixed_csv()

%% Some constants
%paths
PATH_RSA_GLM=fullfile(PATH_ROOT,'RSA','trial_21','Part0-a_GLM_estimate_t_value','child');mkdir(PATH_RSA_GLM);
PATH_PREPROCESS_DATA=fullfile(PATH_ROOT,'preprocessed_data_with_reslice_trial2');
PATH_BEH_DATA=fullfile(PATH_ROOT,'EprimeData_raw');
%files
FILE_VALID_RUN=fullfile(PATH_ROOT,'Run_inclusion_info','inclusive_runs_indexes_new_June_10.csv');
FILE_GLM_TEMPLATE=fullfile(PATH_RSA_GLM,'template.mat');
%Other constants
ONSET_VAR={'fix_onset', 'fraccomp_sti_onset','fraccomp_resp_onset'};
TIME_VAR={'fix_onset','fix_duration','fraccomp_sti_onset','fraccomp_resp_onset','fraccomp_resp_RT','fraccomp_sti_duration'};
VOLUMES=strtrim(cellstr(num2str([6:110]', ',%d')));%Create the trimed indexex(6~110) for indexing EIP .nii file
COV_LIST={'Dummy','Resp','Fixation'};
RSA_ID_TYPE='trial_id_RSA_discrete_18'; %use the discrete version of RSA ID
AMOUNT_COV=numel(COV_LIST); % Amount of covariates
AMOUNT_RSA_ID=18; % amount of unique RSA id

%% Read in run inclusion index info
run_inclusion_index=read_mixed_csv_to_table(FILE_VALID_RUN);
%% Derive subjects with valid runs
subject_list=unique(run_inclusion_index.sub_id);
% subject_list=subject_list(~ismember(subject_list,{'XFC305'})); %305 have
% missing run (already solved)
%% Specify the GLM to extract the beta series
error_first_matrix=cell(size(subject_list,1),1);
for id=1:numel(subject_list)
    % Specify the output path
    path_output=fullfile(PATH_RSA_GLM,subject_list{id});
    
    % Skip the subject if already done
    if ~isempty(dir(fullfile(path_output,'spmT*')))
        strcat('Already done: ',subject_list{id})
        continue;
    end
    
    % Load in the template matlabbatch
    clear matlabbatch
    load(FILE_GLM_TEMPLATE)
    matlabbatch{1}.spm.stats.fmri_spec.dir = {path_output};
    % Read the inclusion run indexs and use them to index runs in the following while loop
    run_rows=strcmp(run_inclusion_index.sub_id,subject_list{id});
    run_valid=table2cell(run_inclusion_index(run_rows,'run_num'));
    run_valid=regexp(run_valid,'(?<=r)\d+','match');
    run_valid=str2double([run_valid{:}]');
    if(isnan(run_valid)==1)
        error_first_matrix{id,1}=strcat('Sub',subject_list{id},' no valid runs');
    else
        % SPM require 'sess' sumber to be consecutive(i.e.,1,2,3,...), this index will differ from 'r', which refer to the true run number
        mkdir(path_output);
        
        %Loop through each valid run
        for run=1:length(run_valid)
            r=run_valid(run);% The index refer to the real run number of the task (so it might not be consecutive due to excluded runs)
            batch_index=run;% %Because SPM require 'sess' sumber to be consecutive(i.e.,1,2,3,...), this index will differ from 'r', which refer to the true run number
            
            % Insert the correct images for this run
            %(EPI (normalized but UNSMOOTHED images -- retain the detailed pattern))
            scans=fullfile(PATH_PREPROCESS_DATA,subject_list{id},...
                strcat('warrun_',num2str(r),'.nii',VOLUMES));
            matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).scans = scans;
            
            % Read in run info (Onset and Values)
            subject_folder_name=regexprep(subject_list{id},'df','DevFracs');
            
            csv=fullfile(PATH_BEH_DATA,subject_folder_name,...
                char(strcat(subject_list{id},"_Run",num2str(r),"_tidy.csv")));
            run_info=read_mixed_csv_to_table(csv);
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
            for v=1:length(ONSET_VAR)
                run_info{:,ONSET_VAR{v}}=run_info{:,ONSET_VAR{v}}-t0;
            end
            
            for v=1:length(TIME_VAR) %Covert from ms to second
                run_info{:,TIME_VAR{v}}=run_info{:,TIME_VAR{v}}/1000;
            end
            % Get the row indexes of non-dummy and hummy trials
            row_dummy=(strcmp(run_info.trial_order,'NA')|strcmp(run_info.trial_order,'NULL'));
            row_non_dummy=~row_dummy;
            % Get all design matix variables: e.g., onset times and distance value
            trial_id_RSA=run_info{row_non_dummy,RSA_ID_TYPE};
            trial_onset=run_info{row_non_dummy,'fraccomp_sti_onset'};
            trial_duration=run_info{row_non_dummy,'fraccomp_sti_duration'};
            % Deal with the missing onset time in df1002 run 1.
            % Use interpolation to fill it up.
            if (strcmp(subject_list{id},'df1002') && run==1)
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
            
            % Trail-by-Trail preictors (from RSA_id = 1 to 18)
            
            % Check if the unique trials amount is 18. If not, stop running
            if(length(unique(trial_id_RSA))~=AMOUNT_RSA_ID)
                error(char(strcat("Unique trials amount is incorrect! - ","sub id: ", subject_list{id},"; run: ", num2str(run))))
                return
            end
            
            for trial_id_RSA_index=1:length(unique(trial_id_RSA)) % Each unique trial for an iteration (24x1)
                % The index(36x1) of trials that belong to the same RSA ID
                trials_id_RSA_indices=strcmp(trial_id_RSA,num2str(trial_id_RSA_index));
                % Trail-by-Trail preictors
                matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(trial_id_RSA_index).name = num2str(trial_id_RSA_index);
                matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(trial_id_RSA_index).onset = trial_onset(trials_id_RSA_indices);
                matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(trial_id_RSA_index).duration = trial_duration(trials_id_RSA_indices);
                
                matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(trial_id_RSA_index).tmod = 0;
                matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(trial_id_RSA_index).pmod = struct('name', {}, 'param', {}, 'poly', {});
                matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(trial_id_RSA_index).orth = 0;
                
            end
            
            % Task-related Covariates (e.g., dummy, responses, fixation)
            for cov=1:size(COV_LIST,2)
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
            % Remove the 'Previous Miss response' parametric modulator if there were no such a thing in this run
            if (sum(fix_previous_miss)>0)
                matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(trial_id_RSA_index+3).pmod(1).name='Previous_Miss';
                matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(trial_id_RSA_index+3).pmod.param=fix_previous_miss;
                matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(trial_id_RSA_index+3).pmod.poly = 1;
            else
                matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(trial_id_RSA_index+3).pmod=struct([]);
            end
            % Motion Parameters ---------------
            motion_parameters_file=fullfile(PATH_PREPROCESS_DATA,...
                subject_list{id},['rp_run_',num2str(r),'.txt']);
            matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).multi_reg = {motion_parameters_file};
            
            % Remove redundant predictors --------------
            remove_index = trial_id_RSA_index + AMOUNT_COV + 1;
            while (1) % Run until broke by the out of bound error
                try
                    matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(remove_index)=[];
                    % Note that no need to "remove_index=remove_index+1",
                    % since the cell will collapse to fill the hole
                    % automatically
                catch ME
                    break % Break the while loop when exceeding the boundary
                end
            end
            
        end% End of the run loop ---------------------------------------
        % Remove redundant sessions if less than 6 valid runs
        remove_index=length(run_valid)+1; %remove the index which equals to amount+1
        while(batch_index < 6)
            batch_index=batch_index+1;
            matlabbatch{1}.spm.stats.fmri_spec.sess(remove_index)=[];
        end
        % Contrast Manager
        % -----------------------------------------------------------------------------------------
        matlabbatch{3}.spm.stats.con.spmmat(1) = ...
            cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
        % Each unique trial for an iteration (24x1) (from RSA_id = 1 to 24)
        for trial_id_RSA_index=1:length(unique(trial_id_RSA))
            % The index(36x1) of trials that belong to the same RSA ID
            trials_id_RSA_indices=strcmp(trial_id_RSA,num2str(trial_id_RSA_index));
            matlabbatch{3}.spm.stats.con.consess{trial_id_RSA_index}.tcon.name = num2str(trial_id_RSA_index);
            % Note the logic of padding here (the order has already been adjusted to line up with RSA ID)
            matlabbatch{3}.spm.stats.con.consess{trial_id_RSA_index}.tcon.weights = [zeros(1,trial_id_RSA_index-1) 1]; 
            matlabbatch{3}.spm.stats.con.consess{trial_id_RSA_index}.tcon.sessrep = 'replsc';
        end
        matlabbatch{3}.spm.stats.con.delete = 1;
        
        % Save the matlabbatch --------------------------------------------------------------------
        save(fullfile(path_output,'first_level_batch.mat'),'matlabbatch');
        % run the batch
%         spm_jobman('run',matlabbatch);
    end % If-else clause of "not valid run"
end

%% Run the GLM 
error_first_matrix=cell(size(subject_list,1),1);
for id=1:numel(subject_list)
    try
        % Specify the SPM file
        path_output=fullfile(PATH_RSA_GLM,subject_list{id});
        
        %Skip the subject if already done
        if ~isempty(dir(fullfile(path_output,'spmT*')))
            strcat('Already done: ',subject_list{id})
            continue;
        end        
        
        % Load in the matlabbatch
        SPM_file=fullfile(path_output,'first_level_batch');
        clear matlabbatch
        load(SPM_file)
        
        % run the batch
        spm_jobman('run',matlabbatch);
    catch ME
        error_first_matrix{id,1}=strcat('id:',num2str(id),':',ME.message);
    end
end
