%% Part1-a: DCM - Extract time series in VOIs (for distance effect DCM)
addpath(char("D:\GoogleDrive\Lambda_code\m_file\LAMBDA\tool_code"));% enable read_mixed_csv()
%spm_jobman('initcfg') -- before running the batch

path='D:\Yun-Shiuan_LAMBDA\DCM\Part0_First_GLM_for_DCM'; % To get the subject list with valid runs 
cd(path);
subject_list=cellstr(ls('*')); % subject list (which has valid first level contrast, i.e., not being excluded)
subject_list=subject_list(~cellfun(@isempty,regexp(subject_list,'(?<=df)\d+','match')));% include only subject name

%% Read in run inclusion index info
run_inclusion_index=cellfun(@(x) regexprep(x,'"',''),...
                read_mixed_csv("D:\Yun-Shiuan_LAMBDA\Run_inclusion_info\inclusive_runs_indexes.csv",','),'un',0);
run_inclusion_index=table(run_inclusion_index(2:end,2),run_inclusion_index(2:end,3),...
                    'VariableNames',{'sub_id','run_num'});
%% Batch for VOI Extraction
% 8 regions invloving in distance sensitivity
VOI_name={"R_V1","L_V1",...
    "R_IPS","L_IPS",...
    "R_SMA","L_SMA",...
    "L_M1","R_VLPFC"}';
peak_coordinates={...
    [21 -100 -4], [-24,-97,-7],...
    [45 -40 47],[-42 -40 41],...
    [6 17 50],[-6 11 50],...
    [-42 -1 32],[48 35 23]};
common_mask='D:\Yun-Shiuan_LAMBDA\template\for_this_study\binary_overlapped_WM_GM_common_mask_n40_with_signal.nii,1';

for i=1:size(subject_list,1)
    % Read in the inclusion run indexes and use them to index runs in the following for loop
    run_rows=strcmp(run_inclusion_index.sub_id,subject_list{i});        
    run_valid=table2cell(run_inclusion_index(run_rows,'run_num'));
    run_valid=regexp(run_valid,'(?<=r)\d+','match');
    run_valid=str2double([run_valid{:}]');
    
    % Make directory (for VOI output)
    outputdir=strcat('D:\Yun-Shiuan_LAMBDA\DCM\Part1_VOIs_for_DCM\',subject_list{i});
    mkdir(outputdir);
    outputdir_batch=strcat(outputdir,'\VOI_batch');
    mkdir(outputdir_batch);
    outputdir_data=strcat(outputdir,'\VOI_data');
    mkdir(outputdir_data);

    for run=1:length(run_valid)
        r=run_valid(run);% The index refer to the real run number of the task (so it might not be consecutive due to exclusive runs)
        for v=1:size(VOI_name,1)
            clear matlabbatch
            % Make batch for VOI time series extraction (per subject run)
            matlabbatch{1}.spm.util.voi.spmmat = {strcat('D:\Yun-Shiuan_LAMBDA\DCM\Part0_First_GLM_for_DCM\',subject_list{i},'\SPM.mat')};
            matlabbatch{1}.spm.util.voi.adjust = 1;
            matlabbatch{1}.spm.util.voi.session = 1;
            matlabbatch{1}.spm.util.voi.name = char(strcat(VOI_name{v},'_run',num2str(r)));
            matlabbatch{1}.spm.util.voi.roi{1}.sphere.centre = peak_coordinates{v};
            matlabbatch{1}.spm.util.voi.roi{1}.sphere.radius = 8;
            matlabbatch{1}.spm.util.voi.roi{1}.sphere.move.fixed = 1;
            matlabbatch{1}.spm.util.voi.roi{2}.mask.image = {common_mask};
            matlabbatch{1}.spm.util.voi.roi{2}.mask.threshold = 0.9;
            matlabbatch{1}.spm.util.voi.expression = 'i1 & i2';
            save(strcat(outputdir_batch,'\',VOI_name{v},'_run',num2str(r)),'matlabbatch') 
            
            %Run the batch
            spm_jobman('run',matlabbatch);
            sourcefiles=strcat('D:\Yun-Shiuan_LAMBDA\DCM\Part0_First_GLM_for_DCM\',subject_list{i},'\');
            movefile(strcat(sourcefiles,'\VOI*'),outputdir_data);
        end
    end
end
%Check if the eigen time series is similar to the averaged time series
%They are actually quite different!!
% mean_time_series=spm_summarise('D:\Yun-Shiuan_LAMBDA\preprocessed_data_without_reslice\df1002\swarun_1.nii',...
%     struct('def','sphere', 'spec',8, 'xyz',[45 -40 47]'),...
%     @mean);
% mean_time_series=mean_time_series(6:end);
% mean_time_series=(mean_time_series-mean(mean_time_series))/std(mean_time_series);