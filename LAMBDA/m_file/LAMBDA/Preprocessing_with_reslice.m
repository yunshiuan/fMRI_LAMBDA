%% Part 1-b - Preprocessing with reslice
%spm_jobman('initcfg') -- before running the batch
%Enable the usage of helper functions----------
spm_jobman('initcfg');

PATH_ROOT='D:\Yun-Shiuan_LAMBDA';
PATH_TOOL_CODE='D:\GoogleDrive\Lambda_code\m_file\LAMBDA\tool_code';
addpath(PATH_TOOL_CODE);% enable read_mixed_csv()

% Constants------------------------------------
PATH_RAW_DATA=fullfile(PATH_ROOT,'raw_data');
PATH_PROCESSED_DATA=fullfile(PATH_ROOT,'preprocessed_data_with_reslice_trial2');
PATH_PROCESSED_MAT=fullfile(PATH_ROOT,'preprocessing_mat','with_resliced_trial2');
FILE_VALID_RUN=fullfile(PATH_ROOT,'Run_inclusion_info','inclusive_runs_indexes_new_April_3.csv');
VALUE_EPI_SLICE_START=6;
VALUE_EPI_SLICE_END=110;
STRING_T1_NAME='T1_rescaled.nii'; % For child: 'T1_rescaled.nii'; For Adult: 'T1.nii'
% Read in run inclusion index info 
% and derive subjects with valid runs
run_inclusion_index=read_mixed_csv_to_table(FILE_VALID_RUN);
subject_list=unique(run_inclusion_index.sub_id);

% subject_list=subject_list(~ismember(subject_list,{'XFC305'})); %305 have
% missing run (Already sovled)

%Child: Exclude df1001(no T1) ,df1026(missing T1) ,and df1054(incomplete)
%(note: df1014  has abnormal amount of runs)
%Adult:
%(1)XFC305: Only 5 runs -Not recorded on the log--> Go to check in the external drive.
%(3)XFC313: 7 runs -> Should exclude Run4(Scan9).

%% Preprocessing: Create .mat file for preprocessing
for id=1:length(subject_list)
    %Some Constants
    path_processed_result=fullfile(PATH_PROCESSED_DATA,subject_list{id});
    file_save_mat_result=fullfile(PATH_PROCESSED_MAT,['pre_',subject_list{id},'.mat']);
    volumes=strtrim(cellstr(num2str((VALUE_EPI_SLICE_START:VALUE_EPI_SLICE_END)', ',%d')));
    
    epi_run_1=cellstr(fullfile(PATH_RAW_DATA,subject_list{id},'nii_raw',strcat('run_',num2str(1),'.nii',volumes)));
    epi_run_2=cellstr(fullfile(PATH_RAW_DATA,subject_list{id},'nii_raw',strcat('run_',num2str(2),'.nii',volumes)));
    epi_run_3=cellstr(fullfile(PATH_RAW_DATA,subject_list{id},'nii_raw',strcat('run_',num2str(3),'.nii',volumes)));
    epi_run_4=cellstr(fullfile(PATH_RAW_DATA,subject_list{id},'nii_raw',strcat('run_',num2str(4),'.nii',volumes)));
    epi_run_5=cellstr(fullfile(PATH_RAW_DATA,subject_list{id},'nii_raw',strcat('run_',num2str(5),'.nii',volumes)));
    epi_run_6=cellstr(fullfile(PATH_RAW_DATA,subject_list{id},'nii_raw',strcat('run_',num2str(6),'.nii',volumes)));
    T1=fullfile(PATH_RAW_DATA,subject_list{id},'nii_raw',[STRING_T1_NAME,',1']);

    %The main body for creating batches
    clear matlabbatch
    matlabbatch{1}.spm.spatial.realign.estwrite.data = {epi_run_1
                                                        epi_run_2
                                                        epi_run_3
                                                        epi_run_4
                                                        epi_run_5
                                                        epi_run_6 }';
    %%
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 0;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 7;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
    
    matlabbatch{2}.spm.temporal.st.scans{1}(1) = cfg_dep('Realign: Estimate & Reslice: Resliced Images (Sess 1)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{1}, '.','rfiles'));
    matlabbatch{2}.spm.temporal.st.scans{2}(1) = cfg_dep('Realign: Estimate & Reslice: Resliced Images (Sess 2)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{2}, '.','rfiles'));
    matlabbatch{2}.spm.temporal.st.scans{3}(1) = cfg_dep('Realign: Estimate & Reslice: Resliced Images (Sess 3)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{3}, '.','rfiles'));
    matlabbatch{2}.spm.temporal.st.scans{4}(1) = cfg_dep('Realign: Estimate & Reslice: Resliced Images (Sess 4)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{4}, '.','rfiles'));
    matlabbatch{2}.spm.temporal.st.scans{5}(1) = cfg_dep('Realign: Estimate & Reslice: Resliced Images (Sess 5)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{5}, '.','rfiles'));
    matlabbatch{2}.spm.temporal.st.scans{6}(1) = cfg_dep('Realign: Estimate & Reslice: Resliced Images (Sess 6)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{6}, '.','rfiles'));
    matlabbatch{2}.spm.temporal.st.nslices = 38;
    matlabbatch{2}.spm.temporal.st.tr = 2;
    matlabbatch{2}.spm.temporal.st.ta = 1.94736842105263;
    matlabbatch{2}.spm.temporal.st.so = [1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 35 37 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38];
    matlabbatch{2}.spm.temporal.st.refslice = 1;
    matlabbatch{2}.spm.temporal.st.prefix = 'a';
    
    matlabbatch{3}.spm.spatial.coreg.estimate.ref(1) = cfg_dep('Slice Timing: Slice Timing Corr. Images (Sess 1)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
    matlabbatch{3}.spm.spatial.coreg.estimate.source = {T1};
    matlabbatch{3}.spm.spatial.coreg.estimate.other(1) = cfg_dep('Slice Timing: Slice Timing Corr. Images (Sess 1)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
    matlabbatch{3}.spm.spatial.coreg.estimate.other(2) = cfg_dep('Slice Timing: Slice Timing Corr. Images (Sess 2)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{2}, '.','files'));
    matlabbatch{3}.spm.spatial.coreg.estimate.other(3) = cfg_dep('Slice Timing: Slice Timing Corr. Images (Sess 3)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{3}, '.','files'));
    matlabbatch{3}.spm.spatial.coreg.estimate.other(4) = cfg_dep('Slice Timing: Slice Timing Corr. Images (Sess 4)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{4}, '.','files'));
    matlabbatch{3}.spm.spatial.coreg.estimate.other(5) = cfg_dep('Slice Timing: Slice Timing Corr. Images (Sess 5)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{5}, '.','files'));
    matlabbatch{3}.spm.spatial.coreg.estimate.other(6) = cfg_dep('Slice Timing: Slice Timing Corr. Images (Sess 6)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{6}, '.','files'));
    matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
    matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
    matlabbatch{4}.spm.spatial.normalise.estwrite.subj.vol = {T1};
    matlabbatch{4}.spm.spatial.normalise.estwrite.subj.resample = {T1};
    matlabbatch{4}.spm.spatial.normalise.estwrite.eoptions.biasreg = 0.0001;
    matlabbatch{4}.spm.spatial.normalise.estwrite.eoptions.biasfwhm = 60;
    matlabbatch{4}.spm.spatial.normalise.estwrite.eoptions.tpm = {'C:\Program Files\MATLAB\R2017a\spm12\tpm\TPM.nii'};
    matlabbatch{4}.spm.spatial.normalise.estwrite.eoptions.affreg = 'mni';
    matlabbatch{4}.spm.spatial.normalise.estwrite.eoptions.reg = [0 0.001 0.5 0.05 0.2];
    matlabbatch{4}.spm.spatial.normalise.estwrite.eoptions.fwhm = 0;
    matlabbatch{4}.spm.spatial.normalise.estwrite.eoptions.samp = 3;
    matlabbatch{4}.spm.spatial.normalise.estwrite.woptions.bb = [-78 -112 -70
                                                                 78 76 85];
    matlabbatch{4}.spm.spatial.normalise.estwrite.woptions.vox = [2 2 2];
    matlabbatch{4}.spm.spatial.normalise.estwrite.woptions.interp = 4;
    matlabbatch{4}.spm.spatial.normalise.estwrite.woptions.prefix = 'w';
    matlabbatch{5}.spm.spatial.normalise.write.subj.def(1) = cfg_dep('Normalise: Estimate & Write: Deformation (Subj 1)', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','def'));
    matlabbatch{5}.spm.spatial.normalise.write.subj.resample(1) = cfg_dep('Coregister: Estimate: Coregistered Images', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','cfiles'));
    matlabbatch{5}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                                                              78 76 85];
    matlabbatch{5}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
    matlabbatch{5}.spm.spatial.normalise.write.woptions.interp = 7;
    matlabbatch{5}.spm.spatial.normalise.write.woptions.prefix = 'w';
    matlabbatch{6}.spm.spatial.smooth.data(1) = cfg_dep('Normalise: Write: Normalised Images (Subj 1)', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
    matlabbatch{6}.spm.spatial.smooth.fwhm = [8 8 8];
    matlabbatch{6}.spm.spatial.smooth.dtype = 0;
    matlabbatch{6}.spm.spatial.smooth.im = 0;
    matlabbatch{6}.spm.spatial.smooth.prefix = 's';
    
    matlabbatch{7}.cfg_basicio.file_dir.file_ops.file_move.files(1) = cfg_dep('Realign: Estimate & Reslice: Realignment Param File (Sess 1)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{1}, '.','rpfile'));
    matlabbatch{7}.cfg_basicio.file_dir.file_ops.file_move.files(2) = cfg_dep('Realign: Estimate & Reslice: Resliced Images (Sess 1)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{1}, '.','rfiles'));
    matlabbatch{7}.cfg_basicio.file_dir.file_ops.file_move.files(3) = cfg_dep('Realign: Estimate & Reslice: Realignment Param File (Sess 2)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{2}, '.','rpfile'));
    matlabbatch{7}.cfg_basicio.file_dir.file_ops.file_move.files(4) = cfg_dep('Realign: Estimate & Reslice: Resliced Images (Sess 2)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{2}, '.','rfiles'));
    matlabbatch{7}.cfg_basicio.file_dir.file_ops.file_move.files(5) = cfg_dep('Realign: Estimate & Reslice: Realignment Param File (Sess 3)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{3}, '.','rpfile'));
    matlabbatch{7}.cfg_basicio.file_dir.file_ops.file_move.files(6) = cfg_dep('Realign: Estimate & Reslice: Resliced Images (Sess 3)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{3}, '.','rfiles'));
    matlabbatch{7}.cfg_basicio.file_dir.file_ops.file_move.files(7) = cfg_dep('Realign: Estimate & Reslice: Realignment Param File (Sess 4)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{4}, '.','rpfile'));
    matlabbatch{7}.cfg_basicio.file_dir.file_ops.file_move.files(8) = cfg_dep('Realign: Estimate & Reslice: Resliced Images (Sess 4)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{4}, '.','rfiles'));
    matlabbatch{7}.cfg_basicio.file_dir.file_ops.file_move.files(9) = cfg_dep('Realign: Estimate & Reslice: Realignment Param File (Sess 5)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{5}, '.','rpfile'));
    matlabbatch{7}.cfg_basicio.file_dir.file_ops.file_move.files(10) = cfg_dep('Realign: Estimate & Reslice: Resliced Images (Sess 5)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{5}, '.','rfiles'));
    matlabbatch{7}.cfg_basicio.file_dir.file_ops.file_move.files(11) = cfg_dep('Realign: Estimate & Reslice: Realignment Param File (Sess 6)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{6}, '.','rpfile'));
    matlabbatch{7}.cfg_basicio.file_dir.file_ops.file_move.files(12) = cfg_dep('Realign: Estimate & Reslice: Resliced Images (Sess 6)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{6}, '.','rfiles'));
    matlabbatch{7}.cfg_basicio.file_dir.file_ops.file_move.files(13) = cfg_dep('Realign: Estimate & Reslice: Mean Image', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rmean'));
    matlabbatch{7}.cfg_basicio.file_dir.file_ops.file_move.files(14) = cfg_dep('Slice Timing: Slice Timing Corr. Images (Sess 1)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
    matlabbatch{7}.cfg_basicio.file_dir.file_ops.file_move.files(15) = cfg_dep('Slice Timing: Slice Timing Corr. Images (Sess 2)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{2}, '.','files'));
    matlabbatch{7}.cfg_basicio.file_dir.file_ops.file_move.files(16) = cfg_dep('Slice Timing: Slice Timing Corr. Images (Sess 3)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{3}, '.','files'));
    matlabbatch{7}.cfg_basicio.file_dir.file_ops.file_move.files(17) = cfg_dep('Slice Timing: Slice Timing Corr. Images (Sess 4)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{4}, '.','files'));
    matlabbatch{7}.cfg_basicio.file_dir.file_ops.file_move.files(18) = cfg_dep('Slice Timing: Slice Timing Corr. Images (Sess 5)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{5}, '.','files'));
    matlabbatch{7}.cfg_basicio.file_dir.file_ops.file_move.files(19) = cfg_dep('Slice Timing: Slice Timing Corr. Images (Sess 6)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{6}, '.','files'));
    matlabbatch{7}.cfg_basicio.file_dir.file_ops.file_move.files(20) = cfg_dep('Coregister: Estimate: Coregistered Images', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','cfiles'));
    matlabbatch{7}.cfg_basicio.file_dir.file_ops.file_move.files(21) = cfg_dep('Normalise: Estimate & Write: Deformation (Subj 1)', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','def'));
    matlabbatch{7}.cfg_basicio.file_dir.file_ops.file_move.files(22) = cfg_dep('Normalise: Estimate & Write: Normalised Images (Subj 1)', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
    matlabbatch{7}.cfg_basicio.file_dir.file_ops.file_move.files(23) = cfg_dep('Normalise: Write: Normalised Images (Subj 1)', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
    matlabbatch{7}.cfg_basicio.file_dir.file_ops.file_move.files(24) = cfg_dep('Smooth: Smoothed Images', substruct('.','val', '{}',{6}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));    
    matlabbatch{7}.cfg_basicio.file_dir.file_ops.file_move.action.moveto = {path_processed_result};
    %Make directory for the preprocessed data
    mkdir(path_processed_result);
    %Save the preprocessing batch
    save(file_save_mat_result,'matlabbatch');
end

%% Run precrocessing
%% NOTE: Need to modify  df1014 (has abnormal amount of runs)
error_preprocessing=cell(length(subject_list),1);
for id=1:numel(subject_list)%(floor(length(subject_list)/2)) %split this into two loops
        try
             clear matlabbatch
             load(strcat(PATH_PROCESSED_MAT,'\pre_',subject_list{id},'.mat'));
             spm_jobman('run',matlabbatch);
             path_processed_result=fullfile(PATH_PROCESSED_DATA,subject_list{id});
             copyfile(fullfile(path_processed_result,'T1.nii'),...
                      fullfile(PATH_RAW_DATA,subject_list{id},'nii_raw'));
             strcat('Finish:',subject_list{id})
        catch ME
             error_preprocessing{id,1}=strcat('Sub',subject_list{id},': ',ME.message)
        end
end
save('D:\Yun-Shiuan_LAMBDA\Adult\preprocessed_data_with_reslice\error_1st_part.mat','error_preprocessing')
%% Create mat for Check Registration 
PATH_CHECK_REG=fullfile(PATH_PROCESSED_DATA,'Check_Reg');
mkdir(PATH_CHECK_REG)
for id=1:numel(subject_list)
    clear matlabbatch
    matlabbatch{1}.spm.util.checkreg.data = {
                                             'C:\Program Files\MATLAB\R2017a\spm12\tpm\TPM.nii,1',...
                                             'C:\Program Files\MATLAB\R2017a\spm12\tpm\TPM.nii,2',...
                                             fullfile(PATH_RAW_DATA,subject_list{id},'\nii_raw\T1_rescaled.nii,1'),...
                                             fullfile(PATH_PROCESSED_DATA,subject_list{id},'\wT1_rescaled.nii,1'),...
                                             'C:\Program Files\MATLAB\R2017a\spm12\tpm\TPM.nii,3',...
                                             fullfile(PATH_PROCESSED_DATA,subject_list{id},'\warrun_6.nii,105'),...
                                             fullfile(PATH_RAW_DATA,subject_list{id},'\nii_raw\run_6.nii,105'),...
                                             fullfile(PATH_PROCESSED_DATA,subject_list{id},'\rrun_6.nii,105')...
                                             }';
     save(fullfile(PATH_CHECK_REG,[subject_list{id},'.mat']),'matlabbatch')
end

%% Convert all motion plots from .txt to .png
% addpath('D:\GoogleDrive\Lambda_code\m_file\LAMBDA\tool_code');
plotMotionFromTxt(PATH_PROCESSED_DATA)