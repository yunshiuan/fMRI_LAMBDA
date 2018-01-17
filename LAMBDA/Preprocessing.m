%% Part 1- Preprocessing
%spm_jobman('initcfg') -- before running the batch
path='D:\Yun-Shiuan_LAMBDA\raw_data';
cd(path);
subject_list=cellstr(ls('*')); % subject list
subject_list=subject_list(~cellfun(@isempty,regexp(subject_list,'(?<=df)\d+','match')));% include only subject name
subject_list=subject_list(~ismember(subject_list,{'df1001','df1054'})); %Exclude df1001(no T1) and df1054(incomplete)
%(note: df1014 & de1018 has abnormal amount of runs)

%% Rescaling T1 image by deviding 100,000 (Still not sure if this is appropriate)
for i=1:length(subject_list) %Error:df 1005(diff. data type)%df 1026(no motion corrected T1)
        if ismember(subject_list{i},{'df1026'})
        warning(strcat(subject_list{i},' does not have T1 image.'))
    else
matlabbatch='';
matlabbatch{1}.spm.util.imcalc.input = {char(strcat('D:\Yun-Shiuan_LAMBDA\raw_data\',subject_list{i},'\nii_raw\T1.nii,1'))};
matlabbatch{1}.spm.util.imcalc.output = 'T1_rescaled';
matlabbatch{1}.spm.util.imcalc.outdir = {char(strcat('D:\Yun-Shiuan_LAMBDA\raw_data\',subject_list{i},'\nii_raw'))};
matlabbatch{1}.spm.util.imcalc.expression = 'i1/100000';
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = -7;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
spm_jobman('run',matlabbatch);
        end
end

%% Preprocessing: Create .mat file for preprocessing
for i=2:length(subject_list);
    if ismember(subject_list{i},{'df1026'})
        warning(strcat(subject_list{i},' should be treated seperately, due to missing T1.'))
    else
%Load template
load('D:\Yun-Shiuan_LAMBDA\preprocessing_mat\without_resliced\template_lambda_estimate_no_reslice.mat');
%Replace EPI info.
volumes=strtrim(cellstr(num2str([6:110]', ',%d')));
epi_run_1=cellstr(strcat(path,'\',subject_list{i},'\nii_raw\run_',num2str(1),'.nii',volumes));
epi_run_2=cellstr(strcat(path,'\',subject_list{i},'\nii_raw\run_',num2str(2),'.nii',volumes));
epi_run_3=cellstr(strcat(path,'\',subject_list{i},'\nii_raw\run_',num2str(3),'.nii',volumes));
epi_run_4=cellstr(strcat(path,'\',subject_list{i},'\nii_raw\run_',num2str(4),'.nii',volumes));
epi_run_5=cellstr(strcat(path,'\',subject_list{i},'\nii_raw\run_',num2str(5),'.nii',volumes));
epi_run_6=cellstr(strcat(path,'\',subject_list{i},'\nii_raw\run_',num2str(6),'.nii',volumes));
T1=strcat('D:\Yun-Shiuan_LAMBDA\raw_data\',subject_list{i},'\nii_raw\T1_rescaled.nii,1');
matlabbatch{1}.spm.temporal.st.scans = {epi_run_1,epi_run_2,epi_run_3,epi_run_4,epi_run_5,epi_run_6};

matlabbatch{3}.spm.spatial.coreg.estimate.source = {T1};
matlabbatch{4}.spm.spatial.normalise.estwrite.subj.vol = {T1};
matlabbatch{4}.spm.spatial.normalise.estwrite.subj.resample = {T1};

processed_result_dir=strcat('D:\Yun-Shiuan_LAMBDA\preprocessed_data_without_reslice\',subject_list{i});
matlabbatch{7}.cfg_basicio.file_dir.file_ops.file_move.action.moveto = {processed_result_dir};
save(strcat('D:\Yun-Shiuan_LAMBDA\preprocessing_mat\without_resliced\pre_',subject_list{i},'.mat'),'matlabbatch');
    end
end


%% Run precrocessing
%% NOTE: Need to modify 'df1004'
error_preprocessing=cell(length(subject_list),1);
for i=31:length(subject_list) %split this into two loops
     if ismember(subject_list{i},{'df1026'})
        warning(strcat(subject_list{i},' does not have T1 image.'))
     else
        try
         load(strcat('D:\Yun-Shiuan_LAMBDA\preprocessing_mat\without_resliced\pre_',subject_list{i},'.mat'));
         spm_jobman('run',matlabbatch);
         processed_result_dir=strcat('D:\Yun-Shiuan_LAMBDA\preprocessed_data_without_reslice\',subject_list{i});
         copyfile(strcat( processed_result_dir,'\T1_rescaled.nii'),strcat('D:\Yun-Shiuan_LAMBDA\raw_data\',subject_list{i},'\nii_raw'))
         strcat('Finish:',num2str(i))

        catch ME
            error_preprocessing{i,1}=strcat('Sub',subject_list{i},': ',ME.message)
        end
      end
end

%% Create mat for Check Registration 
for i=1:length(subject_list);
matlabbatch='';
   if ismember(subject_list{i},{'df1026'})
        warning(strcat(subject_list{i},' does not have T1 image.'))
   else
        matlabbatch{1}.spm.util.checkreg.data = {
                                                 char(strcat('C:\Program Files\MATLAB\R2017a\spm12\tpm\TPM.nii,1')),...
                                                 char(strcat('C:\Program Files\MATLAB\R2017a\spm12\tpm\TPM.nii,2')),...
                                                 char(strcat('D:\Yun-Shiuan_LAMBDA\preprocessed_data_without_reslice\',subject_list{i},'\T1_rescaled.nii,1')),...
                                                 char(strcat('D:\Yun-Shiuan_LAMBDA\preprocessed_data_without_reslice\',subject_list{i},'\wT1_rescaled.nii,1')),...
                                                 char(strcat('C:\Program Files\MATLAB\R2017a\spm12\tpm\TPM.nii,3')),...
                                                 char(strcat('D:\Yun-Shiuan_LAMBDA\preprocessed_data_without_reslice\',subject_list{i},'\warun_6.nii,105')),...
                                                 char(strcat('D:\Yun-Shiuan_LAMBDA\raw_data\',subject_list{i},'\nii_raw\run_6.nii,105')),...
                                                 char(strcat('D:\Yun-Shiuan_LAMBDA\preprocessed_data_without_reslice\',subject_list{i},'\arun_6.nii,105'))...
                                                 }';
         save(strcat('D:\Yun-Shiuan_LAMBDA\preprocessed_data_without_reslice\Check_Reg\',subject_list{i},'.mat'),...
             'matlabbatch')
   end
end

%% Convert all motion plots from .txt to .png
addpath('D:\GoogleDrive\Lambda_code\m_file\LAMBDA\tool_code');
% plotMotionFromTxt('D:\Yun-Shiuan_LAMBDA\preprocessed_data_without_reslice')