%% (Part8-a)Second_level:Extract Betas from ROI(defined by prob)
%addpath('/usr/local/spm12/')
%spm_jobman('initcfg') -- before running the batch
path='/home/.bml/Data/Bank6/ADM-YunShiuan';
cd (strcat(path,'/Second_ROI'));
%% Betas of Prob. (Extract from ROI of Prob)
% chose 40's con_001s
imgfilename=spm_select(Inf,'image','Select files...',{},strcat(path,'/First_design_matrix_with_contrastncovariate'),'.*con_0001.nii',1);
imgfilename=imgfilename([1:8,10:16,18:20,22:size(imgfilename,1)],:); %exclude S031,S040, S044
%[filename,dirs]=spm_select('List',strcat(path,'/First_design_matrix_with_contrastncovariate'),'.*con_0001.nii'); 
beta_prob_left_striatum = spm_summarise(imgfilename, struct('def', 'sphere', 'spec', 8, 'xyz', [-8, 4, -6]'), @mean);
beta_prob_right_striatum = spm_summarise(imgfilename, struct('def', 'sphere', 'spec', 8, 'xyz', [12, 8, -10]'), @mean);
%% Betas of PxM
% chose 40's con_003s
imgfilename=spm_select(Inf,'image','Select files...',{},strcat(path,'/First_design_matrix_with_contrastncovariate'),'.*con_0003.nii',1);
imgfilename=imgfilename([1:8,10:16,18:20,22:size(imgfilename,1)],:); %exclude S031,S040, S044
%[filename,dirs]=spm_select('List',strcat(path,'/First_design_matrix_with_contrastncovariate'),'.*con_0001.nii'); 
beta_pxm_right_striatum = spm_summarise(imgfilename, struct('def', 'sphere', 'spec', 8, 'xyz', [8, 4, 0]'), @mean);
%% Betas of Gain.(Extract from ROI of Gain)
% chose 40's con_005s
imgfilename=spm_select(Inf,'image','Select files...',{},strcat(path,'/First_design_matrix_with_contrastncovariate'),'.*con_0005.nii',1);
imgfilename=imgfilename([1:8,10:16,18:20,22:size(imgfilename,1)],:); %exclude S031,S040, S044
%[filename,dirs]=spm_select('List',strcat(path,'/First_design_matrix_with_contrastncovariate'),'.*con_0001.nii'); 
beta_gain_left_striatum = spm_summarise(imgfilename, struct('def', 'sphere', 'spec', 8, 'xyz', [-16, 6, 0]'), @mean);
beta_gain_right_striatum = spm_summarise(imgfilename, struct('def', 'sphere', 'spec', 8, 'xyz', [14, -8, 4]'), @mean);

%% Output as csv
dlmwrite('betas_Prob_lr_PxM_r_Gain_lr.csv',...
[beta_prob_left_striatum,beta_prob_right_striatum,beta_pxm_right_striatum...
beta_gain_left_striatum,beta_gain_right_striatum],'precision',9);

