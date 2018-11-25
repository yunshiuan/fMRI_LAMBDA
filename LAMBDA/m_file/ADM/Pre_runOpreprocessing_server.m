%% (Part2)Run all the batches!
%spm_jobman('initcfg') -- before running the batch
%mat = ls(strcat(path,'/pre_mat'))%Failed
path='/bml/Data/Bank6/ADM-YunShiuan';
cd (strcat(path,'/pre_mat_3x3x3'));;% Select 3x3x3 (normarlized voxel) version

id=csvread(strcat(path,'/young_idlist.csv')); %read only the young's IDs
young_name=strcat('pre2_',num2str(id),'.mat');
young_name=cellstr(young_name); %convert char to cell to enable choosing string by index

%Save 43 participants brain data to a big cell. (Each row per participant.)

for i = 1:length(young_name); 
    matlabbatch_all{i,1}=load(young_name{i,1}); 
end
%addpath('/usr/local/spm12/') -- before running via SPM
%spm
for i=1:43
spm_jobman('run', matlabbatch_all{i,1}.matlabbatch)
end