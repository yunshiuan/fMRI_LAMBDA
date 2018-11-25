%% (Part3)First_level:Create Design Matrix for Each person!
%% This Scripts To be ignored! ---> Skip directly to Part5!!
%addpath('/usr/local/spm12/')
%spm_jobman('initcfg') -- before running the batch
%mat = ls(strcat(path,'/pre_mat'))%Failed
path='/home/.bml/Data/Bank6/ADM-YunShiuan';
cd (strcat(path,'/first_design_matrix'));

%% Load template
load('./S023/template_creating_design_matrix.mat');
template=matlabbatch;
scan=num2cell(1:218);
scan=cellfun(@(x) num2str(x),scan,'un',0);
%% Set Onset, Prob., and Mag. to Correct numbers
id=csvread(strcat(path,'/young_idlist.csv')); %read only the young's IDs
young_name=strcat('S0',num2str(id));
young_name=cellstr(young_name); %convert char to cell to enable choosing string by index

for i=1:length(young_name);
    %% Read the trials information of each person
        A2={0,0,0,0,0};% save info to "A2"
for j=1:5;
    fid = fopen(char(strcat(path,'/Run_info(csv)/',young_name(i),'/r',num2str(j),'.csv')));  %open file
 %read all contents into data as a char array 
    A = textscan(fid,'%s','Delimiter', ',');
    A2{1,j}=reshape(A{1},8,90)';
end
%% Derive Onset, Prob., and Mag.
prob_1=str2double(A2{1,1}(1:2:length(A2{1,1}),3));
prob_2=str2double(A2{1,2}(1:2:length(A2{1,1}),3));
prob_3=str2double(A2{1,3}(1:2:length(A2{1,1}),3));
prob_4=str2double(A2{1,4}(1:2:length(A2{1,1}),3));
prob_5=str2double(A2{1,5}(1:2:length(A2{1,1}),3));
prob_mean=mean(mean([prob_1,prob_2,prob_3,prob_4,prob_5]));
prob_center=cellfun(@(x) (x-prob_mean),{prob_1,prob_2,prob_3,prob_4,prob_5},'un',0);

mag_1=str2double(A2{1,1}(1:2:length(A2{1,1}),4));
mag_2=str2double(A2{1,2}(1:2:length(A2{1,1}),4));
mag_3=str2double(A2{1,3}(1:2:length(A2{1,1}),4));
mag_4=str2double(A2{1,4}(1:2:length(A2{1,1}),4));
mag_5=str2double(A2{1,5}(1:2:length(A2{1,1}),4));
mag_mean=mean(mean([mag_1,mag_2,mag_3,mag_4,mag_5]));
mag_center=cellfun(@(x) (x-mag_mean),{mag_1,mag_2,mag_3,mag_4,mag_5},'un',0);

probxmag=cell2mat(prob_center).*cell2mat(mag_center);
probxmag=mat2cell(probxmag,45,[1 1 1 1 1]);

onset_1=str2double(A2{1,1}(1:2:length(A2{1,1}),1));
onset_2=str2double(A2{1,2}(1:2:length(A2{1,1}),1));
onset_3=str2double(A2{1,3}(1:2:length(A2{1,1}),1));
onset_4=str2double(A2{1,4}(1:2:length(A2{1,1}),1));
onset_5=str2double(A2{1,5}(1:2:length(A2{1,1}),1));
onset=cellfun(@(x) (x./1000),{onset_1,onset_2,onset_3,onset_4,onset_5},'un',0);

%% Substitute Onset/Prob/Mag to the template
template{1}.spm.stats.fmri_spec.dir = {char(strcat('/home/.bml/Data/Bank6/ADM-YunShiuan/first_design_matrix/',young_name(i)))};
for s=1:5;
template{1}.spm.stats.fmri_spec.sess(s).scans = (strcat('/home/.bml/Data/Bank6/ADM-YunShiuan/IMG_nii/',young_name(i),'/swadm',num2str(s),'.nii,',scan))';
template{1}.spm.stats.fmri_spec.sess(s).cond.onset = [onset{1,s}];
template{1}.spm.stats.fmri_spec.sess(s).cond.pmod(1).name = 'prob';
template{1}.spm.stats.fmri_spec.sess(s).cond.pmod(1).param = [prob_center{1,s}];
template{1}.spm.stats.fmri_spec.sess(s).cond.pmod(2).name = 'mag';
template{1}.spm.stats.fmri_spec.sess(s).cond.pmod(2).param = [mag_center{1,s}];
template{1}.spm.stats.fmri_spec.sess(s).cond.pmod(3).name = 'probxmag';
template{1}.spm.stats.fmri_spec.sess(s).cond.pmod(3).param = [probxmag{1,s}];
end
%Save .mat
matlabbatch=template;
%Save matlabbatch.mat (the file which could be open by GUI)
save(char(strcat(path,'/first_design_matrix/',young_name(i),'/design_matrix_',num2str(id(i)))),'matlabbatch');
 %Create "SPM.mat" (the design matrix which could be checked)
spm_jobman('run', matlabbatch);
end