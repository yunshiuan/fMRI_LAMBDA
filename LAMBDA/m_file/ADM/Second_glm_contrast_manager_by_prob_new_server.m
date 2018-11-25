%% (Part6-c)Second_level: Second Level GLM Excluded SVS(see the pure effect of prob etc.)
%UNUSED: significance of cluster is not higher
%addpath('/usr/local/spm12/')
%spm_jobman('initcfg') -- before running the batch
path='/bml/Data/Bank6/ADM-YunShiuan';
cd (strcat(path,'/Second_ROI_defined_by_prob_new_server'));

%% Get Young ID, gender, age, and Values
%id
%id=csvread(strcat(path,'/young_idlist.csv')); %read only the young's IDs
%young_name=strcat('S0',num2str(id));
%young_name=cellstr(young_name); %convert char to cell to enable choosing string by index
info=load(strcat(path,'/id_age_gender_sec_hed_sti.mat'));
index=logical(cell2mat(cellfun(@(x) ~(isequal(x,'S031')|isequal(x,'S040')|isequal(x,'S044')),info.a(:,1),'un',0)));%exclude i=(S031)(S040)21(S044) who have excessive motions

info.a=info.a(index,:);
%id
young_name=info.a(:,1);
%age
age=info.a(:,2);
%gender(1=M;0=F)
gender=info.a(:,3);
%security
security=info.a(:,4);
%hedonism
hedonism=info.a(:,5);
%stimulation
stimulation=info.a(:,6);


%% Gather contrasts (Prob,Mag,PxM...)
con_index={'prob','mag','pxm',...
    'accept','gain','loss',...
    'reject','ungain','unloss'};
contrast=cellfun(@(x) strcat('/bml/Data/Bank6/ADM-YunShiuan/First_design_matrix_with_contrastncovariate/',...
    young_name(1:end),'/con_000',num2str(x),'.nii,1'),{1,2,3,4,5,6,7,8,9},'un',0);

%% Edit Batch info (Prob,Mag,PxM,Accept,Gain,Loss,Reject,Ungain,Unloss)
for i=1:length(con_index);
load(strcat(path,'/Second_ROI_defined_by_prob_new_server/Second_ROI_template_value2prob.mat'));
%%directory which files will be created
matlabbatch{1}.spm.stats.factorial_design.dir = ...
    {char(strcat('/bml/Data/Bank6/ADM-YunShiuan/Second_ROI_defined_by_prob_new_server/',con_index(i)))};
%Scans(contrasts) to be tested
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans=contrast{i};

%Covariate
%Gender
matlabbatch{1}.spm.stats.factorial_design.cov(1).c = str2num(cell2mat(gender)); %vector of each
matlabbatch{1}.spm.stats.factorial_design.cov(1).cname = 'Gender';
%Age(Excluded)
%matlabbatch{1}.spm.stats.factorial_design.cov(2).c = str2num(cell2mat(age)); %vector
%matlabbatch{1}.spm.stats.factorial_design.cov(2).cname = 'Age';
%Security
% matlabbatch{1}.spm.stats.factorial_design.cov(2).c = str2num(cell2mat(hedonism));
% matlabbatch{1}.spm.stats.factorial_design.cov(2).cname = 'Hedonism';
% %Hedonism
% matlabbatch{1}.spm.stats.factorial_design.cov(3).c = str2num(cell2mat(security));
% matlabbatch{1}.spm.stats.factorial_design.cov(3).cname = 'Security';

% explicit mask
matlabbatch{1}.spm.stats.factorial_design.masking.em = {'/bml/Data/Bank6/ADM-YunShiuan/mask_grey_white.nii,1'};

%Contrast Manager
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name=strcat(con_index{i},'_positive');
matlabbatch{3}.spm.stats.con.consess{2}.tcon.name=strcat(con_index{i},'_negative');
matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights=[-1];
% matlabbatch{3}.spm.stats.con.consess{3}.tcon.name='hedonism_positive';
% matlabbatch{3}.spm.stats.con.consess{3}.tcon.weights=[0 0 1];
% matlabbatch{3}.spm.stats.con.consess{4}.tcon.name='hedonism_negative';
% matlabbatch{3}.spm.stats.con.consess{4}.tcon.weights=[0 0 -1];
% matlabbatch{3}.spm.stats.con.consess{5}.tcon.name='security_positive';
% matlabbatch{3}.spm.stats.con.consess{5}.tcon.weights=[0 0 0 1];
% matlabbatch{3}.spm.stats.con.consess{6}.tcon.name='security_negative';
% matlabbatch{3}.spm.stats.con.consess{6}.tcon.weights=[0 0 0 -1];


%Save matlabbatch.mat (the file which could be open by GUI)
save(char(strcat(path,'/Second_ROI_defined_by_value2prob_new_server/',con_index(i),'/Second_glm')),'matlabbatch');

%Create "SPM.mat" 
spm_jobman('run', matlabbatch);
end
