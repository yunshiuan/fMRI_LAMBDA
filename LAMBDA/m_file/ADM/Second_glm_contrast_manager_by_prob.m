%% (Part6-a)Second_level:Define ROI where Prob/Mag.etc. induce BOLD signal
%addpath('/usr/local/spm12/')
%spm_jobman('initcfg') -- before running the batch
path='/home/.bml/Data/Bank6/ADM-YunShiuan';
cd (strcat(path,'/Second_ROI_defined_by_prob'));

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
%gender
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
contrast=cellfun(@(x) strcat('/home/.bml/Data/Bank6/ADM-YunShiuan/First_design_matrix_with_contrastncovariate/',...
    young_name(1:end),'/con_000',num2str(x),'.nii,1'),{1,2,3,4,5,6,7,8,9},'un',0);

%% Edit Batch info (Prob,Mag,PxM,Accept,Gain,Loss,Reject,Ungain,Unloss)
for i=1:length(con_index);
load('Second_ROI_template.mat');
%%directory which files will be created
matlabbatch{1}.spm.stats.factorial_design.dir = ...
    {char(strcat('/home/.bml/Data/Bank6/ADM-YunShiuan/Second_ROI_defined_by_prob/',con_index(i)))};
%Scans(contrasts) to be tested
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans=contrast{i};

%Covariate
%Gender
matlabbatch{1}.spm.stats.factorial_design.cov(1).c = str2num(cell2mat(gender)); %vector of each
matlabbatch{1}.spm.stats.factorial_design.cov(1).cname = 'Gender';
%Age(Excluded)
%matlabbatch{1}.spm.stats.factorial_design.cov(2).c = str2num(cell2mat(age)); %vector
%matlabbatch{1}.spm.stats.factorial_design.cov(2).cname = 'Age';

% explicit mask
matlabbatch{1}.spm.stats.factorial_design.masking.em = {'/home/.bml/Data/Bank6/ADM-YunShiuan/mask_grey_white.nii,1'};

%Contrast Manager
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name=con_index{i};
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name=strcat(con_index{i},'_positive');
matlabbatch{3}.spm.stats.con.consess{2}.tcon.name=strcat(con_index{i},'_negative');
matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights=[-1];
%Save matlabbatch.mat (the file which could be open by GUI)
save(char(strcat(path,'/Second_ROI_defined_by_prob/',con_index(i),'/Second_ROI')),'matlabbatch');

%Create "SPM.mat" 
spm_jobman('run', matlabbatch);
end