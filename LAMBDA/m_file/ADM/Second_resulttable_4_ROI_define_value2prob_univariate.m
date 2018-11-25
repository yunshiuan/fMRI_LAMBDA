%% (Part7-a)Second_level:Writing Peaks as csv. for defining ROI later(Betas:value-->prob)--Univariate Peaks
% (1)Two values in the same GLM model
% (2)Two values in seperate GLM models
%addpath('/usr/local/spm12/')
%spm_jobman('initcfg') -- before running the batch
path='/bml/Data/Bank6/ADM-YunShiuan';
cd (strcat(path,'/Second_ROI_defined_by_value2prob_new_server/normalise_3x3x3/two_SVS_in_same_GLM'));

con_index={'prob','mag','pxm',...
    'accept','gain','loss',...
    'reject','ungain','unloss','choice'};

%% (1)Two values in the same GLM model
%% significant level: 0.005; cluster size=15;
for i=1:length(con_index); % 10 contrasts
    for j=1:6;% intercept*2+Hedonism*2+Security*2
matlabbatch{1}.spm.stats.results.spmmat = {strcat('/bml/Data/Bank6/ADM-YunShiuan/',...
'Second_ROI_defined_by_value2prob_new_server/normalise_3x3x3/two_SVS_in_same_GLM/',con_index{i},'/SPM.mat')};
matlabbatch{1}.spm.stats.results.conspec.titlestr = '';
matlabbatch{1}.spm.stats.results.conspec.contrasts = j;
matlabbatch{1}.spm.stats.results.conspec.threshdesc = 'none';
matlabbatch{1}.spm.stats.results.conspec.thresh = 0.005;
matlabbatch{1}.spm.stats.results.conspec.extent = 15;
matlabbatch{1}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{1}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{1}.spm.stats.results.units = 1;
matlabbatch{1}.spm.stats.results.export{1}.ps = true;
matlabbatch{1}.spm.stats.results.export{2}.csv = true;
spm_jobman('run', matlabbatch);
    end;
end;

%% (2)Two values in seperate GLM models
cd (strcat(path,'/Second_ROI_defined_by_value2prob_new_server/normalise_3x3x3/two_SVS_in_seperate_GLM/'));
SVS={'hedonism','security'};
%% significant level: 0.005; cluster size=15;
for i=1:length(con_index); % 9 contrasts
    for j=1:4;% intercept*2+Value*2
        for v=1:2; % two values
matlabbatch={};
filelocate=char(strcat(path,'/Second_ROI_defined_by_value2prob_new_server/normalise_3x3x3/two_SVS_in_seperate_GLM/',con_index(i),'/',SVS(v)));            
matlabbatch{1}.spm.stats.results.spmmat = {strcat(filelocate,'/SPM.mat')};
matlabbatch{1}.spm.stats.results.conspec.titlestr = '';
matlabbatch{1}.spm.stats.results.conspec.contrasts = j;
matlabbatch{1}.spm.stats.results.conspec.threshdesc = 'none';
matlabbatch{1}.spm.stats.results.conspec.thresh = 0.005;
matlabbatch{1}.spm.stats.results.conspec.extent = 15;
matlabbatch{1}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{1}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{1}.spm.stats.results.units = 1;
matlabbatch{1}.spm.stats.results.export{1}.ps = true;
matlabbatch{1}.spm.stats.results.export{2}.csv = true;
spm_jobman('run', matlabbatch);
        end;
    end;
end;
