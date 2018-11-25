%% (Part7-b)Second_level:Writing Peaks as csv. for defining ROI later(Betas:value-->prob)--Conjunciton peaks
%addpath('/usr/local/spm12/')
%spm_jobman('initcfg') -- before running the batch
path='/bml/Data/Bank6/ADM-YunShiuan';
cd (strcat(path,'/Second_ROI_defined_by_value2prob_new_server/conjunction'));

con_index={'prob','mag','pxm',...
    'accept','gain','loss',...
    'reject','ungain','unloss'};
matlabbatch={};
%% conjunction analysis : significant level: 0.005; cluster size=15;
for i=1:length(con_index); % 9 contrasts
    for j=1:2;% intercept*2
        for v=3:6; %Hedonism*2+Security*2 (contrast number 3 to 6)
matlabbatch{1}.spm.stats.results.spmmat = {strcat('/bml/Data/Bank6/ADM-YunShiuan/Second_ROI_defined_by_value2prob_new_server/conjunction/',con_index{i},'/SPM.mat')};
matlabbatch{1}.spm.stats.results.conspec.titlestr = 'abc';
matlabbatch{1}.spm.stats.results.conspec.contrasts = [j v];
matlabbatch{1}.spm.stats.results.conspec.threshdesc = 'none';
matlabbatch{1}.spm.stats.results.conspec.thresh = 0.005;
matlabbatch{1}.spm.stats.results.conspec.extent = 15;
matlabbatch{1}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{1}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{1}.spm.stats.results.units = 1;
matlabbatch{1}.spm.stats.results.export{1}.csv = true;
spm_jobman('run', matlabbatch);
        end;
    end;
end;
