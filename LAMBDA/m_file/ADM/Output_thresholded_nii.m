%% Part Appendix-1: Output nii files with k size threshold
%addpath('/usr/local/spm12/')
%spm_jobman('initcfg') -- before running the batch
%% GLM Approach
cd('/bml/Data/Bank6/ADM-YunShiuan/Second_ROI_defined_by_value2prob_new_server/normalise_3x3x3/two_SVS_in_same_GLM')
cond_list={'prob','mag','pxm',...
    'accept','gain','loss',...
    'reject','ungain','unloss','choice'};

for f=1:length(cond_list)
    file_name={strcat('/bml/Data/Bank6/ADM-YunShiuan/Second_ROI_defined_by_value2prob_new_server/normalise_3x3x3/two_SVS_in_same_GLM/',...
        cond_list{f},'/SPM.mat')};
matlabbatch{1}.spm.stats.results.spmmat = file_name;
matlabbatch{1}.spm.stats.results.conspec.titlestr = '';
matlabbatch{1}.spm.stats.results.conspec.contrasts = Inf;
matlabbatch{1}.spm.stats.results.conspec.threshdesc = 'none';
matlabbatch{1}.spm.stats.results.conspec.thresh = 0.005;
matlabbatch{1}.spm.stats.results.conspec.extent = 20;
matlabbatch{1}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{1}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{1}.spm.stats.results.units = 1;
matlabbatch{1}.spm.stats.results.export{1}.tspm.basename = 'k20';
spm_jobman('run',matlabbatch)
end
%% Degree Centrality
cd('/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_4_Regression_degree_cen_by_SVS_by_R/whole_brain_approach/lm_result_standardized_with_SPM')
file_name={'/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_4_Regression_degree_cen_by_SVS_by_R/whole_brain_approach/lm_result_standardized_with_SPM/SPM.mat'};
matlabbatch{1}.spm.stats.results.spmmat = file_name;
matlabbatch{1}.spm.stats.results.conspec.titlestr = '';
matlabbatch{1}.spm.stats.results.conspec.contrasts = Inf;
matlabbatch{1}.spm.stats.results.conspec.threshdesc = 'none';
matlabbatch{1}.spm.stats.results.conspec.thresh = 0.005;
matlabbatch{1}.spm.stats.results.conspec.extent = 15;
matlabbatch{1}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{1}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{1}.spm.stats.results.units = 1;
matlabbatch{1}.spm.stats.results.export{1}.tspm.basename = 'k15';
spm_jobman('run',matlabbatch)
