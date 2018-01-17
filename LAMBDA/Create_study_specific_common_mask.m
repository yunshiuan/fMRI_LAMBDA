%% Create a common mask based on first level results fisrt, to ensure there's signal within the mask for each particiapants

% (Set this as explicit mask for second level analysis)
%Create the common mask (n=40) by ImageCalculator
%NOTE:
%overlapped_WM_GM_common_mask_n40_with_signal 
            %with floating number indicating probabiltity
%binary_overlapped_WM_GM_common_mask_n40_with_signal 
            %only 1 and 0
% Have NOT exclude 'df1021' (the one being partially cut)
masks=strcat('D:\Yun-Shiuan_LAMBDA\First_level_parametric_design_matrix_and_result\',...
             subject_list,'\mask.nii,1');
expression='i1';
for i=2:size(masks,1)
    expression=[expression '.*i' num2str(i)];
end
matlabbatch='';
matlabbatch{1}.spm.util.imcalc.input = masks;
matlabbatch{1}.spm.util.imcalc.output = 'D:\Yun-Shiuan_LAMBDA\template\for_this_study\common_mask_n40_with_signal.nii';
matlabbatch{1}.spm.util.imcalc.outdir = {''};
matlabbatch{1}.spm.util.imcalc.expression = expression;
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = -7;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
%save('D:\Yun-Shiuan_LAMBDA\template\common_mask_n40.mat','matlabbatch');
%spm_jobman('run',matlabbatch);

% Overlap the common mask with grey&white matter
matlabbatch='';
matlabbatch{1}.spm.util.imcalc.input = {
                                        'D:\Yun-Shiuan_LAMBDA\template\for_this_study\common_mask_n40_with_signal.nii,1'
                                        'D:\Yun-Shiuan_LAMBDA\template\mask_grey_white.nii,1'
                                        };
matlabbatch{1}.spm.util.imcalc.output = 'D:\Yun-Shiuan_LAMBDA\template\for_this_study\overlapped_WM_GM_common_mask_n40_with_signal';
matlabbatch{1}.spm.util.imcalc.outdir = {''};
matlabbatch{1}.spm.util.imcalc.expression = 'i1.*i2';
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = -7;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
spm_jobman('run',matlabbatch);
% % Smooth the overlapped mask
% % Do NOT smooth --> It is menaingless if one is going to convert it to binary later.
% matlabbatch='';
% matlabbatch{1}.spm.spatial.smooth.data = {'D:\Yun-Shiuan_LAMBDA\template\for_this_study\overlapped_WM_GM_common_mask_n40_with_signal.nii,1'};
% matlabbatch{1}.spm.spatial.smooth.fwhm = [8 8 8];
% matlabbatch{1}.spm.spatial.smooth.dtype = 0;
% matlabbatch{1}.spm.spatial.smooth.im = 0;
% matlabbatch{1}.spm.spatial.smooth.prefix = 'smooth_';
% spm_jobman('run',matlabbatch);

% Set a cutoff criteria
matlabbatch='';
matlabbatch{1}.spm.util.imcalc.input = {'D:\Yun-Shiuan_LAMBDA\template\for_this_study\overlapped_WM_GM_common_mask_n40_with_signal.nii,1'};
matlabbatch{1}.spm.util.imcalc.output = 'D:\Yun-Shiuan_LAMBDA\template\for_this_study\binary_overlapped_WM_GM_common_mask_n40_with_signal';
matlabbatch{1}.spm.util.imcalc.outdir = {''};
matlabbatch{1}.spm.util.imcalc.expression = 'i1 > 0.65';
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = -7;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
spm_jobman('run',matlabbatch);
