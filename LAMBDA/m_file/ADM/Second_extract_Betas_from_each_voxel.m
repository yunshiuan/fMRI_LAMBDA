%% (Part8-c)Second_level:Extract Betas from each voxel
%addpath('/usr/local/spm12/')
%spm_jobman('initcfg') -- before running the batch
path='/bml/Data/Bank6/ADM-YunShiuan';
cd (strcat(path,'/Second_ROI'));

%read all betas from data(waited to be masked)
imgfilename=spm_select(Inf,'image','Select files...',{},strcat(path,'/First_design_matrix_with_contrastncovariate'),'.*con_0001.nii',1);
imgfilename=imgfilename([1:8,10:16,18:20,22:size(imgfilename,1)],:); %exclude S031,S040, S044
V=spm_vol(imgfilename); %read betas from cons0001.nii --select 40 participants
Y=spm_read_vols(V);    %contain Betas
%size of Y: 79(x)x95(y)x79(z)x40(people)

%define index of mask with voxels interested
V_mask=spm_vol(spm_select);
Y_mask=spm_read_vols(V_mask);
mask_index=find(Y_mask);%index of the mask
[mask_x,mask_y,mask_z]=ind2sub([79,95,79],mask_index);

%extract Beta of each voxel within mask
%for
for i=1:length(mask_index);
 aaa=Y(mask_x(i), mask_y(i), mask_z(i), :);
 squeeze(aaa);
end
%vetorization
bbb = [mask_x(1:3), mask_y(1:3), mask_z(1:3)]';
data = spm_get_data(V, bbb);
%identify region from mni space and then resize to research size

%How to creat mask from mni space?
%addpath('/usr/local/spm12/toolbox/wfu_pickatlas/')
%wfu_pickatlas
%and then save the mask
%Humanatlas
%select:aal
