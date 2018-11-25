%% Part 5-B: Degree Centrality: Extract Residuals for whole brian (grey matter) voxel-wise
%addpath('/usr/local/spm12/')
%spm_jobman('initcfg') -- before running the batch

path='/bml/Data/Bank6/ADM-YunShiuan';
cd(path);
addpath(strcat(path,'/Scripts(m)_files'));% to enable using the function "cor2mni"

id=csvread(strcat(path,'/young_idlist.csv')); %read only the young's IDs
young_name=strcat('S0',num2str(id));
young_name=cellstr(young_name); %convert char to cell to enable choosing string by index
young_name_40=cellstr(strvcat(young_name{~ismember(young_name,{'S031' 'S040' 'S044'}')})); %Exclude S031,S040,S044

%% Load in the MNI coordinatess of the S023 to get MNI indexes of grey matter
cd('/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_2_partial_out_to_remain_only_residuals/residual_result_after_partialing')
%Load in mask of white and grey matter for overlaying latter
%% Create grey mask from Hvard-Oxford Atlas which will be overlaid with the whole brain of S023
 V_mask=spm_vol(strcat(path,'/HarvardOxford-sub-maxprob-thr50-1mm.nii'));% Need transfrom to mni space(because con_0001.nii is in mni space)
 Y_mask=spm_read_vols(V_mask);
 T_mask=V_mask.mat; %obtain the transformation rule to enable transform to mni space
 c_mask=find(Y_mask~=0 & Y_mask~=1 & Y_mask~=12 & Y_mask~=3 & Y_mask~=14 & Y_mask~=8);% select the voxel which mask value==2|13 (Grey cerebral matter)
 [x y z]=ind2sub(size(Y_mask),c_mask);
 Y_mni_mask=cor2mni([x y z],T_mask);%transform to mni space corrdinates
 %% Get the MNI coordinates of S023's preprocessed image (only in cerebral grey matter region)
cd('/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_2_partial_out_to_remain_only_residuals/residual_result_after_partialing/S023')
V_template=spm_vol('Res_0001.nii');
Y_template=spm_read_vols(V_template);
T_template=V_template.mat; %obtain the transformation rule to enable transform to mni space
c_template=find(Y_template);% select tall voxels in the template
[x y z]=ind2sub(size(Y_template),c_template);
Y_mni_template=cor2mni([x y z],T_template); %transform to mni space corrdinates
%Check if the coordinate locates inside grey_white mask
Y_mni_template_grey_index=ismember(Y_mni_template,Y_mni_mask,'rows');
Y_mni_template_grey=Y_mni_template(Y_mni_template_grey_index,:);

%Write the image to display
% vol=V_template;
% file_name='/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_3_R_and_Z_Matrix_of_160_ROI/whole_brain_approach/1_raw_residual_whole_brain/cerebral_grey_mask.nii';
% vol.fname=file_name;
% %convert Y_mni_mask to matrix form(voxle-based coordinates)
% img=NaN([53,63,52]);
% Y_cor=mni2cor([round(Y_mni_template_grey(:,1)),round(Y_mni_template_grey(:,2)),...
%     round(Y_mni_template_grey(:,3))],T_template);
% img(sub2ind([53,63,52], Y_cor(:,1) , Y_cor(:,2), Y_cor(:,3)))=1;
% spm_write_vol(vol,img);

%% Collect Residuals in the cerebral grey matter of 40 people
% Use the below section instead.
% This section of code is way too time-consuming!
% beta_whole_brain={};
% n=size(Y_mni_template_grey,1);
% parfor (id=1:size(young_name_40,1),16);
% [imgfilename,dirs]=spm_select('ExtFPListRec',...%FORMAT [files,dirs] = spm_select('ExtFPListRec',direc,filt,frames)
%            strcat(path,'/degree_centrality/Part_2_partial_out_to_remain_only_residuals/',...
%            'residual_result_after_partialing/',young_name_40{id}),...
%            'Res_\d+',1);%In MNI space (already normalize)
% imgfilename=cellstr(imgfilename);       
% for c=1:n; %use the constant "n" to enable parfor to work
% %% Extract voxel signals from the overlaid ROI
% voxel_index=Y_mni_template_grey(c,:);
% beta_whole_brain{id,c}=... %for all cerebal grey matter voxels
%        spm_summarise(imgfilename,voxel_index',@mean);
% [id c]
% end
% end




%% Another Approach (Don't use spm_summarise. Instead, directly index the 4D matrix)
%% Using spm_summarise is way too time-consuming...
% Read in all 3D matrix per person and convert to a 4D matrix
beta_whole_brain={};
n=size(Y_mni_template_grey,1);
parfor (id=1:size(young_name_40,1),16);
[imgfilename,dirs]=spm_select('ExtFPListRec',...%FORMAT [files,dirs] = spm_select('ExtFPListRec',direc,filt,frames)
           strcat(path,'/degree_centrality/Part_2_partial_out_to_remain_only_residuals/',...
           'residual_result_after_partialing/',young_name_40{id}),...
           'Res_\d+',1);
V_1090=spm_vol(imgfilename);
[Y_1090 XYZ]=spm_read_vols(V_1090);
XYZ=XYZ';
voxel_mni_index=ismember(XYZ,Y_mni_template_grey,'rows');
voxel_mni=XYZ(voxel_mni_index,:);
voxel_cor=mni2cor(voxel_mni,T_template);
for c=1:n; %use the constant "n" to enable parfor to work
beta_whole_brain{id,c}=[reshape(Y_1090(voxel_cor(c,1),voxel_cor(c,2),voxel_cor(c,3),1:size(Y_1090,4)),1,size(Y_1090,4))]';
[id c]
end
end
%Save As .mat 
save('beta_whole_cerebral_grey_matter.mat','beta_whole_brain','-v7.3');
