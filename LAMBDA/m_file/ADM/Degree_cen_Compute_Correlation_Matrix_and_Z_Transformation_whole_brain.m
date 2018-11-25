%% Part 6-B: Degree Centrality (whole brain voxel-wise): Compute Correlation Matrix(R),degree of each ROI(D) and perform Z transformation(Z)
%addpath('/usr/local/spm12/')
%spm_jobman('initcfg') -- before running the batch
%% No cerebellum anyways (Not included in the cerebral grey matter mask)
path='/bml/Data/Bank6/ADM-YunShiuan';
cd(path);
addpath(strcat(path,'/Scripts(m)_files'));% to enable using the function "cor2mni"

id=csvread(strcat(path,'/young_idlist.csv')); %read only the young's IDs
young_name=strcat('S0',num2str(id));
young_name=cellstr(young_name); %convert char to cell to enable choosing string by index
young_name_40=cellstr(strvcat(young_name{~ismember(young_name,{'S031' 'S040' 'S044'}')})); %Exclude S031,S040,S044

%% Load in the Residual in the the whole cerebral grey matter of 40 people
cd('/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_3_R_and_Z_Matrix_of_160_ROI/whole_brain_approach/1_raw_residual_whole_brain')
load('beta_whole_cerebral_grey_matter.mat');

cd('/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_3_R_and_Z_Matrix_of_160_ROI/whole_brain_approach/2_R_D_and_Z_matrix')
%% Compute R&D matrix (Don't save R matrix because it is way too big)
%Conver cell to 1090 0bservations x 33221 voxels
% collect_R_matrix={};
collect_D_matrix={};
for id=1:size(beta_whole_brain,1);
% load('beta_whole_cerebral_grey_matter.mat');% Really Bad practice, in order to reduce memory usage
% beta_whole_brain=beta_whole_brain(id,:);
    try
df=cell2mat(beta_whole_brain(id,:)); % 1090 0bservations x 33221 voxels
% collect_R_matrix{1,id}=single(corrcoef(df));
R=corrcoef(df);
D=sum(R>0.25 & R~=1,1);
collect_D_matrix{id,1}=D;
    catch ME
        fprintf(strcat('Error in : ',num2str(id),';', ME.message));
    end;
id
end
collect_D_matrix=cell2mat(collect_D_matrix); % 40id x 33221 ROI matrix
save('collect_D_matrix_whole_brain_id_40.mat','collect_D_matrix','-v7.3');

% %Compute by hand (corrcoef will exceed memory limit)
% df=cell2mat(beta_whole_brain(id,:));
% An=df(:,i);
% Bn=df(:,j);
% An=bsxfun(@minus,A,mean(A,1));
% Bn=bsxfun(@minus,B,mean(B,1));
% An=bsxfun(@times,An,1./sqrt(sum(An.^2,1)));
% Bn=bsxfun(@times,Bn,1./sqrt(sum(Bn.^2,1)));
% C=sum(An.*Bn,1);
% R[i,j]=C;
% R[j,i]=C;

% Check how many positive and negative correlation for each person
% pn_ratio={};
% for id=1:size(beta_160_ROI,1);
% pn_ratio{1,id}=sum(sum(collect_R_matrix{id}>0 & collect_R_matrix{id}~=1))/sum(sum(collect_R_matrix{id}<0));
% end
% mean(cell2mat(pn_ratio))% positive-negative ratio mean=0.95

%% Derive Z matrix (strandardize across 160 ROI within the person)
collect_Z_matrix=zscore(collect_D_matrix,1,2); %Standardise across row(160 ROI)
%mean(collect_Z_matrix,2)
%std(collect_Z_matrix,0,2)
save('collect_Z_matrix_whole_brain_id_40','collect_Z_matrix');

%% Smooth Z matrix before lm with SVS
%% Write Z matrix as nii. first
% Create grey mask from Hvard-Oxford Atlas which will be overlaid with the whole brain of S023
 V_mask=spm_vol('/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_3_R_and_Z_Matrix_of_160_ROI/whole_brain_approach/1_raw_residual_whole_brain/cerebral_grey_mask.nii');% Need transfrom to mni space(because con_0001.nii is in mni space)
 Y_mask=spm_read_vols(V_mask);
 T_mask=V_mask.mat; %obtain the transformation rule to enable transform to mni space
 c_mask=find(Y_mask==1);% select the voxel which mask value==2|13 (Grey cerebral matter)
 [x y z]=ind2sub(size(Y_mask),c_mask);
 Y_mni_mask=cor2mni([x y z],T_mask);%transform to mni space corrdinates
 % Get the MNI coordinates of S023's preprocessed image (only in cerebral grey matter region)
cd('/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_2_partial_out_to_remain_only_residuals/residual_result_after_partialing/S023')
V_template=spm_vol('Res_0001.nii');
Y_template=spm_read_vols(V_template);
T_template=V_template.mat; %obtain the transformation rule to enable transform to mni space
c_template=find(Y_template);% select tall voxels in the template
[x y z]=ind2sub(size(Y_template),c_template);
Y_mni_template=cor2mni([x y z],T_template); %transform to mni space corrdinates
%Check if the coordinate locates inside grey_white mask
Y_template_grey_index=ismember(Y_mni_template,Y_mni_mask,'rows');
Y_mni_template_grey=Y_mni_template(Y_template_grey_index,:);
Y_cor_template=[x y z];
Y_cor_template_grey=Y_cor_template(Y_template_grey_index,:);

%Write nii
cd('/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_3_R_and_Z_Matrix_of_160_ROI/whole_brain_approach/2_R_D_and_Z_matrix/Z_matrix_unsmoothed_nii')
T_template=V_template.mat; %obtain the transformation rule to enable transform to mni space
for id=1:size(young_name_40,1)
img=NaN([53,63,52]);
vol=V_template;
img(sub2ind([53,63,52], Y_cor_template_grey(:,1) , Y_cor_template_grey(:,2), Y_cor_template_grey(:,3)))=collect_Z_matrix(id,:)';
file_name=char(strcat('/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_3_R_and_Z_Matrix_of_160_ROI/whole_brain_approach/2_R_D_and_Z_matrix/',...
    'Z_matrix_unsmoothed_nii/Z_matrix',young_name_40{id},'.nii'));
vol.fname=file_name;
spm_write_vol(vol,img);
end

%% Smooth Z matrix.nii
%%
file_name_list={};
for id=1:size(young_name_40,1)
file_name_list{id,1}=char(strcat('/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_3_R_and_Z_Matrix_of_160_ROI/whole_brain_approach/2_R_D_and_Z_matrix/',...
    'Z_matrix_unsmoothed_nii/Z_matrix',young_name_40{id},'.nii'));
end
matlabbatch{1}.spm.spatial.smooth.data = file_name_list;
%%
matlabbatch{1}.spm.spatial.smooth.fwhm = [8 8 8];
matlabbatch{1}.spm.spatial.smooth.dtype = 0;
matlabbatch{1}.spm.spatial.smooth.im = 0;
matlabbatch{1}.spm.spatial.smooth.prefix = 'smoothed_';
spm_jobman('run',matlabbatch);

%% Read Smoothed nii. in and save as a single .mat file (40 x 35289 matrix)
cd('/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_3_R_and_Z_Matrix_of_160_ROI/whole_brain_approach/2_R_D_and_Z_matrix/Z_matrix_smoothed_nii')
smoothed_file_name_list=char(strcat('/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_3_R_and_Z_Matrix_of_160_ROI/whole_brain_approach/2_R_D_and_Z_matrix/',...
    'Z_matrix_smoothed_nii/smoothed_Z_matrix',young_name_40,'.nii'));
V_smoothed_list=spm_vol(smoothed_file_name_list);
[Y_smoothed_list XYZ]=spm_read_vols(V_smoothed_list);
XYZ=XYZ';
voxel_mni_index=ismember(XYZ,Y_mni_template_grey,'rows');
voxel_mni=XYZ(voxel_mni_index,:);
voxel_cor=mni2cor(voxel_mni,T_template);

collect_Z_matrix_smoothed={};
n=size(Y_mni_template_grey,1);
parfor (id=1:size(young_name_40,1),8)
for c=1:n; 
collect_Z_matrix_smoothed{id,c}=Y_smoothed_list(voxel_cor(c,1),voxel_cor(c,2),voxel_cor(c,3),id);
[id c]
end
end
collect_Z_matrix_smoothed=cell2mat(collect_Z_matrix_smoothed); % 1090 0bservations x 33221 voxels
save(strcat('/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_3_R_and_Z_Matrix_of_160_ROI/whole_brain_approach/2_R_D_and_Z_matrix/',...
'smoothed_collect_Z_matrix_whole_brain_id_40'),'collect_Z_matrix_smoothed');

