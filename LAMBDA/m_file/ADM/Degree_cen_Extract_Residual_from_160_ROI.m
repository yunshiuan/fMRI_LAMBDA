%% Part 5-A: Degree Centrality: Extract Residuals from 160 ROIs
%addpath('/usr/local/spm12/')
%spm_jobman('initcfg') -- before running the batch
path='/bml/Data/Bank6/ADM-YunShiuan';
cd(path);
addpath(strcat(path,'/Scripts(m)_files'));% to enable using the function "cor2mni"

id=csvread(strcat(path,'/young_idlist.csv')); %read only the young's IDs
young_name=strcat('S0',num2str(id));
young_name=cellstr(young_name); %convert char to cell to enable choosing string by index
young_name_40=cellstr(strvcat(young_name{~ismember(young_name,{'S031' 'S040' 'S044'}')})); %Exclude S031,S040,S044

%% Load in the MNI coordinatess of 160 ROI and Create the 8x8x8 Sphere
cd('/bml/Data/Bank6/ADM-YunShiuan/degree_centrality')

%Load in mask of white and grey matter for overlaying latter
%% Create grey_white mask which will be overlaid with ROI to create a combined ROI for extracting Beta later
 V_mask=spm_vol(strcat(path,'/smooth_mask_grey_white.nii'));% Need transfrom to mni space(because con_0001.nii is in mni space)
 Y_mask=spm_read_vols(V_mask);
 Y_mask=Y_mask>0.8;%Set criteria >0.8 as brain(white/grey) region
 T_mask=V_mask.mat; %abtain the transformation rule to enable transform to mni space
 c_mask=find(Y_mask==1);% select the voxel which mask value==1
 [x y z]=ind2sub(size(Y_mask),c_mask);
 Y_mni_mask=cor2mni([x y z],T_mask);%transform to mni space corrdinates
 %% Get the MNI coordinates and make shpere around it
load('Dosenbach_Science_160ROIs_Center.mat')
% Get all coordinate 
ROI_160_center=Dosenbach_Science_160ROIs_Center;
%Add ID for each center(1~160)
ROI_160_center=[[1:size(ROI_160_center,1)]' ROI_160_center];

%% Collect Residuals in the 160 ROI of 40 people
beta_160_ROI={};
n=size(ROI_160_center,1);
parfor (id=1:size(young_name_40,1),16);
[imgfilename,dirs]=spm_select('ExtFPListRec',...%FORMAT [files,dirs] = spm_select('ExtFPListRec',direc,filt,frames)
           strcat(path,'/degree_centrality/Part_2_partial_out_to_remain_only_residuals/',...
           'residual_result_after_partialing/',young_name_40{id}),...
           'Res_\d+',1);%In MNI space (already normalize)
imgfilename=cellstr(imgfilename);       

for c=1:n; %use the constant "n" to enable parfor to work
%% Make Sphere
center=ROI_160_center(c,2:4);
radius=5;
Dist_sq = (bsxfun(@minus,center(1,1),Y_mni_mask(:,1)).^2+...
    bsxfun(@minus,center(1,2),Y_mni_mask(:,2)).^2+...
    bsxfun(@minus,center(1,3),Y_mni_mask(:,3)).^2);
sphere_index=(find(Dist_sq<radius^2));%indexes for the coordinates in the sphere

%% The Critical Overlaying Step
ROI_overlaid=Y_mni_mask(sphere_index,:);  %Extract the coordinates in shpere and mask 

%% Extract voxel signals from the overlaid ROI
ROI_overlaid=round(ROI_overlaid);
beta_160_ROI{id,c}=... %160 ROI
       spm_summarise(imgfilename,ROI_overlaid',@mean);
[id c]
end
end

%Save As .mat and c.sv
save('beta_160_ROI.mat','beta_160_ROI');





%% Save the 160 spheres for display.      
for c=1:size(ROI_160_center,1);
%% Make Sphere
center=ROI_160_center(c,2:4);
radius=5;
Dist_sq = (bsxfun(@minus,center(1,1),Y_mni_mask(:,1)).^2+...
    bsxfun(@minus,center(1,2),Y_mni_mask(:,2)).^2+...
    bsxfun(@minus,center(1,3),Y_mni_mask(:,3)).^2);
sphere_index=(find(Dist_sq<radius^2));%indexes for the coordinates in the sphere

%% The Critical Overlaying Step
ROI_overlaid=Y_mni_mask(sphere_index,:);  %Extract the coordinates in shpere and mask 
vol=V_mask;
file_name=char(strcat(path,'/degree_centrality/Part_3_R_and_Z_Matrix_of_160_ROI/raw_residual_160_ROI/160_ROI_nii/roi_',...
    num2str(center),'.nii'));
vol.fname=file_name;
%convert Y_mni_mask to matrix form(voxle-based coordinates)
img=NaN([121,145,121]);
Y_cor_roi=mni2cor([round(ROI_overlaid(:,1)),round(ROI_overlaid(:,2)), round(ROI_overlaid(:,3))],T_mask);
img(sub2ind([121,145,121], Y_cor_roi(:,1) , Y_cor_roi(:,2), Y_cor_roi(:,3)))=1;
spm_write_vol(vol,img);
c
end

%% Write a single nii with all 160 ROI
ROI_overlaid_all=[];
for c=1:size(ROI_160_center,1);
%% Make Sphere
center=ROI_160_center(c,2:4);
radius=5;
Dist_sq = (bsxfun(@minus,center(1,1),Y_mni_mask(:,1)).^2+...
    bsxfun(@minus,center(1,2),Y_mni_mask(:,2)).^2+...
    bsxfun(@minus,center(1,3),Y_mni_mask(:,3)).^2);
sphere_index=(find(Dist_sq<radius^2));%indexes for the coordinates in the sphere
%% The Critical Overlaying Step
ROI_overlaid=Y_mni_mask(sphere_index,:);  %Extract the coordinates in shpere and mask 
ROI_overlaid_all=[ROI_overlaid_all;ROI_overlaid];%Collect 160 ROI indexes
c
end
vol=V_mask;
file_name=char(strcat(path,'/degree_centrality/Part_3_R_and_Z_Matrix_of_160_ROI/raw_residual_160_ROI/160_ROI_nii/',...
    'all_160_ROI.nii'));
vol.fname=file_name;
%convert Y_mni_mask to matrix form(voxle-based coordinates)
img=NaN([121,145,121]);
Y_cor_roi=mni2cor([round(ROI_overlaid_all(:,1)),round(ROI_overlaid_all(:,2)), round(ROI_overlaid_all(:,3))],T_mask);
img(sub2ind([121,145,121], Y_cor_roi(:,1) , Y_cor_roi(:,2), Y_cor_roi(:,3)))=1;
spm_write_vol(vol,img);


%No Need to write as csv file, since R Could read .mat anyways.

% dlmwrite(strcat('beta_ROI_',num2str(c),...
%     'ID',young_name_40{ID}),...
%       [beta_160_ROI{id,c}],'precision',9)  