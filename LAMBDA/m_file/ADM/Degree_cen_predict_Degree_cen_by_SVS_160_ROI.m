%% Part 7-A: Degree Centrality: Regress Degree Centrality by SVS(D~SVS+sex) (160 ROI)
%% Two Mode to choose: Exclude Crebelum or NOT
%% Other Two Mode to Choose: Standardized or Unstandardized
%% T map and P map are generated in R enviroment"/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_4_Regression_degree_cen_by_SVS_by_R/R_script_and_RData/regression_deree_cen_by_SVS.R"
%addpath('/usr/local/spm12/')
%spm_jobman('initcfg') -- before running the batch
path='/bml/Data/Bank6/ADM-YunShiuan';
cd(path);
addpath(strcat(path,'/Scripts(m)_files'));% to enable using the function "cor2mni"

id=csvread(strcat(path,'/young_idlist.csv')); %read only the young's IDs
young_name=strcat('S0',num2str(id));
young_name=cellstr(young_name); %convert char to cell to enable choosing string by index
young_name_40=cellstr(strvcat(young_name{~ismember(young_name,{'S031' 'S040' 'S044'}')})); %Exclude S031,S040,S044

%% Load in data-->Do it in the R enviroment to utilize ggplot and "label4MRI"
%SVS and sex
% load('id_age_gender_sec_hed_sti.mat');
% id=a(:,1);
% id_valid_index=~ismember(id,{'S031' 'S040' 'S044'}');
% age=str2double(a(id_valid_index,2));
% gender=str2double(a(id_valid_index,3));
% hedonism=str2double(a(id_valid_index,5));
% security=str2double(a(id_valid_index,4));
% %D matrix 
% cd('/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_3_R_and_Z_Matrix_of_160_ROI/R_D_and_Z_matrix');
% load('collect_D_matrix_ROI_160_id_40.mat');
%% Regression: Degree Centrality~Security+Hedonism+sex (each per ROI)-->Done in the R enviroment to utilize ggplot
%% Save T map and P map as nii. file (load the .mat generated in R)
%Read in the lm result file
cd('/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_4_Regression_degree_cen_by_SVS_by_R/lm_result_unstandardized/without_cerebelum')
load('collect_lm_160_ROI_by_2_SVS.mat');

%Read in 160 ROI location
load('/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Dosenbach_Science_160ROIs_Center.mat')
% Get all coordinate 
ROI_160_center=Dosenbach_Science_160ROIs_Center;
ROI_160_center=[[1:size(ROI_160_center,1)]' ROI_160_center]; %Add ID for each center(1~160)

%% Exclude cerebelum (use it only if needed)
load('/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_4_Regression_degree_cen_by_SVS_by_R/lm_result_unstandardized/without_cerebelum/cerebelum_aal_name.mat');
cerebelum_index=ismember(aal_name,cerebelum_aal_name);
ROI_160_center=ROI_160_center(~cerebelum_index,:); %(142 ROI remained)


%Get the grey_white mask
%% Create grey_white mask which will be overlaid with ROI to create a combined ROI for extracting Beta later
 V_mask=spm_vol('/bml/Data/Bank6/ADM-YunShiuan/smooth_mask_grey_white.nii');% Need transfrom to mni space(because con_0001.nii is in mni space)
 Y_mask=spm_read_vols(V_mask);
 Y_mask=Y_mask>0.8;%Set criteria >0.8 as brain(white/grey) region
 T_mask=V_mask.mat; %abtain the transformation rule to enable transform to mni space
 c_mask=find(Y_mask==1);% select the voxel which mask value==1
 [x y z]=ind2sub(size(Y_mask),c_mask);
 Y_mni_mask=cor2mni([x y z],T_mask);%transform to mni space corrdinates
 %% T map
 t_map={hedonism_t,security_t};
 map_name={'hedonism_t_positive','hedonism_t_negative'...
           'security_t_positive','security_t_negative'};
   
 order=1;
for SVS=1:2; % two SVS
    for sign=1:2; %one SVS, two signs
ROI_overlaid_all=[];
img=NaN([121,145,121]);

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
Y_cor_roi=mni2cor([round(ROI_overlaid(:,1)),round(ROI_overlaid(:,2)), round(ROI_overlaid(:,3))],T_mask);
if sign==1 % positive sign
img(sub2ind([121,145,121], Y_cor_roi(:,1) , Y_cor_roi(:,2), Y_cor_roi(:,3)))= t_map{SVS}(c);
else %negative map
img(sub2ind([121,145,121], Y_cor_roi(:,1) , Y_cor_roi(:,2), Y_cor_roi(:,3)))= -1*t_map{SVS}(c);
end 

[SVS sign c]
end
assignin('base', map_name{order},img);

file_name=char(strcat('/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_4_Regression_degree_cen_by_SVS_by_R/lm_result_unstandardized/without_cerebelum/',...
    'lm_160_ROI_',map_name{order},'.nii'));
vol=V_mask;%Refresh the Vol info
vol.fname=file_name;
%convert Y_mni_mask to matrix form(voxle-based coordinates)
spm_write_vol(vol,img);
order=order+1;
    end
end
%% Overall Degree Centrality T Map
cd('/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_3_R_and_Z_Matrix_of_160_ROI/2_R_D_and_Z_matrix/exclude_cerebelum');
load('collect_Z_matrix_ROI_160_id_40.mat');%Load in 'collect_Z_matrix' (Standardized Appraoch)
%load('collect_D_matrix_ROI_160_id_40.mat');%Load in 'collect_D_matrix' (Unstandardized Approach)
% overall_t_map=mean(collect_Z_matrix,1)./(std(collect_Z_matrix,1));
[h,p,ci,stats]=ttest(collect_Z_matrix);
%[h,p,ci,stats]=ttest(collect_D_matrix); (Unstandardized Approach)
overall_t_map=stats.tstat;
ROI_overlaid_all=[];
img=NaN([121,145,121]);
cd('/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_4_Regression_degree_cen_by_SVS_by_R/lm_result_unstandardized/without_cerebelum/')

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
Y_cor_roi=mni2cor([round(ROI_overlaid(:,1)),round(ROI_overlaid(:,2)), round(ROI_overlaid(:,3))],T_mask);
img(sub2ind([121,145,121], Y_cor_roi(:,1) , Y_cor_roi(:,2), Y_cor_roi(:,3)))= overall_t_map(c);
c
end
file_name=char(strcat('/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_4_Regression_degree_cen_by_SVS_by_R/lm_result_unstandardized/without_cerebelum/',...
    'overall_160_ROI_degree_cen_t.nii'));
vol=V_mask;%Refresh the Vol info
vol.fname=file_name;
%convert Y_mni_mask to matrix form(voxle-based coordinates)
spm_write_vol(vol,img);

%% p map
 p_map={hedonism_p,security_p};
 map_name={'hedonism_p','security_p'};

 order=1;
for SVS=1:2; % two SVS
ROI_overlaid_all=[];
img=NaN([121,145,121]);

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
Y_cor_roi=mni2cor([round(ROI_overlaid(:,1)),round(ROI_overlaid(:,2)), round(ROI_overlaid(:,3))],T_mask);
img(sub2ind([121,145,121], Y_cor_roi(:,1) , Y_cor_roi(:,2), Y_cor_roi(:,3)))= p_map{SVS}(c);


[SVS c]
end
assignin('base', map_name{order},img);

file_name=char(strcat('/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_4_Regression_degree_cen_by_SVS_by_R/lm_result_unstandardized/without_cerebelum/',...
    'lm_160_ROI_',map_name{order},'.nii'));
vol=V_mask;%Refresh the Vol info
vol.fname=file_name;
%convert Y_mni_mask to matrix form(voxle-based coordinates)
spm_write_vol(vol,img);
order=order+1;
end
    

%% Overall Degree Centrality P Map
cd('/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_3_R_and_Z_Matrix_of_160_ROI/2_R_D_and_Z_matrix/exclude_cerebelum/');
load('collect_Z_matrix_ROI_160_id_40.mat');%Load in 'collect_Z_matrix'
%load('collect_D_matrix_ROI_160_id_40.mat'); Unstandardized Approach
% overall_t_map=mean(collect_Z_matrix,1)./(std(collect_Z_matrix,1));

[h,p,ci,stats]=ttest(collect_Z_matrix);
%[h,p,ci,stats]=ttest(collect_D_matrix); Unstandardized Appraoch

overall_p_map=p;
ROI_overlaid_all=[];
img=NaN([121,145,121]);
cd('/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_4_Regression_degree_cen_by_SVS_by_R/lm_result_unstandardized/without_cerebelum/')

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
Y_cor_roi=mni2cor([round(ROI_overlaid(:,1)),round(ROI_overlaid(:,2)), round(ROI_overlaid(:,3))],T_mask);
img(sub2ind([121,145,121], Y_cor_roi(:,1) , Y_cor_roi(:,2), Y_cor_roi(:,3)))= overall_p_map(c);
c
end
file_name=char(strcat('/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_4_Regression_degree_cen_by_SVS_by_R/lm_result_unstandardized/without_cerebelum/',...
    'overall_160_ROI_degree_cen_p.nii'));
vol=V_mask;%Refresh the Vol info
vol.fname=file_name;
%convert Y_mni_mask to matrix form(voxle-based coordinates)
spm_write_vol(vol,img);
save('overall_t_p_160_ROI.mat','overall_t_map','overall_p_map')
