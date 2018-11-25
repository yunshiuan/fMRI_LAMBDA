%% Part 7-B(With R): Degree Centrality: Regress Degree Centrality by SVS(D~SVS+sex) (whole cerebral grey matter)
%% T map and P map are generated in R enviroment"/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_4_Regression_degree_cen_by_SVS_by_R/R_script_and_RData/regression_deree_cen_by_SVS.R"
%% Note: Since this is a whole-brian voxels appraoch, one could actually conduct second-level by SPM (See Part7-C).
%% By SPM, tidy result table could be output (peak location, k size, p_FDR etc.).
%% Alternatively, one could use XJview with those nii. generated in this script to output result table.
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
cd('/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_4_Regression_degree_cen_by_SVS_by_R/whole_brain_approach/lm_result_standardized')
load('collect_lm_whole_brain_by_2_SVS.mat');

%% Get the mni coordinates of the residual map within cerebral_grey_matter mask
%% Create grey mask from Hvard-Oxford Atlas which will be overlaid with the whole brain of S023
 V_mask=spm_vol(strcat(path,'/degree_centrality/Part_3_R_and_Z_Matrix_of_160_ROI/whole_brain_approach/1_raw_residual_whole_brain/cerebral_grey_mask.nii'));% Need transfrom to mni space(because con_0001.nii is in mni space)
 Y_mask=spm_read_vols(V_mask);
 T_mask=V_mask.mat; %obtain the transformation rule to enable transform to mni space
 c_mask=find(Y_mask==1);% 
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
Y_template_grey_index=ismember(Y_mni_template,Y_mni_mask,'rows');
Y_mni_template_grey=Y_mni_template(Y_template_grey_index,:);
Y_cor_template=[x y z];
Y_cor_template_grey=Y_cor_template(Y_template_grey_index,:);

cd('/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_4_Regression_degree_cen_by_SVS_by_R/whole_brain_approach/lm_result_standardized')
%% T map
 t_map={hedonism_t,security_t};
 map_name={'hedonism_t_positive','hedonism_t_negative'...
           'security_t_positive','security_t_negative'};   
 order=1;
for SVS=1:2; % two SVS
    for sign=1:2; %one SVS, two signs
img=NaN([53,63,52]);
if sign==1 % positive sign
img(sub2ind([53,63,52], Y_cor_template_grey(:,1) , Y_cor_template_grey(:,2), Y_cor_template_grey(:,3)))= t_map{SVS};
else %negative map
img(sub2ind([53,63,52], Y_cor_template_grey(:,1) , Y_cor_template_grey(:,2), Y_cor_template_grey(:,3)))= -1*t_map{SVS};

[SVS sign]
end
assignin('base', map_name{order},img);

file_name=char(strcat('/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_4_Regression_degree_cen_by_SVS_by_R/whole_brain_approach/lm_result_standardized/',...
    'lm_whole_brain_',map_name{order},'.nii'));
vol=V_template;%Refresh the Vol info
vol.fname=file_name;
%convert Y_mni_mask to matrix form(voxle-based coordinates)
spm_write_vol(vol,img);
order=order+1;
    end
end
%% Overall Degree Centrality T Map
cd('/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_3_R_and_Z_Matrix_of_160_ROI/whole_brain_approach/2_R_D_and_Z_matrix');
load('smoothed_collect_Z_matrix_whole_brain_id_40.mat');%Load in 'collect_Z_matrix' (Standardized Appraoch)
%load('collect_D_matrix_ROI_160_id_40.mat');%Load in 'collect_D_matrix' (Unstandardized Approach)
% overall_t_map=mean(collect_Z_matrix,1)./(std(collect_Z_matrix,1));
[h,p,ci,stats]=ttest(collect_Z_matrix_smoothed);
%[h,p,ci,stats]=ttest(collect_D_matrix); (Unstandardized Approach)
overall_t_map=stats.tstat;
img=NaN([53,63,52]);
img(sub2ind([53,63,52], Y_cor_template_grey(:,1) , Y_cor_template_grey(:,2), Y_cor_template_grey(:,3)))= overall_t_map;
file_name=char(strcat('/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_4_Regression_degree_cen_by_SVS_by_R/whole_brain_approach/lm_result_standardized/',...
    'overall_whole_brain_degree_cen_t.nii'));
vol=V_template;%Refresh the Vol info
vol.fname=file_name;
%convert Y_mni_mask to matrix form(voxle-based coordinates)
spm_write_vol(vol,img);

%% p map
 p_map={hedonism_p,security_p};
 map_name={'hedonism_p','security_p'};
 order=1;
 
for SVS=1:2; % two SVS
img=NaN([53,63,52]);
img(sub2ind([53,63,52], Y_cor_template_grey(:,1) , Y_cor_template_grey(:,2), Y_cor_template_grey(:,3)))= p_map{SVS};

[SVS sign]
assignin('base', map_name{order},img);

file_name=char(strcat('/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_4_Regression_degree_cen_by_SVS_by_R/whole_brain_approach/lm_result_standardized/',...
    'lm_whole_brain_',map_name{order},'.nii'));
vol=V_template;%Refresh the Vol info
vol.fname=file_name;
%convert Y_mni_mask to matrix form(voxle-based coordinates)
spm_write_vol(vol,img);
order=order+1;
    
end   
%% Overall Degree Centrality P Map
cd('/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_3_R_and_Z_Matrix_of_160_ROI/whole_brain_approach/2_R_D_and_Z_matrix');
load('smoothed_collect_Z_matrix_whole_brain_id_40.mat');%Load in 'collect_Z_matrix' (Standardized Appraoch)
%load('collect_D_matrix_ROI_160_id_40.mat');%Load in 'collect_D_matrix' (Unstandardized Approach)
% overall_t_map=mean(collect_Z_matrix,1)./(std(collect_Z_matrix,1));
[h,p,ci,stats]=ttest(collect_Z_matrix_smoothed);
%[h,p,ci,stats]=ttest(collect_D_matrix); (Unstandardized Approach)
overall_p_map=p;
img=NaN([53,63,52]);
img(sub2ind([53,63,52], Y_cor_template_grey(:,1) , Y_cor_template_grey(:,2), Y_cor_template_grey(:,3)))= overall_p_map;
file_name=char(strcat('/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_4_Regression_degree_cen_by_SVS_by_R/whole_brain_approach/lm_result_standardized/',...
    'overall_whole_brain_degree_cen_p.nii'));
vol=V_template;%Refresh the Vol info
vol.fname=file_name;
%convert Y_mni_mask to matrix form(voxle-based coordinates)
spm_write_vol(vol,img);

%% Use XJ view
addpath('/usr/local/spm12/toolbox/xjview9/')
%-->QAQ it doesn't know that the input images is T map

