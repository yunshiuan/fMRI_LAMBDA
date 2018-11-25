%% (Part2)Covert the mixed-effects results by R to nii. format which enable SPM to read
%addpath('/usr/local/spm12/')
%spm_jobman('initcfg') -- before running the batch
path='/bml/Data/Bank6/ADM-YunShiuan';
cd (strcat(path,'/Second_mixed_model_by_R/Betas_derived/new_df_t_beta'));
%% Read 13 csv. file (each contains a volumn of Beta coefficient by mixed model)
Beta_file=dir('*.csv'); % 13 Beta coefficient +1 valid voxel indexes
Beta_file=cellstr(strvcat(Beta_file.name))';%Convert to character array
Beta_file=setdiff(Beta_file,'valid_voxel_index.csv');%Exclude 'valid_voxel_index' from being a type of Beta
Beta_type=cellfun(@(x) regexprep(x,'_beta_t_df.csv',''),Beta_file,'un',0);

Beta_collect={};%A cell to collect all 13 Betas
t_collect={};%A cell to collect all 13 t
df_collect={};%A cell to collect all 13 df
for i=1:length(Beta_file);
df=csvread(Beta_file{i},1,1);%Each should contain 63030 rows
Beta_collect{1,i}=df(:,1);
t_collect{1,i}=df(:,2);
df_collect{1,i}=df(:,3);
end

cellfun('length',Beta_collect)%Check if each file has 63030 rows
arrayfun(@(x) all(df_collect{x} == df_collect{x}(1)),[1:13])%Check if all dfs in the same file are the same

%% Append xovel indexes to each Beta type
voxel_index=csvread('valid_voxel_index.csv',1,1);
voxel_index_linear=arrayfun(@(n) sub2ind([53,63,52],voxel_index(n,1),...
            voxel_index(n,2),voxel_index(n,3)),1:size(voxel_index,1));% Has already been ckecked to match with "c" in the script "Mixed_Extract_Betas"
cd(strcat(path,'/Second_mixed_model_by_R/nii_derived/new_df_t_beta/'))
 for i=1:length(Beta_type);
%% Fill in Beta values into each matrix
img=NaN([53,63,52]);%Creat a matrix with NaN which is going to be filled with Beta values
img(voxel_index_linear)=t_collect{i};
%% Write nii. files
%(Unused: make_nii/save.nii/load.nii)
% nii=make_nii(img,[3 3 3], [27 38.3 24.3],[16]);
%load nii relevant info which will be added to 'nii'
% nii_info=load_nii(strcat(path,'/First_design_matrix_with_contrastncovariate_categorical/S023/normalise_3x3x3/con_0001.nii'));
% nii.hdr=nii_info.hdr;
% save_nii(nii,strcat('../nii_derived/',Beta_type{i},'.nii'));
%(Used:sp_/write_vol/spm_vol)
%Read in vol information for t map (should be the same as the input of LMER)
vol=spm_vol(char(strcat(path,'/First_design_matrix_with_contrastncovariate_categorical/S023/normalise_3x3x3/spmT_0001.nii')));
%Modify the vol information
file_name=char(strcat(path,'/Second_mixed_model_by_R/nii_derived/new_df_t_beta/',...
               Beta_type{i},'.nii'));
vol.fname=file_name;
descrip_name=char(strcat('SPM{T_[',num2str(df_collect{i}(1)),'.0]} - contrast : ', Beta_type{i}));
vol.descrip=descrip_name;
spm_write_vol(vol,img);%need to change the fname (the location where file is going to be written)
end;

%(Used) Test the spm_write_vol(Pass)
% V=spm_vol('/bml/Data/Bank6/ADM-YunShiuan/Second_ROI_defined_by_value2prob_new_server/normalise_3x3x3/two_SVS_in_same_GLM/prob/spmT_0001.nii');
% Y=spm_read_vols(V);
% V.fname='/bml/Data/Bank6/ADM-YunShiuan/Second_ROI_defined_by_value2prob_new_server/normalise_3x3x3/two_SVS_in_same_GLM/prob/test.nii'
% spm_write_vol(V,Y)


% %%(Unused) Test the make_nii function--> Error: The x-axis  is reversed orz (origin is in the down-right corner)
% V=spm_vol('beta_0001.nii');%Specify the nii. files interested
% Y=spm_read_vols(V);
% nii_test=make_nii(Y,[3 3 3], [27 38.3 24.3],[16]);
% nii=load_nii('beta_0001.nii');%Get nii. relevant info(hdr)
% nii_test.hdr=nii.hdr;
% save_nii(nii_test,'Beta_test.nii');
% V_test=spm_vol('Beta_test.nii');
% Y_test=spm_read_vols(V_test);
% %% Although the a-aixe is reversed, it generate the same value by spm_vol
% 
