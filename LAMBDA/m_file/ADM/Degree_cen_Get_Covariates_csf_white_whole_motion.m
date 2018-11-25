%% Part 3: Degree Centrality: Right after Detrend and Bandpas, Get Brain-Relevant Covariates for partialing out in the next step
%Covariates for "Partial out to remain only the Residuals" . They include:
%lateral ventricle activations/ white matter/ whole brain/ head motions, and also their
%first temporal derivatives (one-image-backward difference).

%In addtion, all task-related covariates should be imported(e.g., prob,
%gain,...), if dealing with task EPI data.(I will deal with it in the next part)
path='/bml/Data/Bank6/ADM-YunShiuan';
cd(path);
addpath(strcat(path,'/Scripts(m)_files'));% to enable using the function "cor2mni"


id=csvread(strcat(path,'/young_idlist.csv')); %read only the young's IDs
young_name=strcat('S0',num2str(id));
young_name=cellstr(young_name); %convert char to cell to enable choosing string by index
young_name_40=cellstr(strvcat(young_name{~ismember(young_name,{'S031' 'S040' 'S044'}')})); %Exclude S031,S040,S044
%% Derive Deep cerebral White Matter/ Lateral Ventricle/ 6 motions of each person
%% Derive Deep cerebral White Matter: "collect_white_deep_signal"
% Get the mask of Deep Cerabal White Matter
V_white=spm_vol('/usr/local/spm12/tpm/TPM.nii,2');
Y_white=spm_read_vols(V_white);%voxel based coordinates(need to convert to MNI space later, because EPI has already been normalized to MNI space)
Y_white_deep=find(Y_white>0.95);
[x_white_deep y_white_deep z_white_deep]=ind2sub(size(Y_white),Y_white_deep);%voxel based coordinates of deep white matter

T_white=V_white.mat;%convert to mni space
Y_mni_white_deep=cor2mni([x_white_deep y_white_deep z_white_deep],T_white);
% Y_mni_white_deep=round(Y_mni_white_deep); %Round after ovelaying

% Display the deep white mask (save nii file for spm- display)
% Waring: white matter include cerebelum and brain stem......
% img=NaN([121,145,121]);
%     vol=V_white;
%     file_name=char(strcat(path,'/white_deep.nii'));
%     vol.fname=file_name;
%     vol.n=[1 1];
% 
%       Y_cor_white=mni2cor([Y_mni_white_deep(:,1),Y_mni_white_deep(:,2), Y_mni_white_deep(:,3)],T_white);
%       img(sub2ind([121,145,121], Y_cor_white(:,1) , Y_cor_white(:,2), Y_cor_white(:,3)))=1;
%       spm_write_vol(vol,img);    

%Overlay with Harvard-Oxford Atlas(cerebral White matter) to exclude cerebelum and brain stem
V_hav_ox=spm_vol(strcat(path,'/degree_centrality/Part_2_partial_out_to_remain_only_residuals/covarites/mask_for_extract_white_csf/HarvardOxford-sub-maxprob-thr50-1mm.nii,1'));
Y_hav_ox=spm_read_vols(V_hav_ox);%voxel based coordinates(need to convert to MNI space later, because EPI has already been normalized to MNI space)
Y_cerebral_white=find(Y_hav_ox==1|Y_hav_ox==12);
[x_cerebral_white y_cerebral_white z_cerebral_white]=ind2sub(size(Y_hav_ox),Y_cerebral_white);%voxel based coordinates of cerebral white matter

T_hav_ox=V_hav_ox.mat;%convert to mni space
Y_mni_cerebral_white=cor2mni([x_cerebral_white y_cerebral_white z_cerebral_white],T_hav_ox);
% Y_mni_cerebral_white=round(Y_mni_cerebral_white); %round after overlaying

%Overlaying: Downsize voxel amount from 18549 to 16377(excluding brain stem and cerebelum)
Y_mni_white_deep_cerebral=intersect(round(Y_mni_white_deep),round(Y_mni_cerebral_white),'rows');


% Display the overlaid deep cerebral white mater mask (save nii file for spm- display)
%     img=NaN([121,145,121]);
%     vol=V_white;    
%     file_name=char(strcat(path,'/degree_centrality/Part_2_partial_out_to_remain_only_residuals/covarites/white_cerebral_deep.nii'));
%     vol.fname=file_name;
%     vol.n=[1 1];
% 
%       Y_cor_white_cerebral=mni2cor([Y_mni_white_deep_cerebral(:,1),Y_mni_white_deep_cerebral(:,2), Y_mni_white_deep_cerebral(:,3)],T_white);
%       img(sub2ind([121,145,121], Y_cor_white_cerebral(:,1) ,Y_cor_white_cerebral(:,2), Y_cor_white_cerebral(:,3)))=1;
%       spm_write_vol(vol,img);    

collect_white_deep_signal={};
for id=1:length(young_name_40);
    for run=1:5;
% Get the names of 5 raw EPI images of the person
[imgfilename,dirs]=spm_select('ExtFPListRec',...%FORMAT [files,dirs] = spm_select('ExtFPListRec',direc,filt,frames)
           strcat(path,'/degree_centrality/Part_1_raw_data_and_detrend_and_bandpass/',young_name_40{id},...
           '/wradm',num2str(run),'.nii_after_bandpass'),...
           'Filtered_4DVolume.nii',1:218);%In MNI space (already normalize)
imgfilename=cellstr(imgfilename);       
% SPM_summarise 
signal_white_deep=spm_summarise(imgfilename,Y_mni_white_deep_cerebral',@mean);% a 1090x1 time serie of deep white matter signal
collect_white_deep_signal{id,run}=signal_white_deep;
[id run]
end
end

%% Derive CSF: "collect_csf_signal"
% Get the mask of CSF region
V_csf=spm_vol('/usr/local/spm12/tpm/TPM.nii,3');
Y_csf=spm_read_vols(V_csf);%voxel based coordinates(need to convert to MNI space later, because EPI has already been normalized to MNI space)
Y_csf_deep=find(Y_csf>0.9);
[x_csf_deep y_csf_deep z_csf_deep]=ind2sub(size(Y_csf),Y_csf_deep);%voxel based coordinates of deep white matter

T_csf=V_csf.mat;%convert to mni space
Y_mni_csf_deep=cor2mni([x_csf_deep y_csf_deep z_csf_deep],T_csf);
Y_mni_csf_deep=round(Y_mni_csf_deep);

% Display the deep white mask
%Waring: Falsely covering the eyes and other non-lateral ventricels
% img=NaN([121,145,121]);
%     vol=V_csf;
%     file_name=char(strcat(path,'/csf_deep.nii'));
%     vol.fname=file_name;    
% 
%       Y_cor_csf=mni2cor([Y_mni_csf_deep(:,1),Y_mni_csf_deep(:,2), Y_mni_csf_deep(:,3)],T_csf);
%       img(sub2ind([121,145,121], Y_cor_csf(:,1) , Y_cor_csf(:,2), Y_cor_csf(:,3)))=1;
%       spm_write_vol(vol,img);  

%Overlay with Harvard-Oxford Atlas(cerebral White matter) to exclude eyes and other non-lateral ventricels
Y_lateral_ventricle=find(Y_hav_ox==3|Y_hav_ox==14);
[x_lateral_ventricle y_lateral_ventricle z_lateral_ventricle]=ind2sub(size(Y_hav_ox),Y_lateral_ventricle);%voxel based coordinates of cerebral white matter

T_hav_ox=V_hav_ox.mat;%convert to mni space
Y_mni_lateral_ventricle=cor2mni([x_lateral_ventricle y_lateral_ventricle z_lateral_ventricle],T_hav_ox);
% Y_mni_lateral_ventricle=round(Y_mni_lateral_ventricle); %round after overlaying

%Overlaying: Downsize voxel amount from  16073  to  5711 (exclude eyes and other non-lateral ventricels)
Y_mni_csf_deep_ventricle=intersect(round(Y_mni_csf_deep),round(Y_mni_lateral_ventricle),'rows');


%Display the overlaid deep cerebral white mater mask (save nii file for spm- display)
%     img=NaN([121,145,121]);
%     vol=V_white;    
%     file_name=char(strcat(path,'/degree_centrality/Part_2_partial_out_to_remain_only_residuals/covarites/csf_ventricle_deep.nii'));
%     vol.fname=file_name;
%     vol.n=[1 1];
% 
%       Y_cor_csf_ventricle=mni2cor([Y_mni_csf_deep_ventricle(:,1),Y_mni_csf_deep_ventricle(:,2),Y_mni_csf_deep_ventricle(:,3)],T_white);
%       img(sub2ind([121,145,121], Y_cor_csf_ventricle(:,1) ,Y_cor_csf_ventricle(:,2), Y_cor_csf_ventricle(:,3)))=1;
%       spm_write_vol(vol,img);    

collect_csf_deep_signal={};
for id=1:length(young_name_40);    % Get the names of 5 raw EPI images of the person
    for run=1:5;
[imgfilename,dirs]=spm_select('ExtFPListRec',...%FORMAT [files,dirs] = spm_select('ExtFPListRec',direc,filt,frames)
           strcat(path,'/degree_centrality/Part_1_raw_data_and_detrend_and_bandpass/',young_name_40{id},...
           '/wradm',num2str(run),'.nii_after_bandpass'),...
           'Filtered_4DVolume.nii',1:218);
imgfilename=cellstr(imgfilename);       
% SPM_summarise 
signal_csf_deep=spm_summarise(imgfilename,Y_mni_csf_deep_ventricle',@mean);% a 1090x1 time serie of deep csf signal
collect_csf_deep_signal{id,run}=signal_csf_deep;
[id,run]
    end
end

%% Derive Whole Brain: "collect_whole_signal"

%Use Harvard-Oxford Atlas (excluding cerebellum and brain stem, since we only care about cortical and sub-cortical connectivity) to define whole brain
Y_whole_brain=find(Y_hav_ox~=8 & Y_hav_ox~=0);%excluding brain stem(cerebellum was not included in the first place) and ouside-brain background
[x_whole_brain y_whole_brain z_whole_brain]=ind2sub(size(Y_hav_ox),Y_whole_brain);%voxel based coordinates of cerebral white matter

T_hav_ox=V_hav_ox.mat;%convert to mni space
Y_mni_whole_brain=cor2mni([x_whole_brain y_whole_brain z_whole_brain],T_hav_ox);


%Display the whole brain mask (save nii file for spm- display)
%     img=NaN([182,218,182]); %NOTE: the size is corresponding to
%     hav_ox.nii. If TPM.nii is overlaid with hav_ox.nii, then one should
%     use the size of TPM.nii([121,145,121]), since it has lower spatial resolution.
%     vol=V_hav_ox;    
%     file_name=char(strcat(path,'/degree_centrality/Part_2_partial_out_to_remain_only_residuals/covarites/whole_brain.nii'));
%     vol.fname=file_name;
%     vol.n=[1 1];
% 
%       Y_cor_whole_brain=mni2cor([Y_mni_whole_brain(:,1),Y_mni_whole_brain(:,2),Y_mni_whole_brain(:,3)],T_hav_ox);
%       img(sub2ind([182,218,182], Y_cor_whole_brain(:,1) ,Y_cor_whole_brain(:,2), Y_cor_whole_brain(:,3)))=1;
%       spm_write_vol(vol,img);    

collect_whole_signal={};
for id=1:length(young_name_40);    % Get the names of 5 raw EPI images of the person
    for run=1:5;
[imgfilename,dirs]=spm_select('ExtFPListRec',...%FORMAT [files,dirs] = spm_select('ExtFPListRec',direc,filt,frames)
           strcat(path,'/degree_centrality/Part_1_raw_data_and_detrend_and_bandpass/',young_name_40{id},...
           '/wradm',num2str(run),'.nii_after_bandpass'),...
           'Filtered_4DVolume.nii',1:218);
imgfilename=cellstr(imgfilename);       
% SPM_summarise 
signal_whole=spm_summarise(imgfilename,Y_mni_whole_brain',@mean);% a 1090x1 time serie of deep csf signal
collect_whole_signal{id,run}=signal_whole;
[id,run]
    end
end

%% Save as .mat
save(strcat(path,'/degree_centrality/Part_2_partial_out_to_remain_only_residuals/',...
    'covarites/derived_covariate_mat/collect_white_csf_whole_signal.mat'),...
    'collect_white_deep_signal','collect_csf_deep_signal','collect_whole_signal');         

%Check how much does excluding eyes/ brainstem/ cerebelum influence the
%overall averaged signal.---> As for CSF, it matters A LOT!!!
% %cd('/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_2_partial_out_to_remain_only_residuals/covarites')
% % load('collect_white_csf_signal_overlaid.mat')
% % csf_2=collect_csf_deep_signal;
% % white_2=collect_white_deep_signal;
% % load('collect_white_csf_signal_unoverlaid_caution.mat')
% % csf_1=collect_csf_deep_signal;
% % white_1=collect_white_deep_signal;

%% Derive Head Motion Parameters
% Get 6 Head Motion Parameters data

%Run it only if copying is needed
% % Copy the parameter txt file into degree centrality folder
% % for id=1:length(young_name_40);
% % path_old = strcat(path,'/IMG_nii/',young_name_40{id},'/degree_cen_normalise_3x3x3/rp_adm*');
% % path_new = strcat(path,'/degree_centrality/Part_2_partial_out_to_remain_only_residuals/covarites/head_motion_raw/',young_name_40{id},'/');
% % copyfile (path_old,path_new); 
% % end

% Read them in (5 files, each per run. Each file contains 6 head motion parameters.)
%Save them to a big cell
collect_6_motion={};
for id=1:length(young_name_40);
path_new = strcat(path,'/degree_centrality/Part_2_partial_out_to_remain_only_residuals/covarites/head_motion_raw/',young_name_40{id},'/');
D=dir(path_new);% get the 5 file names
D=struct2cell(D);
filename=D(1,:);
file_index=~cell2mat(D(4,:));%get the indexes of files, exluding folders
filename_valid=filename(file_index);

for run=1:length(filename_valid);%(5 files, each per run)
file_locate=strcat(path_new,filename_valid{run});
fileID = fopen(file_locate);
collect_6_motion{id,run} = textscan(fileID,'%f %f %f %f %f %f'); %'%f' for double precision
fclose(fileID);
end
end
%Check if each run of each person has 218 values
% a=cellfun(@(x)strcat(num2str(length(x{1})),num2str(length(x{2})),...
%     num2str(length(x{3})),num2str(length(x{4})),num2str(length(x{5}))),collect_6_motion,'un',0)
save(strcat(path,'/degree_centrality/Part_2_partial_out_to_remain_only_residuals/',...
    'covarites/derived_covariate_mat/collect_6_head_motion.mat'),'collect_6_motion')  

%% Derive first Derivitive terms of White Matter, CSF, and 6 Head Motions
cd('/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_2_partial_out_to_remain_only_residuals/covarites/derived_covariate_mat/')
% load('collect_white_csf_signal.mat');
% load('collect_6_head_motion.mat');
collect_white_deep_signal_1_drv=cellfun(@(x) [0;x(2:end)]-[0;x(1:end-1)],collect_white_deep_signal,'un',0);
collect_csf_deep_signal_1_drv=cellfun(@(x) [0;x(2:end)]-[0;x(1:end-1)],collect_csf_deep_signal,'un',0);
collect_whole_signal_1_drv=cellfun(@(x) [0;x(2:end)]-[0;x(1:end-1)],collect_whole_signal,'un',0);
collect_6_motion_1_drv=cellfun(@(x)...
    arrayfun(@(p)[0;x{p}(2:end)]-[0;x{p}(1:end-1)],[1,2,3,4,5,6],'un',0),collect_6_motion,'un',0);
%Check accuracy
%[collect_csf_deep_signal{1}(2:6) collect_csf_deep_signal{1}(1:5) collect_csf_deep_signal_1_drv{1}(2:6)]
%Check accuracy
% [collect_6_motion{1}{1}(2:6) collect_6_motion{1}{1}(1:5) collect_6_motion_1_drv{1}{1}(2:6)]
save(strcat(path,'/degree_centrality/Part_2_partial_out_to_remain_only_residuals/',...
    'covarites/derived_covariate_mat/collect_6_head_motion_and_white_csf_whole_signal_1_drv.mat'),...
'collect_white_deep_signal_1_drv','collect_csf_deep_signal_1_drv','collect_whole_signal_1_drv','collect_6_motion_1_drv')  

