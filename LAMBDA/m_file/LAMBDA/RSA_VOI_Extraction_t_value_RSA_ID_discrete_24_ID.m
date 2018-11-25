%% Part 0-b-ii: RSA : VOI extraction -t values
% Part0-b: To generate extract t maps within the interested VOIs.
% analysis afterwards.
% ii: use t maps (24 unique trials with discrete distance)
% [Note that session information for each subject has already been collapsed into single t maps]
% instead of beta maps (144  unique trials with continuous distance)
addpath('D:\GoogleDrive\Lambda_code\m_file\LAMBDA\tool_code');% enable read_mixed_csv()

%% Constants
% Parameters
VOI_name={"R_V1","L_V1",...
    "R_IPS","L_IPS"}';
peak_coordinates={...
    [30 -90 -2], [-26,-94,-4],...
    [28 -68 40],[-26 -56 48]};
radius=8; % Set the radius of the shere
amount_RSA_id=24; % amount of unique RSA id

% Files and paths
file_common_mask='D:\Yun-Shiuan_LAMBDA\template\for_this_study\binary_overlapped_WM_GM_common_mask_n40_with_signal_with_reslice.nii,1';
file_valid_run='D:\Yun-Shiuan_LAMBDA\Run_inclusion_info\inclusive_runs_indexes.csv';
path_t_maps='D:\Yun-Shiuan_LAMBDA\RSA\Part0-a_GLM_estimate_t_value_RSA_ID_discrete';
path_VOI_output='D:\Yun-Shiuan_LAMBDA\RSA\Part0-b_VOIs_t_value_extracted_RSA_ID_discrete';
files_t_maps=strcat('spmT_',pad(string([1:24]'),4,'left','0'),'.nii');

%% Read in run inclusion index info (not that only subject info in used for t maps)
run_inclusion_index=cellfun(@(x) regexprep(x,'"',''),...
    read_mixed_csv(file_valid_run,','),'un',0);
run_inclusion_index=table(run_inclusion_index(2:end,2),run_inclusion_index(2:end,3),...
    'VariableNames',{'sub_id','run_num'});
%% Derive subjects with valid runs
subject_list=unique(run_inclusion_index.sub_id);

%% Extract VOIs from Beta series=====================================================================
%%Preprocessing

% Read in the common mask for voxel extraction restriction
V_mask=spm_vol(file_common_mask);% Need transfrom to mni space(because con_0001.nii is in mni space)
Y_mask=spm_read_vols(V_mask);
T_mask=V_mask.mat; %obtain the transformation rule to enable transform to mni space
c_mask=find(Y_mask);% Convert the 3d index to a linear index (which 1 refers to "within-mask")
[x y z]=ind2sub(size(Y_mask),c_mask);
Y_mni_mask=cor2mni([x y z],T_mask);%transform to mni space corrdinates

%% Loop over subjects
for i=1:size(subject_list,1)
    
    % The output direction for VOI output per subject run
    path_VOI_output_id=fullfile(path_VOI_output,subject_list{i});
    
    % Skip if already done on the subject run
    if(length(ls([path_VOI_output_id,'/t_map*']))==amount_RSA_id*size(VOI_name,1))
        strcat('Already done! sub_id :',subject_list{i})
    else
        mkdir(path_VOI_output_id);
        %% Loop over t maps
        for t_map=1:size(files_t_maps,1)
            %Read in the beta.nii
            
            file_t_map_image=fullfile(path_t_maps,subject_list{i},files_t_maps{t_map});
            V_t_map=spm_vol(file_t_map_image);
            
            %% Loop over VOIs
            for v=1:size(VOI_name,1)
                % Extract the VOI (overlapped with the common mask)
                center=peak_coordinates{v};
                Dist_sq = (bsxfun(@minus,center(1,1),Y_mni_mask(:,1)).^2+...
                    bsxfun(@minus,center(1,2),Y_mni_mask(:,2)).^2+...
                    bsxfun(@minus,center(1,3),Y_mni_mask(:,3)).^2);
                sphere_index=(find(Dist_sq<radius^2));%indexes for the mni coordinates in the sphere
                %% The Critical Overlaying Step
                ROI_overlaid=Y_mni_mask(sphere_index,:);  %Extract the mni coordinates in shpere and mask
                ROI_overlaid=round(ROI_overlaid);
                %% A vector which each element is a t value of a voxel within the ROI
                t_map_in_ROI=spm_summarise(V_t_map,ROI_overlaid');
                %% Output the beta vector
                betas_name=strcat('t_map_RSA_',num2str(t_map),'_',VOI_name{v});
                file_VOI_output=fullfile(path_VOI_output_id,char(strcat(betas_name,'.mat')));
                save(file_VOI_output,'t_map_in_ROI');
            end
        end
    end
end



