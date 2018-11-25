%% Part 0-b-iii: RSA : VOI extraction -t values
% Part0-b: To generate extract t maps within the interested VOIs.
% analysis afterwards.

% iv: RDA_ID_discrete 18 ID:
% - RSA ID : 1~18 (3 formats x 3 distance levels x 2 directions)
% - The motivation is to increase the statistical efficiency (Conditions in the Trial 4 is too noisy),
% while still keep that same-notation trials have the same statistical efficiency
% as the cross-notation trials (as Trial 4 did).
% [LL: 12 trials = 3 distance levels x 2 direction x 2 repetition]
% [FF: 12 trials = 3 distance levels x 2 direction x 2 repetition]
% [cross-notation: 12 trials = 3 distance levels x 2 direction x 2 [repetition]
% (Having more repetition in a condition entails a smaller degree of random noise.)

% The motivation for having 18 trials is to ensure beta estimate of each condition is subjected to same amount of noise.
% (Having more repetition in a condition entails a smaller degree of random noise.)

% [Note that session information for each subject has already been collapsed into single t maps]
% instead of beta maps (144  unique trials with continuous distance)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%IPS adjusted by emperical results (2018/06/10)
%Update subject list for children (i.e., inclusive_runs_indexes_new_June_10.csv) (2018/6/10)
%Update subject list for children (i.e., inclusive_runs_indexes_new_June_11.csv) (2018/6/10)

PATH_ROOT='D:\Yun-Shiuan_LAMBDA';
PATH_TOOL_CODE='D:\GoogleDrive\Lambda_code\m_file\LAMBDA\tool_code';
addpath(genpath(fullfile(PATH_ROOT,'rsatoolbox')));
addpath(PATH_TOOL_CODE);% enable read_mixed_csv()
%% Constants
% Parameters
VOI_NAME={"R_V1","L_V1",...
          "R_IPS","L_IPS",'R_DLPFC'}';
% New meta-analytic peaks (but IPS adjusted by empirical results)
PEAK_COORDINATES={...
    [14,-84,6], [-9,-89,-6],...
    [28,-64,40],[-26,-56,48],[44,28,22]};
% % New meta-analytic peaks
% PEAK_COORDINATES={...
%     [14,-84,6], [-9,-89,-6],...
%     [26,-64,36],[-26,-64,36],[44,28,22]};
% Old empirical peaks
% PEAK_COORDINATES={...
%     [30,-90,-2], [-26,-94,-4],...
%     [28,-68,40],[-26,-56,48],[48,34,22]};
RADIUS=8; % Set the radius of the shere
AMOUNT_RSA_ID=18; % amount of unique RSA id
VERSION_OUTPUT='trial_22';
% Files and paths
% Files and paths
% Paths
PATH_T_MAPS_CHILD=fullfile(PATH_ROOT,'RSA','trial_21','Part0-a_GLM_estimate_t_value','child');
PATH_VOI_OUTPUT_CHILD=fullfile(PATH_ROOT,'RSA',VERSION_OUTPUT,'Part0-b_VOIs_t_value_extracted','child');
PATH_T_MAPS_ADULT=fullfile(PATH_ROOT,'RSA','trial_7_id_18_include_adult','Part0-a_GLM_estimate_t_value','Adult');
PATH_VOI_OUTPUT_ADULT=fullfile(PATH_ROOT,'RSA',VERSION_OUTPUT,'Part0-b_VOIs_t_value_extracted','adult');
% Files
FILE_COMMON_MASK=fullfile(PATH_ROOT,'template','for_this_study','binary_overlapped_WM_GM_common_mask_n40_with_signal_with_reslice.nii,1');
FILE_VALID_RUN_CHILD=fullfile(PATH_ROOT,'Run_inclusion_info','inclusive_runs_indexes_new_June_11.csv');
FILE_VALID_RUN_ADULT=fullfile(PATH_ROOT,'Adult','Run_inclusion_info','inclusive_runs_indexes.csv');

FILES_T_MAPS=strcat('spmT_',pad(string([1:AMOUNT_RSA_ID]'),4,'left','0'),'.nii');

%% Read in run inclusion index info (not that only subject info in used for t maps)
run_inclusion_index_child=read_mixed_csv_to_table(FILE_VALID_RUN_CHILD);
run_inclusion_index_adult=read_mixed_csv_to_table(FILE_VALID_RUN_ADULT);

%% Derive subjects with valid runs
subject_list_child=unique(run_inclusion_index_child.sub_id);
subject_list_adult=unique(run_inclusion_index_adult.sub_id);
% subject_list=subject_list(~ismember(subject_list,{'XFC305'})); %305 have
% missing run (already solved)
mkdir(PATH_VOI_OUTPUT_CHILD)
mkdir(PATH_VOI_OUTPUT_ADULT)

%% Extract VOIs from Beta series=====================================================================
%%Preprocessing

% Read in the common mask for voxel extraction restriction
V_mask=spm_vol(FILE_COMMON_MASK);% Need transfrom to mni space(because con_0001.nii is in mni space)
Y_mask=spm_read_vols(V_mask);
T_mask=V_mask.mat; %obtain the transformation rule to enable transform to mni space
c_mask=find(Y_mask);% Convert the 3d index to a linear index (which 1 refers to "within-mask")
[x, y, z]=ind2sub(size(Y_mask),c_mask);
Y_mni_mask=cor2mni([x y z],T_mask);%transform to mni space corrdinates

%% Loop over subjects
for age_group_index=1:2 %1 for child, 2for adult
    if age_group_index==1
        path_VOI_output=PATH_VOI_OUTPUT_CHILD;
        subject_list=subject_list_child;
        path_t_maps=PATH_T_MAPS_CHILD;
    else
        path_VOI_output=PATH_VOI_OUTPUT_ADULT;
        subject_list=subject_list_adult;
        path_t_maps=PATH_T_MAPS_ADULT;
    end
    for sub_id=1:numel(subject_list)    
        % The output direction for VOI output per subject run
        path_VOI_output_id=fullfile(path_VOI_output,subject_list{sub_id});
        mkdir(path_VOI_output_id)
        % Skip if already done on the subject run
        if(length(dir([path_VOI_output_id,'/t_map*']))==AMOUNT_RSA_ID*size(VOI_NAME,1))
            strcat('Already done! sub_id :',subject_list{sub_id})
        else
            %% Loop over t maps
            for t_map=1:numel(FILES_T_MAPS)
                %Read in the Tmap.nii

                file_t_map_image=fullfile(path_t_maps,subject_list{sub_id},FILES_T_MAPS{t_map});
                V_t_map=spm_vol(file_t_map_image);

                %% Loop over VOIs
                for v=1:numel(VOI_NAME)
                    % Extract the VOI (overlapped with the common mask)
                    center=PEAK_COORDINATES{v};
                    Dist_sq = (bsxfun(@minus,center(1,1),Y_mni_mask(:,1)).^2+...
                        bsxfun(@minus,center(1,2),Y_mni_mask(:,2)).^2+...
                        bsxfun(@minus,center(1,3),Y_mni_mask(:,3)).^2);
                    sphere_index=(find(Dist_sq<RADIUS^2));%indexes for the mni coordinates in the sphere
                    %% The Critical Overlaying Step
                    ROI_overlaid=Y_mni_mask(sphere_index,:);  %Extract the mni coordinates in shpere and mask
                    ROI_overlaid=round(ROI_overlaid);
                    %% A vector which each element is a t value of a voxel within the ROI
                    t_map_in_ROI=spm_summarise(V_t_map,ROI_overlaid');
                    %% Output the beta vector
                    betas_name=strcat('t_map_RSA_',num2str(t_map),'_',VOI_NAME{v});
                    file_VOI_output=fullfile(path_VOI_output_id,char(strcat(betas_name,'.mat')));
                    save(file_VOI_output,'t_map_in_ROI');
                end
            end
        end
    end
end



