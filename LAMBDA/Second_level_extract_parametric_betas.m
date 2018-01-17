%% Part 2-b - Second Level (Parametric)- Extract parametric first level betas from second level peaks.
% For further visualization in R.
% Extract First level Con_arameters(the parameter estimate, not the t value per se),
% (so that it's the raw parameter estimates, which can be used for simple slope analysis etc.)
% including Cond_Dist_All,
% from Second level results(Age-effect peak)-The corresponding GLM to the
% dependent varialbe,
% e.g., Extract B_Dist_ALL from Gist_all_by_age_peak
% e.g., Extract Dist_LL from Dist_LL_by_age_peak

addpath(char("D:\GoogleDrive\Lambda_code\m_file\LAMBDA\tool_code"));% enable read_mixed_csv()
path='D:\Yun-Shiuan_LAMBDA\First_level_parametric_design_matrix_and_result'; % To get the subject list with valid runs (see below)
cd(path);
subject_list=cellstr(ls('*')); % subject list (which has valid first level contrast, i.e., not being excluded)
subject_list=subject_list(~cellfun(@isempty,regexp(subject_list,'(?<=df)\d+','match')));% include only subject name
id=regexp(subject_list,'(?<=df10)\d+$','match');% convert to pure number for writing csv
id=str2double(cellstr([id{:}]'));
Cond_Parameters_folder_name={'B0_LL' 'B0_LF' 'B0_FF'... 
                   'Dist_All'...
                   'Dist_LL' 'Dist_LF' 'Dist_FF'... 
                   't_LF-LL' 't_LL-LF' ...% Note that the values to be extracted is NOT
                   't_FF-LL' 't_LL-FF' ...% t value, but the contrast value (e.g.,con_001)
                   't_FF-LF' 't_LF-FF'};  % The names here is just to match with the folder name          
Cond_Parameters_real_name={'B0_LL' 'B0_LF' 'B0_FF'... 
                   'Dist_All'...
                   'Dist_LL' 'Dist_LF' 'Dist_FF'... 
                   'B0_LF-LL' 'B0_LL-LF' ...
                   'B0_FF-LL' 'B0_LL-FF' ...
                   'B0_FF-LF' 'B0_LF-FF'};      
second_interested={'pos_Age_p05_k15.csv','neg_Age_p05_k15.csv'};
atlas=spm_atlas('load','C:\Program Files\MATLAB\R2017a\spm12\atlas\AAL.xml');% Load in aal atlas for labeling during write csv

%% Read in Parametric second level peaks coordinate (age-modulating peaks)
% (Extract all peaks and sub peaks)
%% And then extract first levels beta within these peaks (beta which is the DV in second GLM)
% This correspond to 'first_cond' during the iteration(e.g., betas to extract for Dist_All-peaks is then 'Dist_All')

% Read in the common mask which will be overlaid with ROI 
% to create a combined ROI for extracting Beta later
V_mask=spm_vol('D:\Yun-Shiuan_LAMBDA\template\for_this_study\binary_overlapped_WM_GM_common_mask_n40_with_signal.nii'); 
                                                    % Need transfrom to mni space(because con_0001.nii is in mni space)
Y_mask=spm_read_vols(V_mask);
T_mask=V_mask.mat; %abtain the transformation rule to enable transform to mni space
c_mask=find(Y_mask>0.99);% select the voxel which mask value==1
[x y z]=ind2sub(size(Y_mask),c_mask);
Y_mni_mask=cor2mni([x y z],T_mask);%transform to mni space corrdinates

for first_cond=1:size(Cond_Parameters_folder_name,2)
    % Load in the First level image
    img_index=pad(cellstr(num2str([first_cond+1])),4,'left','0');
    [imgfilename,dirs]=spm_select('ExtFPListRec',...%FORMAT [files,dirs] = spm_select('ExtFPListRec',direc,filt,frames)    
                             path, strcat('.*con_',img_index,'.nii'),1);              
     for second_cond=1:size(second_interested,2)
                peaks=read_mixed_csv(strcat("D:\Yun-Shiuan_LAMBDA\Second_level_parametric_new\",...
                         Cond_Parameters_folder_name{first_cond},"\Result_table\",second_interested{second_cond}),',');
                table_peaks=cell2table(peaks(3:end,1:end),...
                        'VariableNames',matlab.lang.makeValidName(cellstr(strcat(peaks(1,:),"_",peaks(2,:)))));
                % Add Main peak/Sub peak info
                table_peaks.main_peak=~cellfun(@isempty,table_peaks.cluster_equivk);
                
                % Extract mni coordinates    
                mni_x=cellfun(@str2num,table_peaks.x_x);
                mni_y=cellfun(@str2num,table_peaks.x_y);
                mni_z=cellfun(@str2num,table_peaks.x_z_mm_);
                mni_xyz=[mni_x,mni_y,mni_z];
                main_peak=table_peaks.main_peak;                     
          for peak=1:size(mni_x,1) 
                    % Create a mask which overlay ROI with grey_white_common_matter
                    % -Calculate the distance and define those cordianates
                    % -with "d<radius" as overlaid coordinate (both mask and center is in mni space(voxel size=1x1x1mm))
                    radius=8;
                    % List all distance between each voxel in the mask to
                    % specific the peak
                    Dist_sq = (bsxfun(@minus,mni_xyz(peak,1),Y_mni_mask(:,1)).^2+...
                               bsxfun(@minus,mni_xyz(peak,2),Y_mni_mask(:,2)).^2+...
                               bsxfun(@minus,mni_xyz(peak,3),Y_mni_mask(:,3)).^2);
                    sphere_index=(find(Dist_sq<radius^2));%indexes for the coordinates in the sphere
                    %% The Critical Overlaying Step
                    ROI_overlaid=Y_mni_mask(sphere_index,:);  %Extract the coordinates in shpere and mask 
                    %% Extract first level beta from the overlaid ROI
                    ROI_overlaid=round(ROI_overlaid);
                    beta_peak_value_cond_roi=spm_summarise(imgfilename,ROI_overlaid',@mean);
                    
                    %% Write these values as .csv
                    % Directory to write files to
                    cd('D:\Yun-Shiuan_LAMBDA\Second_level_parametric_new\First_Betas_extract_from_Second_level_peaks');
                    % Creat come sames for the file name
                    first_level_beta=Cond_Parameters_real_name{first_cond};
                    second_level_peak=regexp(second_interested{second_cond},'^.*(?=_p\d)','match');
                    aal_name=spm_atlas('query',atlas, mni_xyz(peak,:)');
                    
                    filename=char(strcat(first_level_beta,'_from_',second_level_peak,...
                            '_pk_',num2str(peak),...
                            '_k_',table_peaks.cluster_equivk(peak),...
                            '_x',table_peaks.x_x(peak),...
                            '_y',table_peaks.x_y(peak),...
                            '_z',table_peaks.x_z_mm_(peak),...
                            '_aal_', aal_name,...
                            '_clfdr',table_peaks.cluster_p_FDR_corr_(peak),...
                            '_pkpfdr',table_peaks.peak_p_FDR_corr_(peak),...
                            '_pkt',table_peaks.peak_T(peak),...
                            '_pkpunc',table_peaks.peak_p_unc_(peak),...
                            '.csv'));
                    dlmwrite(filename,...
                            [beta_peak_value_cond_roi, id],'precision',9)  
              strcat('first: ',num2str(first_cond),...
                  '; second: ',num2str(second_cond),...
                  '; peak: ',num2str(peak))
          end

    end     
end