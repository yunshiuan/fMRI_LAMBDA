%% (Part8-b)Second_level:Extract Betas from ROI(defined by value2prob)e.g. Beta(Prob.) from Prob. peaks
%addpath('/usr/local/spm12/')
%spm_jobman('initcfg') -- before running the batch
path='/bml/Data/Bank6/ADM-YunShiuan';
cd (strcat(path,'/Second_ROI_defined_by_value2prob_new_server/normalise_3x3x3/two_SVS_in_same_GLM'));
%% Read Peaks info(file name="all_peak")
addpath(strcat(path,'/Scripts(m)_files'));% to enable using the function "read_mixed_csv" 
cond={'prob','mag','pxm',...
    'accept','gain','loss',...
    'reject','ungain','unloss','choice'};
value={'intercept_positive_peak','intercept_negative_peak',...
    'hedonism_positive_peak','hedonism_negative_peak',...
    'security_positive_peak','security_negative_peak'};
all_peak={};

for v=1:6;
    for c_index=1:10;
peak=[];
C=read_mixed_csv(strcat(path,'/Second_ROI_defined_by_value2prob_new_server/normalise_3x3x3/two_SVS_in_same_GLM/',cond{c_index},'/result_table_p005_size15/',value{v},'.csv'),',');
peak=str2double(C(3:end,[4,5,8,9,11:14])); %read : the peak size,p-value,peak locations
all_peak{c_index,v}=peak; %collect all peak to a singel big cell(6 value effects x 9 conditions)
    end
end

%% Extract Betas of Prob (Peaks of 2 intercepts and 4 kind of value effects on Prob)
%GUI:imgfilename2=spm_select(Inf,'image','Select files...',{},strcat(path,'/First_design_matrix_with_contrastncovariate'),'.*con_0001.nii',1);
%code:FORMAT [files,dirs] = spm_select('ExtFPListRec',direc,filt,frames)
beta_peak_value_cond_roi={};
%% Create grey_white mask which will be overlaid with ROI to create a combined ROI for extracting Beta later
 V_mask=spm_vol(strcat(path,'/smooth_mask_grey_white.nii'));% Need transfrom to mni space(because con_0001.nii is in mni space)
 Y_mask=spm_read_vols(V_mask);
 Y_mask=Y_mask>0.8;%Set criteria >0.8 as brain(white/grey) region
 T_mask=V_mask.mat; %obtain the transformation rule to enable transform to mni space
 c_mask=find(Y_mask==1);% select the voxel which mask value==1
 [x y z]=ind2sub(size(Y_mask),c_mask);
 Y_mni_mask=cor2mni([x y z],T_mask);%transform to mni space corrdinates
 
con_num_index=strcat(num2str(0),num2str([1:9]'));
con_num_index=cellstr([con_num_index;num2str([10]')]);
for c_index=1:10; %specify condition(contrast)
        [imgfilename,dirs]=spm_select('ExtFPListRec',...%FORMAT [files,dirs] = spm_select('ExtFPListRec',direc,filt,frames)
            strcat(path,'/First_design_matrix_with_contrastncovariate_parametric/'),strcat('.*con_00',con_num_index{c_index},'.nii'),1);              
        imgfilename_index=find(cell2mat(cellfun(@(x) ~isempty(regexp(x,'3x3x3','match')),cellstr(imgfilename),'un',0)));
        imgfilename=imgfilename(imgfilename_index,:);%only select those with  normarlise_3x3x3
        imgfilename=imgfilename([1:8,10:16,18:size(imgfilename,1)],:); %exclude S031,S040  %exclude S031,S040(S44 has already been skipped in the first level)            
    for v=1:6; %specify value        
          for p=1:size(all_peak{c_index,v},1); %specify which peak in the results table
            if size(all_peak{c_index,v},1)==0 %skip if this condition has no peaks
                {};
            else
            %% Create a mask which overlay ROI with grey_white_matter
                %(Used): Calculate the distance and define those cordianates
              % with "d<radius" as overlaid coordinate (both mask and center is in mni space(voxel size=1x1x1mm))
              center=all_peak{c_index,v}(p,end-2:end);
              radius=8;
              Dist_sq = (bsxfun(@minus,center(1,1),Y_mni_mask(:,1)).^2+...
                      bsxfun(@minus,center(1,2),Y_mni_mask(:,2)).^2+...
                      bsxfun(@minus,center(1,3),Y_mni_mask(:,3)).^2);
              sphere_index=(find(Dist_sq<radius^2));%indexes for the coordinates in the sphere
          %% The Critical Overlaying Step
            ROI_overlaid=Y_mni_mask(sphere_index,:);  %Extract the coordinates in shpere and mask 
             %% Extract voxel signals from the overlaid ROI
            ROI_overlaid=round(ROI_overlaid);
            beta_peak_value_cond_roi{p,v,c_index}=... %peaks, 6 value effects, and 9 cond
                spm_summarise(imgfilename,ROI_overlaid',@mean);
%             beta_peak_value_cond_roi{p,v,c_index}=... %peaks, 6 value effects, and 9 cond
%                 [spm_summarise(imgfilename, struct('def', 'sphere', 'spec', 8, 'xyz', all_peak{c_index,v}(p,end-2:end)'), @mean),...
%                 ];   % cbind the ID info in the end(when writing csv file) 
         p  
        end
        end
         v
    end
         c_index
end
%check accuracy
% nrow={};ncol={};
% for i=1:80;
% for j=1:4
%     for k=1:9;
%     nrow{i,j,k}=size(beta_peak_value_cond_roi{i,j,k},1);
%     ncol{i,j,k}=size(beta_peak_value_cond_roi{i,j,k},2);
%     end
% end
% end
%Found that: S044(i=21) do not have contrast:Ungain and Unloss(since he do not expereince them)

%check accuracy
%beta_peak_value_cond_roi{1,3,9}(:,1)==...
%    spm_summarise(imgfilename, struct('def', 'sphere', 'spec', 8, 'xyz', [50 -30 -10]'), @mean)
%% Output as csv
id=csvread(strcat(path,'/young_idlist.csv')); %read only the young's IDs
id=setdiff(id,[31 40 44]); %Remove ID = 31 40 44
cd(strcat(path,'/Second_ROI_defined_by_value2prob_new_server/normalise_3x3x3/two_SVS_in_same_GLM/Betas_extracted_from_ROI_overlaid'))
for c_index=1:10;
    for v=1:6;
        for p=1:size(all_peak{c_index,v},1); %specify which peak in the results table 
            dlmwrite(strcat('betas_',cond{c_index},'_',value{v},num2str(p),...
            '_size',num2str(all_peak{c_index,v}(p,2)),...
            '_clfdr',num2str(all_peak{c_index,v}(p,1)),...
            '_pkpfdr',num2str(all_peak{c_index,v}(p,3)),...
            '_pkt',num2str(all_peak{c_index,v}(p,4)),...
            '_pkpunc',num2str(all_peak{c_index,v}(p,5)),...
            '_x',num2str(all_peak{c_index,v}(p,end-2)),...
            '_y',num2str(all_peak{c_index,v}(p,end-1)),...
            '_z',num2str(all_peak{c_index,v}(p,end)),'.csv'),...
              [beta_peak_value_cond_roi{p,v,c_index}, id],'precision',9)  
        end
        v
    end
    c_index
end