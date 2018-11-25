%% (Part9-a)Simple Slope Analysis(pxm interaction)-->Extrac Beta(parametric:prob,mag,pxm) from PxM peak(parametric)
%addpath('/usr/local/spm12/')
%spm_jobman('initcfg') -- before running the batch
path='/bml/Data/Bank6/ADM-YunShiuan';
cd (strcat(path,'/Second_ROI_defined_by_value2prob_new_server/'));
%% Retrieve pxm peaks
addpath(strcat(path,'/Scripts(m)_files'));% to enable using the function "read_mixed_csv" 
cond={'prob','mag','pxm',...
    'accept','gain','loss',...
    'reject','ungain','unloss','choice'};
value={'intercept_positive_peak','intercept_negative_peak',...
    'hedonism_positive_peak','hedonism_negative_peak',...
    'security_positive_peak','security_negative_peak'};
betatype={'prob','mag','pxm'};
all_peak={};
%Specify the v2pxm peaks location
for v=3:6;
    for c_index=3;
peak=[];
C=read_mixed_csv(strcat(path,'/Second_ROI_defined_by_value2prob_new_server/normalise_3x3x3/two_SVS_in_same_GLM/',cond{c_index},'/result_table_p005_size15/',value{v},'.csv'),',');
peak=str2double(C(3:end,[12:14])); %read : the peak size,p-value,peak locations
all_peak{c_index,v}=peak; 
    end
end
%% Extract Prob and Mag Betas from the peaks
beta_peak_value_cond_roi={};
%% Create grey_white mask which will be overlaid with ROI to create a combined ROI for extracting Beta later
 V_mask=spm_vol(strcat(path,'/smooth_mask_grey_white.nii'));% Need transfrom to mni space(because con_0001.nii is in mni space)
 Y_mask=spm_read_vols(V_mask);
 Y_mask=Y_mask>0.8;%Set criteria >0.8 as brain(white/grey) region
 T_mask=V_mask.mat; %abtain the transformation rule to enable transform to mni space
 c_mask=find(Y_mask==1);% select the voxel which mask value==1
 [x y z]=ind2sub(size(Y_mask),c_mask);
 Y_mni_mask=cor2mni([x y z],T_mask);%transform to mni space corrdinates
 
for c_index=3; %specify condition(pxm peaks location)
    for beta= 1:3; % 1 for prob; 2 for mag; 3 for pxm (extract the irst level Betas in the v2p location)
        [imgfilename,dirs]=spm_select('ExtFPListRec',...%FORMAT [files,dirs] = spm_select('ExtFPListRec',direc,filt,frames)
            strcat(path,'/First_design_matrix_with_contrastncovariate_parametric'),strcat('.*con_000',num2str(beta),'.nii'),1);
            imgfilename_index=find(cell2mat(cellfun(@(x) ~isempty(regexp(x,'3x3x3','match')),cellstr(imgfilename),'un',0)));
           imgfilename=imgfilename(imgfilename_index,:);%only select those with  normarlise_3x3x3
           imgfilename=imgfilename([1:8,10:16,18:size(imgfilename,1)],:); %exclude S031,S040  %exclude S031,S040(S44 has already been skipped in the first level)            
    for v=3:6; %specify value type of peak region(v2p peaks)       
        for p=1:size(all_peak{c_index,v},1); %speall_peakcify which peak in the results table
            if size(all_peak{c_index,v},1)==0 %skip if this condition has no peaks
                {};
            else
         %% Create a mask which overlay ROI with grey_white_matter
                %(Used): Calculate the distance and define those cordianates
              % with "d<radius" as overlaid coordinate (both mask and center is in mni space(voxel size=1x1x1mm))
              center=all_peak{c_index,v}(p,end-2:end);%specify the peak locations
              radius=8;
              Dist_sq = (bsxfun(@minus,center(1,1),Y_mni_mask(:,1)).^2+...
                      bsxfun(@minus,center(1,2),Y_mni_mask(:,2)).^2+...
                      bsxfun(@minus,center(1,3),Y_mni_mask(:,3)).^2);...
              sphere_index=(find(Dist_sq<radius^2));%indexes for the coordinates in the sphere
          %% The Critical Overlaying Step
            ROI_overlaid=Y_mni_mask(sphere_index,:);  %Extract the coordinates in shpere and mask 
          %% Extract voxel signals from the overlaid ROI
            ROI_overlaid=round(ROI_overlaid);
            beta_peak_value_cond_roi{p,v,c_index,beta}=... %peaks, 4 value effects, and 1 cond
                spm_summarise(imgfilename,ROI_overlaid',@mean);
%             beta_peak_value_cond_roi{p,v,c_index,beta}=... %peaks, 2 value effects, and 1 cond
%                 [spm_summarise(imgfilename, struct('def', 'sphere', 'spec', 8, 'xyz', all_peak{c_index,v}(p,end-2:end)'), @mean),...
%                 ];   % cbind the ID info in tall_peakhe end(when writing csv file) 
         p  
            end
        end
    end
         beta
         v
    end
         c_index
end
%% Output as csv
id=csvread(strcat(path,'/young_idlist.csv')); %read only the young's IDs
id=setdiff(id,[31 40 44]); %Remove ID = 31 40 44
cd(strcat(path,'/Second_ROI_defined_by_value2prob_new_server/normalise_3x3x3/two_SVS_in_same_GLM/Betas_of_p_m_extracted_from_pm_peaks_overlaid'))
for c_index=3;
    for v=3:6;
       for beta= 1:3; % 1 for prob; 2 for mag; 3 for pxm
        for p=1:size(all_peak{c_index,v},1); %specify which peak in the results table 
            dlmwrite(strcat('betas_',betatype{beta},'_extracted_from_pm_peaks',...
            '_x',num2str(all_peak{c_index,v}(p,end-2)),...
            '_y',num2str(all_peak{c_index,v}(p,end-1)),...
            '_z',num2str(all_peak{c_index,v}(p,end)),'.csv'),...
              [beta_peak_value_cond_roi{p,v,c_index,beta}, id],'precision',9)  
            end
        v
        end
    beta
    c_index
    end
end

%% Save all Overlaid ROI for displaying
cd(strcat(path,'/Second_ROI_defined_by_value2prob_new_server/normalise_3x3x3/two_SVS_in_same_GLM/Betas_of_p_m_extracted_from_pm_peaks_overlaid/ROI_overlaid_nii_file'))
for c_index=3; %specify condition(pxm peaks location)
    for beta= 1:15; % 15 conditions
           for v=3:6; %specify v2p peaks      
        for p=1:size(all_peak{c_index,v},1); %specify which peak in the results table
            if size(all_peak{c_index,v},1)==0 %skip if this condition has no peaks
                {};
            else
            %% Create a mask which overlay ROI with grey_white_matter
%         %Y_mni_mask %the mask of grey and white matter;
              %(Used): Calculate the distance and define those cordianates
              % with "d<radius" as overlaid coordinate (both mask and center is in mni space(voxel size=1x1x1mm))
              center=all_peak{c_index,v}(p,end-2:end);
              radius=8;
              Dist_sq = (bsxfun(@minus,center(1,1),Y_mni_mask(:,1)).^2+...
                      bsxfun(@minus,center(1,2),Y_mni_mask(:,2)).^2+...
                      bsxfun(@minus,center(1,3),Y_mni_mask(:,3)).^2);...
              sphere_index=(find(Dist_sq<radius^2));%indexes for the coordinates in the sphere
              %% The Critical Overlaying Step
              ROI_overlaid=Y_mni_mask(sphere_index,:);  %Extract the coordinates in shpere and mask 
              %Save as nii. for display to check if the ROI_overlaid is right-->Correct
                 vol=V_mask;
                 file_name=char(strcat(path,'/Second_ROI_defined_by_value2prob_new_server/normalise_3x3x3/two_SVS_in_same_GLM/Betas_of_p_m_extracted_from_pm_peaks_overlaid/ROI_overlaid_nii_file/roi_',...
                                num2str(center),'.nii'));
                 vol.fname=file_name;
                 %convert Y_mni_mask to matrix form(voxle-based coordinates)
                 img=NaN([121,145,121]);
                 Y_cor_roi=mni2cor([round(ROI_overlaid(:,1)),round(ROI_overlaid(:,2)), round(ROI_overlaid(:,3))],T_mask);
                 img(sub2ind([121,145,121], Y_cor_roi(:,1) , Y_cor_roi(:,2), Y_cor_roi(:,3)))=1;
                 spm_write_vol(vol,img);                         
         p 
            end
        end
    end
         beta
         v
    end
         c_index
end


%% Get the HRF transformed Mag./Prob. Level (H/M/L)--> NO NEED!!!(Since HRF(cxProb)=cxHRF(Prob), no need to calculate HRF(c))
% path='/bml/Data/Bank6/ADM-YunShiuan';
% cd (strcat(path,'/First_design_matrix_with_contrastncovariate/'));
% %Get young id
% info=load(strcat(path,'/id_age_gender_sec_hed_sti.mat'));
% index=logical(cell2mat(cellfun(@(x) ~(isequal(x,'S031')|isequal(x,'S040')|isequal(x,'S044')),info.a(:,1),'un',0)));%exclude i=(S031)(S040)21(S044) who have excessive motions
% info.a=info.a(index,:);
% %id
% young_name=info.a(:,1);
% %Acces the HRF transformed value from each person's design matrix
% %Raw Mag, Prob, PxM
% %% Get the raw Mag/Prob value
% %Get the raw value to enable matching with design matrix value
%   A2={0,0,0,0,0};% save info to "A2"
% for j=1:5;
%     fid = fopen(char(strcat(path,'/Run_info(csv)/S023/r',num2str(j),'.csv')));  %open file
%  %read all contents into data as a char array 
%     A = textscan(fid,'%s','Delimiter', ',');
%     A2{1,j}=reshape(A{1},8,90)';
% end
% prob_1=str2double(A2{1,1}(1:2:length(A2{1,1}),3));
% prob_2=str2double(A2{1,2}(1:2:length(A2{1,1}),3));
% prob_3=str2double(A2{1,3}(1:2:length(A2{1,1}),3));
% prob_4=str2double(A2{1,4}(1:2:length(A2{1,1}),3));
% prob_5=str2double(A2{1,5}(1:2:length(A2{1,1}),3));
% %prob_mean=mean(mean([prob_1,prob_2,prob_3,prob_4,prob_5]));
% %prob_center=cellfun(@(x) (x-prob_mean),{prob_1,prob_2,prob_3,prob_4,prob_5},'un',0);
% prob=[prob_1;prob_2;prob_3;prob_4;prob_5];
% 
% mag_1=str2double(A2{1,1}(1:2:length(A2{1,1}),4));
% mag_2=str2double(A2{1,2}(1:2:length(A2{1,1}),4));
% mag_3=str2double(A2{1,3}(1:2:length(A2{1,1}),4));
% mag_4=str2double(A2{1,4}(1:2:length(A2{1,1}),4));
% mag_5=str2double(A2{1,5}(1:2:length(A2{1,1}),4));
% %mag_mean=mean(mean([mag_1,mag_2,mag_3,mag_4,mag_5]));
% %mag_center=cellfun(@(x) (x-mag_mean),{mag_1,mag_2,mag_3,mag_4,mag_5},'un',0);
% mag=[mag_1;mag_2;mag_3;mag_4;mag_5];
% 
% probxmag=cell2mat(prob_center).*cell2mat(mag_center);
% probxmag=mat2cell(probxmag,45,[1 1 1 1 1]);
% %Checking: The raw centered Mag.(Identical to the centered value above)
% %SPM.Sess(1).U(1).P(2).P
% %% (Failed)%% Approach 1 : Retrive the HRF transformed value from design matrix
% % cd (strcat(path,'/First_design_matrix_with_contrastncovariate/'));
% % %Identify the Mag=L/M/H trial
% % % figure()
% % % plot(SPM.xX.X(1:218,3)) % The HRF signal of Mag
% % find(mag==7,1)%Mag=L trial 12-->
% % find(mag==56,1)%Mag=M trial 3-->HRF= 1.0450
% % find(mag==105,1)%Mag=H trial 32-->
% % %Fail to capture all the 45 trials HRF signal
% 
% % n_image=218;
% % trial_index=[];
% % for i = 1:n_image;
% % hrf_this=SPM.xX.X(i,3);
% % hrf_next1=SPM.xX.X(i+1,3);
% % hrf_next2=SPM.xX.X(i+2,3);
% % hrf_next3=SPM.xX.X(i+3,3);
% % hrf_next4=SPM.xX.X(i+4,3);
% % hrf_next5=SPM.xX.X(i+5,3);
% % hrf_next6=SPM.xX.X(i+6,3);
% % delta_this_next1=hrf_next1-hrf_this;
% % delta_next1_next2=hrf_next2-hrf_next1;
% % delta_next2_next3=hrf_next3-hrf_next2;
% % delta_next3_next4=hrf_next4-hrf_next3;
% % delta_next4_next5=hrf_next5-hrf_next4;
% % delta_next5_next6=hrf_next6-hrf_next5;
% %     if (delta_next2_next3*delta_next3_next4<0)&&...
% %         (delta_this_next1*delta_next1_next2>0)&&(delta_next1_next2*delta_next2_next3>0)&&...
% %     (delta_next3_next4*delta_next4_next5>0)&&(delta_next4_next5*delta_next5_next6>0)
% %     trial_index=[trial_index i+3];
% %     end
% % end
% 
% %% (Used) Approach 2: Create fake BOLD signal and transform by HRF
% %Confirm that c as a constant, then，HRF(axc)=c x HRF(a).
% %Therefore, HRF(Prob x Mag), given Mag=c, then HRF(Prob x Mag)= 
% %c x HRF(Prob)!! --> Could use HRF(Mag=c) to do simple slope analysis
% % bf = spm_get_bf;
% % ran= 10*rand(5,1);
% % reg = [zeros(2,1); ran; zeros(2,1); ran; zeros(2,1)];
% % U.u = reg;
% % U.name = {'reg'};
% % convreg = spm_Volterra(U, bf.bf);
% % figure()
% % plot(convreg)
% % 
% % reg10 = [zeros(2,1); 10*ran; zeros(2,1); 10*ran; zeros(2,1)];
% % U10.u=reg10;
% % U10.name = {'reg'};
% % convreg10 = spm_Volterra(U10, bf.bf);
% % figure()
% % plot(convreg10)
% %Calculate the specified Mag Level
%  bf = spm_get_bf;%Time =2
%  mag_level_raw=[7 56 105];
%  mag_level_hrf=[];
%  
%  for i=1:3;
%  reg = [zeros(2,1); mag_level_raw(i)-mag_mean; zeros(20,1)];
%  U.u = reg;
%  U.name = {'reg'};
%  convreg = spm_Volterra(U, bf.bf);
% %  figure()
% %  plot(convreg)
%  mag_hrf=max(convreg);
%  mag_level_hrf=[mag_level_hrf mag_hrf];
%  end
% 
% %   reg = [zeros(2,1); 55-mag_mean; zeros(20,1)];
% %  U.u = reg;
% %  U.name = {'reg'};
% %  convreg = spm_Volterra(U, bf.bf);
% %  figure()
% %  plot(convreg)
% %  mag_hrf=max(convreg);
% %-0.0719