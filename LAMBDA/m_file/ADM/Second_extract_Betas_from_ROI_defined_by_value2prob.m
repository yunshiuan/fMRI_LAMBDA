%% (Part8-b)Second_level:Extract Betas from ROI(defined by value2prob)
%addpath('/usr/local/spm12/')Second_ROI_defined_by_value2prob
%spm_jobman('initcfg') -- before running the batch
path='/home/.bml/Data/Bank6/ADM-YunShiuan';
cd (strcat(path,'/Second_ROI_defined_by_value2prob/'));
%% Read Peaks info(file name="all_peak")
addpath(strcat(path,'/Scripts(m)_files'));% to enable using the function "read_mixed_csv" 
cond={'prob','mag','pxm',...
    'accept','gain','loss',...
    'reject','ungain','unloss'};
value={'hedonism_positive_peak','hedonism_negative_peak',...
    'security_positive_peak','security_negative_peak'};
all_peak={};

for v=1:4;
    for c_index=1:9;
peak=[];
C=read_mixed_csv(strcat(path,'/Second_ROI_defined_by_value2prob/',cond{c_index},'/result_table_p005_size15/',value{v},'.csv'),',');
peak=str2double(C(3:end,[5 11:14])); %read : the peak size,p-value,peak locations
all_peak{c_index,v}=peak; %collect all peak to a singel big cell(4 value effects x 9 conditions)
    end
end

%% Extract Betas of Prob (Peaks of 4 kind of value effects on Prob)
%GUI:imgfilename2=spm_select(Inf,'image','Select files...',{},strcat(path,'/First_design_matrix_with_contrastncovariate'),'.*con_0001.nii',1);
%code:FORMAT [files,dirs] = spm_select('ExtFPListRec',direc,filt,frames)
beta_peak_value_cond_roi={};
for c_index=1:9; %specify condition(contrast)
        [imgfilename,dirs]=spm_select('ExtFPListRec',...%FORMAT [files,dirs] = spm_select('ExtFPListRec',direc,filt,frames)
            strcat(path,'/First_design_matrix_with_contrastncovariate'),strcat('.*con_000',num2str(c_index),'.nii'),1);
        if ((c_index~=8)&(c_index~=9))
        imgfilename=imgfilename([1:8,10:16,18:20,22:size(imgfilename,1)],:); %exclude S031,S040, S044            
        else %sincei=21(S044 do not have condition ungain and unloss, and was already excluded during spm_select)
        imgfilename=imgfilename([1:8,10:16,18:size(imgfilename,1)],:); %exclude S031,S040          
        end     
    for v=1:4; %specify value
        
        for p=1:size(all_peak{c_index,v},1); %specify which peak in the results table
            beta_peak_value_cond_roi{p,v,c_index}=... %peaks, 4 value effects, and 9 cond
                [spm_summarise(imgfilename, struct('def', 'sphere', 'spec', 8, 'xyz', all_peak{c_index,v}(p,3:5)'), @mean),...
                str2num(imgfilename(:,83:84))];   % cbind the ID info  
         p  
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
cd(strcat(path,'/Second_ROI_defined_by_value2prob/Betas_extracted_from_ROI'))
for c_index=1:9;
    for v=1:4;
        for p=1:size(all_peak{c_index,v},1); %specify which peak in the results table 
            dlmwrite(strcat('betas_',cond{c_index},'_',value{v},num2str(p),'_size',...
                num2str(all_peak{c_index,v}(p,1)),'_p',num2str(all_peak{c_index,v}(p,2)),'_x',...
                num2str(all_peak{c_index,v}(p,3)),'_y',num2str(all_peak{c_index,v}(p,4)),'_z',...
                num2str(all_peak{c_index,v}(p,5)),'.csv'),...
               beta_peak_value_cond_roi{p,v,c_index},'precision',9)  
        end
        v
    end
    c_index
end