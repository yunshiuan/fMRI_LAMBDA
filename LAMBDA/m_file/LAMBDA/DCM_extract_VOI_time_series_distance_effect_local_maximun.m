%% Part1-a-ii: DCM - Extract time series in VOIs (for distance effect DCM)
%% NOTE: Different VOI for all participants (Movement of center per each participant's local maximun)
%%       First, create sphere based on 2nd level overall negative distsance sensitivity result peaks (collapsed age).
%%       Then adjust the center of the sphere to each's nearest local maximun shown by 1st level N-F contrast (which ignores notation and mimics overall distance sensitivity)

addpath(char("D:\GoogleDrive\Lambda_code\m_file\LAMBDA\tool_code"));% enable read_mixed_csv()
%spm_jobman('initcfg') -- before running the batch

path='D:\Yun-Shiuan_LAMBDA\DCM\Part0_First_GLM_for_DCM_single_run'; % To get the subject list with valid runs 
cd(path);
subject_list=cellstr(ls('*')); % subject list (which has valid first level contrast, i.e., not being excluded)
subject_list=subject_list(~cellfun(@isempty,regexp(subject_list,'(?<=df)\d+','match')));% include only subject name

%% Read in run inclusion index info
run_inclusion_index=cellfun(@(x) regexprep(x,'"',''),...
                read_mixed_csv("D:\Yun-Shiuan_LAMBDA\Run_inclusion_info\inclusive_runs_indexes.csv",','),'un',0);
run_inclusion_index=table(run_inclusion_index(2:end,2),run_inclusion_index(2:end,3),...
                    'VariableNames',{'sub_id','run_num'});
%% Some constants
% 8 regions invloving in distance sensitivity
VOI_name={"R_V1","L_V1",...
    "R_IPS","L_IPS",...
    "R_SMA","L_SMA",...
    "L_M1","R_VLPFC"}';
peak_coordinates={...
    [21 -100 -4], [-24,-97,-7],...
    [45 -40 47],[-42 -40 41],...
    [6 17 50],[-6 11 50],...
    [-42 -1 32],[48 35 23]};
common_mask='D:\Yun-Shiuan_LAMBDA\template\for_this_study\binary_overlapped_WM_GM_common_mask_n40_with_signal.nii,1';
%% Batch for VOI Extraction
for i=1:size(subject_list,1)
    % Read in the inclusion run indexes and use them to index runs in the following for loop
    run_rows=strcmp(run_inclusion_index.sub_id,subject_list{i});        
    run_valid=table2cell(run_inclusion_index(run_rows,'run_num'));
    run_valid=regexp(run_valid,'(?<=r)\d+','match');
    run_valid=str2double([run_valid{:}]');
    
    % Make directory (for VOI output)
    outputdir=strcat('D:\Yun-Shiuan_LAMBDA\DCM\Part1_VOIs_local_maximun_for_DCM\',subject_list{i});
    mkdir(outputdir);
    outputdir_batch=strcat(outputdir,'\VOI_batch');
    mkdir(outputdir_batch);
    outputdir_data=strcat(outputdir,'\VOI_data');
    mkdir(outputdir_data);

    for run=1:length(run_valid)
        r=run_valid(run);% The index refer to the real run number of the task (so it might not be consecutive due to exclusive runs)
        for v=1:size(VOI_name,1)
            clear matlabbatch
            % Make batch for VOI time series extraction (per subject run)
            matlabbatch{1}.spm.util.voi.spmmat = {strcat('D:\Yun-Shiuan_LAMBDA\DCM\Part0_First_GLM_for_DCM_single_run\',...
                         subject_list{i},'\run',num2str(r),'\SPM.mat')};           
            matlabbatch{1}.spm.util.voi.adjust = 1;
            matlabbatch{1}.spm.util.voi.session = 1;
            matlabbatch{1}.spm.util.voi.name = char(strcat(VOI_name{v},'_run',num2str(r)));
            % Common mask in which every participant has signal
            matlabbatch{1}.spm.util.voi.roi{1}.mask.image = {'D:\Yun-Shiuan_LAMBDA\template\for_this_study\binary_overlapped_WM_GM_common_mask_n40_with_signal.nii,1'};
            matlabbatch{1}.spm.util.voi.roi{1}.mask.threshold = 0.9;
            
            % Fisrt level result for each participant
            matlabbatch{1}.spm.util.voi.roi{2}.spm.spmmat = {strcat('D:\Yun-Shiuan_LAMBDA\DCM\Part0_First_GLM_for_DCM_multiple_run\',...
                                                                subject_list{i},'\SPM.mat')};
            matlabbatch{1}.spm.util.voi.roi{2}.spm.contrast = 8; % The N_F T contrast
            matlabbatch{1}.spm.util.voi.roi{2}.spm.conjunction = 1;
            matlabbatch{1}.spm.util.voi.roi{2}.spm.threshdesc = 'none';
            matlabbatch{1}.spm.util.voi.roi{2}.spm.thresh = 0.05;
            matlabbatch{1}.spm.util.voi.roi{2}.spm.extent = 0;
            matlabbatch{1}.spm.util.voi.roi{2}.spm.mask = struct('contrast', {}, 'thresh', {}, 'mtype', {});
            % Sphere which center is from 2nd level result. The center
            % position is goin to be adjusted according to nearest local maximun from the 
            % 1st level result above.
            matlabbatch{1}.spm.util.voi.roi{3}.sphere.centre = peak_coordinates{v};
            matlabbatch{1}.spm.util.voi.roi{3}.sphere.radius = 8;
            matlabbatch{1}.spm.util.voi.roi{3}.sphere.move.local.spm = 2;
            matlabbatch{1}.spm.util.voi.roi{3}.sphere.move.local.mask = '';
            
            % VOI defined by the intersection (overlay center-adjusted shpere and 1ST level significant cluster shape)
            matlabbatch{1}.spm.util.voi.expression = 'i1 & i2 & i3';
            % Save the bacth
            save(strcat(outputdir_batch,'\',VOI_name{v},'_run',num2str(r)),'matlabbatch') 
            
            %Run the batch
            spm_jobman('run',matlabbatch);
            
            % Move the exracted data
            sourcefiles=strcat('D:\Yun-Shiuan_LAMBDA\DCM\Part0_First_GLM_for_DCM_single_run\',subject_list{i},'\run',num2str(r));
            movefile(strcat(sourcefiles,'\VOI*'),outputdir_data);
        end
    end
end

%% Check the VOI amount (see if everyone has 8 VOI per run)
% Read in Gender and Age to see if the valid amount correlates with age
demographic=cellfun(@(x) regexprep(x,'"',''),...
                read_mixed_csv("D:\Yun-Shiuan_LAMBDA\demographic\tidy_demographic.csv",','),'un',0);
table_demopraphic=cell2table(demographic(2:end,2:end),...
                         'VariableNames',cellstr(demographic(1,2:end)));
table_demopraphic=sortrows(table_demopraphic,'sub_id','ascend'); % Need to sort bi subID to make sure the order is ascending
valid_Gender_Dummy=str2double(table2array(table_demopraphic(ismember(lower(table_demopraphic.sub_id),subject_list),'gender_dummy')));
valid_Age_Precise=str2double(table2array(table_demopraphic(ismember(lower(table_demopraphic.sub_id),subject_list),'scan_age_precise')));

check_amount={};
for i=1:size(subject_list,1)
    % Read in the inclusion run indexes and use them to index runs in the following for loop
    run_rows=strcmp(run_inclusion_index.sub_id,subject_list{i});        
    run_valid=table2cell(run_inclusion_index(run_rows,'run_num'));
    run_valid=regexp(run_valid,'(?<=r)\d+','match');
    run_valid=str2double([run_valid{:}]');
    
    run_valid_amount=size(run_valid,1);
    VOI_amount_predicted=size(peak_coordinates,2);
    amount_predicted=run_valid_amount* VOI_amount_predicted;
    
    amount_created=size(ls(strcat('D:\Yun-Shiuan_LAMBDA\DCM\Part1_VOIs_local_maximun_for_DCM\',subject_list{i},'\VOI_data\*mat')),1);
    VOI_amount_created=amount_created/run_valid_amount;
    
    check_amount{i,1}=subject_list{i};
    check_amount{i,2}=amount_predicted;
    check_amount{i,3}=amount_created;
    check_amount{i,4}=VOI_amount_created;
end
table_check_amount=cell2table(check_amount,'VariableNames',{'ID','amount_predicted','amount_created','VOI_per_run_created'});
table_check_amount.age=valid_Age_Precise;
sortrows(table_check_amount,'VOI_per_run_created','descend') % Doesn't seem to be an age-specific issue

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Remained to ......
%(1)See if the divergence is correlated with age
%(2)See if "masking" the thresholded SPM 1st level result could help (to limit the local maximun search space)
%(3)After resolving the out-of-scope issue by (2), check if the centers locate in reasonable places.
%%%%%%%If not,
%%%%%%%(3-1) Reconsider if local maximun adjustment is needed
%%%%%%%(3-2) Try a different local maximun finding algorithm
%% Check the VOI positions
% Collect adjusted centers (40 ids x 8 VOIs x 3 mni coordinates)
cd('D:\Yun-Shiuan_LAMBDA\DCM\Part1_VOIs_local_maximun_for_DCM')
warning_check=cell(length(subject_list),1);
collect_centers=nan(size(subject_list,1),size(VOI_name,1),3); % 40 ids x 8 VOIs x 3 mni coordinates
for v=1:size(VOI_name,1)
    for i=1:size(subject_list,1)
        VOI_dir=strcat('D:\Yun-Shiuan_LAMBDA\DCM\Part1_VOIs_local_maximun_for_DCM\',subject_list{i},'\VOI_data\');
        cd(VOI_dir)
        mat_list=cellstr(ls('*mat'));
        mat_name=mat_list(~cellfun(@isempty,regexp(mat_list,VOI_name{v},'match'))); %list VOIs from all runs of the specific ROI for this subject 
        if isempty(mat_name)
            warning_check{i,1}=strcat('None valid VOI was found.');
        else
            load(mat_name{1},'xY') %Only need to load in any run (center locations for every runs per subject are identical)
            collect_centers(i,v,1)=xY.xyz(1);
            collect_centers(i,v,2)=xY.xyz(2);
            collect_centers(i,v,3)=xY.xyz(3);
        end
    end
end

valid_Age_Precise=str2double(table2array(table_demopraphic(ismember(lower(table_demopraphic.sub_id),subject_list),'scan_age_precise')));

% Visulaize
valid_grade=str2double(table2array(table_demopraphic(ismember(lower(table_demopraphic.sub_id),subject_list),'grade')));
for v=1:size(VOI_name,1)
    figure
    scatter3(collect_centers((valid_grade==5),v,1),collect_centers((valid_grade==5),v,2),collect_centers((valid_grade==5),v,3),...
        30,'blue','filled')
    hold on
    scatter3(collect_centers((valid_grade==2),v,1),collect_centers((valid_grade==2),v,2),collect_centers((valid_grade==2),v,3),...
        30,[255/256,153/256,102/256],'filled')
    title({char(strcat(VOI_name{v},' : ',num2str(peak_coordinates{v}))),...
           char(strcat('Participants with valid coordinates within the brain area : ',...
                num2str(sum(~isnan(collect_centers(:,v,1)))))),...
           char(strcat('Valid: 5th: ',num2str(sum(~isnan(collect_centers(valid_grade==5,v,1)))),'/',num2str(sum(valid_grade==5)),...
                ' ;2nd: ',num2str(sum(~isnan(collect_centers(valid_grade==2,v,1)))),'/',num2str(sum(valid_grade==2))))})   
    scatter3(peak_coordinates{v}(1),peak_coordinates{v}(2),peak_coordinates{v}(3),100,'r','Marker','x')
%     text(collect_centers(:,v,1)+0.5,collect_centers(:,v,2)+2,collect_centers(:,v,3)+2,cellstr(num2str(valid_grade)))
    legend('5 grader', '2 graders')
    
    file_name=char(strcat('D:\Yun-Shiuan_LAMBDA\DCM\Part1_VOIs_local_maximun_for_DCM\locations_adjusted_by_local_max\',VOI_name{v},'.png'))
    saveas(gcf,file_name)
    close(gcf)
end
% load('VOI_L_IPS_run1_1.mat')
% run1_center=xY.xyz;
% run1_signal=Y;
% load('VOI_L_IPS_run2_1.mat')
% run2_center=xY.xyz;
% run2_signal=Y;
% [run1_center run2_center]
% [run1_signal run2_signal]

% Running job #2
% ------------------------------------------------------------------------
% Running 'Volume of Interest'
% 	SPM computation                 :                        ...done
% 	SPM computation                 :                        ...done
%    Sphere centre moved from [  6  17  50] to [ 48  47  11]
%    VOI saved as D:\Yun-Shiuan_LAMBDA\DCM\Part0_First_GLM_for_DCM_single_run\df1046\run2\VOI_R_SMA_run2_1.mat
% Done    'Volume of Interest'
% Done

% ------------------------------------------------------------------------
% Running job #2
% ------------------------------------------------------------------------
% Running 'Volume of Interest'
% 	SPM computation                 :                        ...done
% 	SPM computation                 :                        ...done
%    Sphere centre moved from [-42 -40  41] to [ -6 -82 -28]
% Warning: Empty region. 
% > In spm_regions (line 155)
%   In spm_run_voi (line 69)
%   In cfg_run_cm (line 29)
%   In cfg_util>local_runcj (line 1688)
%   In cfg_util (line 959)
%   In spm_jobman>fill_run_job (line 469)
%   In spm_jobman (line 247) 
% Done    'Volume of Interest'
% Done
% 

