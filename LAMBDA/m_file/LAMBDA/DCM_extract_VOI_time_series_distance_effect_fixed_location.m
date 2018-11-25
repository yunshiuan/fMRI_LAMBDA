%% Part1-a-i: DCM - Extract time series in VOIs (for distance effect DCM) 
%% NOTE: Fixed VOI for all participants (No movement of center per each participant's local maximun.) 
%%       Rather, purely based on 2nd level overall negative distsance sensitivity result peaks(collapsed age) to create sphere for extraction.

addpath(char("D:\GoogleDrive\Lambda_code\m_file\LAMBDA\tool_code"));% enable read_mixed_csv()
%spm_jobman('initcfg') -- before running the batch

path='D:\Yun-Shiuan_LAMBDA\DCM\Part0_First_GLM_for_DCM_single_run_with_reslice'; % To get the subject list with valid runs 
cd(path);
subject_list=cellstr(ls('*')); % subject list (which has valid first level contrast, i.e., not being excluded)
subject_list=subject_list(~cellfun(@isempty,regexp(subject_list,'(?<=df)\d+','match')));% include only subject name

%% Read in run inclusion index info
run_inclusion_index=cellfun(@(x) regexprep(x,'"',''),...
                read_mixed_csv("D:\Yun-Shiuan_LAMBDA\Run_inclusion_info\inclusive_runs_indexes.csv",','),'un',0);
run_inclusion_index=table(run_inclusion_index(2:end,2),run_inclusion_index(2:end,3),...
                    'VariableNames',{'sub_id','run_num'});
%% Constants
% 8 regions invloving in distance sensitivity
VOI_name={"R_V1","L_V1",...
    "R_IPS","L_IPS",...
    "R_SMA","L_SMA",...
    "L_M1","R_VLPFC"}';
peak_coordinates={...
    [30 -90 -2], [-26,-94,-4],...
    [28 -68 40],[-26 -56 48],...
    [6 22 46],[-6 12 48],...
    [-44 0 32],[48 34 22]};
common_mask='D:\Yun-Shiuan_LAMBDA\template\for_this_study\binary_overlapped_WM_GM_common_mask_n40_with_signal_with_reslice.nii,1';
%% Batch for VOI Extraction
radius=4; % Set the radius of the shere
for i=1:size(subject_list,1)
    % Read in the inclusion run indexes and use them to index runs in the following for loop
    run_rows=strcmp(run_inclusion_index.sub_id,subject_list{i});        
    run_valid=table2cell(run_inclusion_index(run_rows,'run_num'));
    run_valid=regexp(run_valid,'(?<=r)\d+','match');
    run_valid=str2double([run_valid{:}]');
    
    % Make directory (for VOI output)
    outputdir=strcat('D:\Yun-Shiuan_LAMBDA\DCM\Part1_VOIs_fixed_ROI_for_DCM_with_reslice\Houde_resliced_peak_4mm\',subject_list{i});
    outputdir_batch=strcat(outputdir,'\VOI_batch');
    mkdir(outputdir_batch);
    outputdir_data=strcat(outputdir,'\VOI_data');
    mkdir(outputdir_data);

    for run=1:length(run_valid)
        r=run_valid(run);% The index refer to the real run number of the task (so it might not be consecutive due to exclusive runs)
        subject_mask=strcat(path,'\',subject_list{i},'\run',num2str(r),'\mask.nii');
                     % Subject mask is neccesary for resliced EPI,
                     % because the common mask no longer ensures every sunjects' 
                     % resliced EPI(which has been slightly desplaced) have signals within the mask
        for v=1:size(VOI_name,1)
            clear matlabbatch
            % Make batch for VOI time series extraction (per subject run)
            matlabbatch{1}.spm.util.voi.spmmat = {strcat(path,'\',subject_list{i},'\run',num2str(r),'\SPM.mat')};
            matlabbatch{1}.spm.util.voi.adjust = 1;
            matlabbatch{1}.spm.util.voi.session = 1;
            matlabbatch{1}.spm.util.voi.name = char(strcat(VOI_name{v},'_run',num2str(r)));
            matlabbatch{1}.spm.util.voi.roi{1}.sphere.centre = peak_coordinates{v};
            matlabbatch{1}.spm.util.voi.roi{1}.sphere.radius = radius;
            matlabbatch{1}.spm.util.voi.roi{1}.sphere.move.fixed = 1;
            matlabbatch{1}.spm.util.voi.roi{2}.mask.image = {common_mask};
            matlabbatch{1}.spm.util.voi.roi{2}.mask.threshold = 0.9;
            % For un-resliced common mask only
%             
%             matlabbatch{1}.spm.util.voi.roi{3}.mask.image = {subject_mask}; % Subject mask is neccesary for resliced EPI,
%                                                                             % because the common mask no longer ensures every sunjects' 
%                                                                             % resliced EPI(which has been slightly desplaced) have signals within the mask
%             matlabbatch{1}.spm.util.voi.roi{3}.mask.threshold = 0.9;
            matlabbatch{1}.spm.util.voi.expression = 'i1 & i2';
            save(strcat(outputdir_batch,'\',VOI_name{v},'_run',num2str(r)),'matlabbatch') 
            
            %Run the batch
            spm_jobman('run',matlabbatch);
            sourcefiles=strcat('D:\Yun-Shiuan_LAMBDA\DCM\Part0_First_GLM_for_DCM_single_run_with_reslice\',subject_list{i},'\run',num2str(r));
            movefile(strcat(sourcefiles,'\VOI*'),outputdir_data);
        end
    end
end

%Check if the eigen time series is similar to the averaged time series
%They are actually quite different!!
% mean_time_series=spm_summarise('D:\Yun-Shiuan_LAMBDA\preprocessed_data_without_reslice\df1002\swarun_1.nii',...
%     struct('def','sphere', 'spec',8, 'xyz',[45 -40 47]'),...
%     @mean);
% mean_time_series=mean_time_series(6:end);
% mean_time_series=(mean_time_series-mean(mean_time_series))/std(mean_time_series);

%% Check the unnormalized VOI time series (only for checking)
plot_file='D:\Yun-Shiuan_LAMBDA\DCM\Part1_VOIs_fixed_ROI_for_DCM_with_reslice\Houde_resliced_peak_4mm\summarized_VOI_time_series_unnormalized\';
mkdir(plot_file)
for i=1:size(subject_list,1)
    % Read in the inclusion run indexes and use them to index runs in the following for loop
    run_rows=strcmp(run_inclusion_index.sub_id,subject_list{i});        
    run_valid=table2cell(run_inclusion_index(run_rows,'run_num'));
    run_valid=regexp(run_valid,'(?<=r)\d+','match');
    run_valid=str2double([run_valid{:}]');
    
    outputdir=strcat('D:\Yun-Shiuan_LAMBDA\DCM\Part1_VOIs_fixed_ROI_for_DCM_with_reslice\Houde_resliced_peak_4mm\',subject_list{i});
    outputdir_data=strcat(outputdir,'\VOI_data');
       for run=1:length(run_valid)
            r=run_valid(run);% The index refer to the real run number of the task (so it might not be consecutive due to exclusive runs)
            figure
            std_collect=cell(1,size(VOI_name,1)); % Collect the std of the each time series
            for v=1:size(VOI_name,1)
                VOI_mat=strcat(outputdir_data,'\VOI_',VOI_name{v},'_run',num2str(r),'_1.mat');
                load(VOI_mat)
                plot(Y) %Plot the line
                std_collect{1,v}=std(Y);
                hold on
            end
            legends=strcat(regexprep(cellstr(VOI_name),'_',''),'(',num2str(round(cell2mat(std_collect)',2)),')');
            lgd=legend(legends,'Location','NorthEastOutside');
            title(lgd,'VOI (std)')
            axis([0 120 -4 4])
            plot_file_name=strcat(plot_file,...
                subject_list{i},'_run',num2str(r),'.png')
            saveas(gcf,plot_file_name)
            close(gcf)
       end
end
%% Rescale the extracted VOI time series
for i=1:size(subject_list,1)
    % Read in the inclusion run indexes and use them to index runs in the following for loop
    run_rows=strcmp(run_inclusion_index.sub_id,subject_list{i});        
    run_valid=table2cell(run_inclusion_index(run_rows,'run_num'));
    run_valid=regexp(run_valid,'(?<=r)\d+','match');
    run_valid=str2double([run_valid{:}]');
    
    inputdir=strcat('D:\Yun-Shiuan_LAMBDA\DCM\Part1_VOIs_fixed_ROI_for_DCM_with_reslice\Houde_resliced_peak_4mm\',subject_list{i});
    inputdir_data=strcat(inputdir,'\VOI_data');
    
    outputdir_data=strcat(inputdir,'\VOI_data_normalized');
    mkdir(outputdir_data)
       for run=1:length(run_valid)
            r=run_valid(run);% The index refer to the real run number of the task (so it might not be consecutive due to exclusive runs)
            
            for v=1:size(VOI_name,1)
                VOI_mat=strcat(inputdir_data,'\VOI_',VOI_name{v},'_run',num2str(r),'_1.mat');
                load(VOI_mat)
                xY.u=zscore(xY.u);
                Y=zscore(Y);
                VOI_normalized_mat=strcat(outputdir_data,'\VOI_',VOI_name{v},'_run',num2str(r),'_1.mat');
                save(VOI_normalized_mat,'xY','Y');
            end
       end
end

%% Check the normalized VOI time series (only for checking)
for i=1:size(subject_list,1)
    % Read in the inclusion run indexes and use them to index runs in the following for loop
    run_rows=strcmp(run_inclusion_index.sub_id,subject_list{i});        
    run_valid=table2cell(run_inclusion_index(run_rows,'run_num'));
    run_valid=regexp(run_valid,'(?<=r)\d+','match');
    run_valid=str2double([run_valid{:}]');
    
    outputdir=strcat('D:\Yun-Shiuan_LAMBDA\DCM\Part1_VOIs_fixed_ROI_for_DCM_with_reslice\Houde_resliced_peak_4mm\',subject_list{i});
    outputdir_data=strcat(outputdir,'\VOI_data_normalized');
       for run=1:length(run_valid)
            r=run_valid(run);% The index refer to the real run number of the task (so it might not be consecutive due to exclusive runs)
            figure
            std_collect=cell(1,size(VOI_name,1)); % Collect the std of the each time series
            for v=1:size(VOI_name,1)
                VOI_mat=strcat(outputdir_data,'\VOI_',VOI_name{v},'_run',num2str(r),'_1.mat');
                load(VOI_mat)
                plot(Y) %Plot the line
                std_collect{1,v}=std(Y);
                hold on
            end
            legends=strcat(regexprep(cellstr(VOI_name),'_',''),'(',num2str(round(cell2mat(std_collect)',2)),')');
            lgd=legend(legends,'Location','NorthEastOutside');
            title(lgd,'VOI (std)')
            axis([0 120 -4 4])
            plot_file_name=strcat('D:\Yun-Shiuan_LAMBDA\DCM\Part1_VOIs_fixed_ROI_for_DCM_with_reslice\Houde_resliced_peak_4mm\summarized_VOI_time_series_normalized\',...
                subject_list{i},'_run',num2str(r),'.png')
            saveas(gcf,plot_file_name)
            close(gcf)
       end
end