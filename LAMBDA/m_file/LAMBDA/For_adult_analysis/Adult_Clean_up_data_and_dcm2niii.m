%% Part 0- Clean up & Unzip data& Transfer all dcm. to nii. & Rename EPI images
% For adult. Note that only necessary dcm files are copied from the shared
% folders.
%Enable the usage of helper functions----------
addpath('D:\GoogleDrive\Lambda_code\m_file\LAMBDA\tool_code');

% Constants------------------------------------
PATH_RAW_DATA='D:\Yun-Shiuan_LAMBDA\Adult\raw_data';
FILE_VALID_RUN='D:\Yun-Shiuan_LAMBDA\Adult\Run_inclusion_info\inclusive_runs_indexes.csv';
NUM_RUNS=6;
% Read in run inclusion index info
% and derive subjects with valid runs
run_inclusion_index=read_mixed_csv_to_table(FILE_VALID_RUN);
subject_list=unique(run_inclusion_index.sub_id);

%% EPI: Unzip all tgz to tar + transform dcm to nii
success_list={};
success_mv_list={};
success_dcm2nii_list={};
error_list_sub={};
for id=1:length(subject_list)% Skip df1001 due to the missing T1 image
    %     Should NOT use cd or relative directory in parfor
    %     cd(strcat(path,subject_list{i},'\dicoms'));% Get the file dir. of EPI images
    %     file_interested=cellstr(ls('*.fMRI'));
    %get the list of the interested folder names
    path_this_id=fullfile(PATH_RAW_DATA,subject_list{id});
    path_this_id_linux=...
        regexprep(regexprep(path_this_id,'\','/'),'D:','/cygdrive/d');
    path_this_id_dcm_all=fullfile(PATH_RAW_DATA,subject_list{id},'dicoms');
    
    path_this_id_collect_nii=fullfile(path_this_id,'nii_raw');
    path_this_id_collect_nii_linux=...
        regexprep(regexprep(path_this_id_collect_nii,'\','/'),'D:','/cygdrive/d');
    
    all_files=cellstr(ls(path_this_id_dcm_all));
    
    %The folder names for EPI data
    file_interested = all_files(~cellfun(@isempty,regexp(all_files,'\.fMRI$','match')));
    
    if (numel(file_interested)~=6) % Make sure there is 6 scans for each person
        warning(strcat('Subject ',num2str(id),' does not have 6 runs.'));
        % Check if nii_raw already exist
    elseif (sum(~cellfun(@isempty,regexp(cellstr(ls(fullfile(PATH_RAW_DATA,subject_list{id},'nii_raw'))),'run_\d+.nii','match'))))==NUM_RUNS
        warning(strcat('Subject ',num2str(id),' have already completed dcm2nii.'));
    else
        
        mkdir(path_this_id_collect_nii);
        
        parfor (f=1:NUM_RUNS,4)%Process all EPI images of this subject
            strcat("start- id: ",num2str(id),"; run: ",num2str(f))
            path_this_id_dcm_this_run=fullfile(path_this_id_dcm_all,file_interested{f});
            path_path_this_id_dcm_this_run_linux=...
                regexprep(regexprep(path_this_id_dcm_this_run,'\','/'),'D:','/cygdrive/d');
            index_file_real_run=char(regexp(path_this_id_dcm_this_run,'\d+(?=.fMRI)','match'));
            index_file_real_run=strip(index_file_real_run,'left','0');
            %% For unzipping to dicoms
            %             file_zip=cellstr(ls(strcat(path_interested_subject_run)));
            %             file_tgz=file_zip(~cellfun(@isempty,regexp(file_zip,'tgz$','match')));
            %             file_tgz_full=strcat(path_interested_subject_run,'\',file_tgz);%Full path of the tgz file
            %
            %             % The 'untar' function provided by MATLAB is
            %             % freaking slow and error-prone! (It miss one file, remained it unzipped!!WTF)
            %             % untar(char(file_tgz_full),char(strcat(path_interested_subject_run,'\dicom')));%unzip all dicoms into the folder 'dicom'
            %
            %             %Unzip: Instead, directly Call the bash function ' tar -zxvf {tarball.tgz'
            %             % via bash terminal
            %             %Change directory to linux style(('/cygdrive/m/Yun-Shiuan_LAMBDA/raw_data/df1003/dicoms/00005.fMRI'))
            %             file_tgz_full_linux=regexprep(regexprep(file_tgz_full,'\','/'),'M:','/cygdrive/m');
            %             folder_target_linux=regexprep(regexprep(char(strcat(path_interested_subject_run)),'\','/'),'D:','/cygdrive/d');
            %             %
            %             %                             bash_unzip=char(strcat("C:\cygwin64\bin\bash.exe --login -c ",...
            %             %                                 '"',...
            %             %                                 "cd ",folder_target_linux,"&&",...
            %             %                                 "tar -zxvf ",file_tgz_full_linux,'"'));
            %             %                             dos(bash_unzip);
            %             %                                 %Check if the amount of unzipped slices are 4180
            %             %                                 number_slice=length(cellstr(ls(path_interested_subject_run)))-4;%(minus json, tgz, '.','..'.
            %             %                                 success_list{i,f}=char(strcat('Subject',num2str(i),...
            %             %                                      ' Run',num2str(f),' unzipped successfully,with slices amount:',...
            %             %                                      num2str(number_slice)))
            
            %% Prepare for: Convert dcm to nii (gather them in an independent folder nii_raw)
            bash_pre_dcm2nii=char(strcat("C:\cygwin64\bin\bash.exe --login -c ",...
                '"',...
                "mkdir -p ",path_path_this_id_dcm_this_run_linux,"/dicom_raw","&&",...
                "mv ",path_path_this_id_dcm_this_run_linux,"/*.dcm ",path_path_this_id_dcm_this_run_linux,"/dicom_raw",...
                '"'));
            dos(bash_pre_dcm2nii);
            %% dcm2niix [Not using linux]
            % Create files for dcm2nii
            %(should implement dcm2nii
            % via window cmd code since this is a windoew OS
            % and dcm2niix.exe is also Window-type
            bash_dcm2nii=char(strcat("C:\Users\ychuang26\Desktop\mricrogl\dcm2niix -z n",...
                " -o ",'"',path_this_id_collect_nii,'"'," ",...
                '"',path_this_id_dcm_this_run,"/dicom_raw",'"'));
            dos(bash_dcm2nii);
            %Check if the amount moved is
            %correct
            number_mv=length(cellstr(ls(char(strcat(path_this_id_dcm_this_run,"\dicom_raw")))))-2;%(minus '.','..'.)
            success_mv_list{id,f}=char(strcat('Subject',num2str(id),...
                ' Run',num2str(f),' moved successfully,with slices amount:',...
                num2str(number_mv)));
            %Check if the amount tranformed is
            %correct
            number_dcm2nii=sum(~cellfun(@isempty,...
                regexp(cellstr(ls(path_this_id_collect_nii)),...
                char(strcat('fMRI_',index_file_real_run,'.nii')))));%if there is an nii. which match the run number, then output 1.
            success_dcm2nii_list{id,f}=char(strcat('Subject',num2str(id),...
                ' Run',num2str(f),' dcm2nii successfully,with nii amount:',...
                num2str(number_dcm2nii)));
        end
        
    end
end
cd(PATH_RAW_DATA);
save('EPI_zipping_record.mat','success_list');
save('EPI_moving_record_40_54.mat','success_mv_list');
save('EPI_dcm2nii_record_40_54.mat','success_dcm2nii_list');
save('EPI_error_list.mat','error_list_sub');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Rename the EPI.nii such that it correpsond to real order of task
for id=5%1:length(subject_list)
    
    %List out all old file names along with new names
    path_this_id=fullfile(PATH_RAW_DATA,subject_list{id});
    path_this_id_collect_nii=fullfile(path_this_id,'nii_raw');
    
    all_files=cellstr(ls(strcat(path_this_id,'\nii_raw')));
    file_interested=all_files(~cellfun(@isempty,regexp(all_files,'fMRI_\d+.nii$','match')));
    if(length(file_interested)~=6)
        warning(strcat('Subject ',num2str(id),' does not have 6 runs.'));
    else
        %Create new and unified name for EPI images
        file_interested_index=regexp(file_interested,'(?<=fMRI_)\d+(?=.nii)','match');
        file_interested_table=table(file_interested,...
            str2num(strvcat([file_interested_index{:}]')),...
            'VariableNames',{'filename_nii' 'index_scan'});
        file_interested_table.index_run=...
            file_interested_table.index_scan-...
            min(file_interested_table.index_scan)+1;
        file_interested_table.new_name=...
            cellstr(strcat('run_',num2str(file_interested_table.index_run),'.nii'));
        %Sort table to facilitate monitoring progess
        file_interested_table=sortrows(file_interested_table,'index_run','ascend');
        
        %Renaming one-by-one
        for f=1:6
            oldname=strcat(path_this_id,'\nii_raw\',file_interested_table.filename_nii{f});
            newname=strcat(path_this_id,'\nii_raw\',file_interested_table.new_name{f});
            movefile(oldname,newname);
        end
    end
end

%% T1: Unzip all tgz to tar + transform dcm to nii=======================================================
success_list={};
success_mv_list={};
success_dcm2nii_list={};
error_list_sub={};
for id=1:length(subject_list)% Skip df1001 due to the missing T1 image
    %     Should NOT use cd or relative directory in parfor
    %     cd(strcat(path,subject_list{i},'\dicoms'));% Get the file dir. of EPI images
    %     file_interested=cellstr(ls('*.fMRI'));
    %get the list of the interested folder names
    path_this_id=fullfile(PATH_RAW_DATA,subject_list{id});
    path_this_id_linux=...
        regexprep(regexprep(path_this_id,'\','/'),'D:','/cygdrive/d');
    path_this_id_dcm_all=fullfile(PATH_RAW_DATA,subject_list{id},'dicoms');
    
    path_this_id_collect_nii=fullfile(path_this_id,'nii_raw');
    path_this_id_collect_nii_linux=...
        regexprep(regexprep(path_this_id_collect_nii,'\','/'),'D:','/cygdrive/d');
    
    all_files=cellstr(ls(path_this_id_dcm_all));
    
    %The folder names for EPI data
    file_interested = all_files(~cellfun(@isempty,regexp(all_files,'\.BRAVO','match')));
    for f=1%:numel( file_interested) %Process The T1 image (only f=1)
        
        if (numel(file_interested)~=1) % Make sure there is 1 T1 for each person
            warning(strcat('Subject ',num2str(id),' does not have 1 T1 image.'));
            % Check if nii_raw already exist
        elseif (sum(~cellfun(@isempty,regexp(cellstr(ls(fullfile(PATH_RAW_DATA,subject_list{id},'nii_raw'))),'BRAVO.*.nii','match'))))==1
            warning(strcat('Subject ',num2str(id),' have already completed dcm2nii.'));
        else
            
            mkdir(path_this_id_collect_nii);
            
            %Process the T1 image of this subject
            strcat("start- id: ",num2str(id),"; T1: ",num2str(f))
            path_this_id_dcm_this_run=fullfile(path_this_id_dcm_all,file_interested{f});
            path_path_this_id_dcm_this_run_linux=...
                regexprep(regexprep(path_this_id_dcm_this_run,'\','/'),'D:','/cygdrive/d');
            
            %% For unzipping to dicoms
            if exist([path_this_id_dcm_this_run,'\dicom_raw'],'dir')==7
                warning(strcat('Subject ',num2str(id),' have already zipped.'));
            else
                file_zip=cellstr(ls(strcat(path_this_id_dcm_this_run)));
                file_tgz=file_zip(~cellfun(@isempty,regexp(file_zip,'tgz$','match')));
                file_tgz_full=strcat(path_this_id_dcm_this_run,'\',file_tgz);%Full path of the tgz file
                
                % The 'untar' function provided by MATLAB is
                % freaking slow and error-prone! (It miss one file, remained it unzipped!!WTF)
                % untar(char(file_tgz_full),char(strcat(path_interested_subject_run,'\dicom')));%unzip all dicoms into the folder 'dicom'
                
                %Unzip: Instead, directly Call the bash function ' tar -zxvf {tarball.tgz'
                % via bash terminal
                %Change directory to linux style(('/cygdrive/m/Yun-Shiuan_LAMBDA/raw_data/df1003/dicoms/00005.fMRI'))
                file_tgz_full_linux=regexprep(regexprep(file_tgz_full,'\','/'),'D:','/cygdrive/d');
                path_this_id_dcm_this_run_linux=regexprep(regexprep(char(strcat(path_this_id_dcm_this_run)),'\','/'),'D:','/cygdrive/d');
                % While unzipping,gather them in an independent folder nii_raw.
                bash_unzip=char(strcat("C:\cygwin64\bin\bash.exe --login -c ",...
                    '"',...
                    "cd ",path_this_id_dcm_this_run_linux,"&&",...
                    "mkdir -p ",path_this_id_dcm_this_run_linux,"/dicom_raw","&&",...
                    "tar -zxvf ",file_tgz_full_linux,"&&",...
                    "mv ",path_this_id_dcm_this_run_linux,"/*BRAVO*.dcm ",path_this_id_dcm_this_run_linux,"/dicom_raw",...
                    '"'));
                dos(bash_unzip);
                %Rename the T1 dcms so that it coudl be consider as having
                %valid names by the dcm2niix.exe
                t1_dcm_old_name=cellstr(ls([path_this_id_dcm_this_run,'\dicom_raw']));
                t1_dcm_old_name=t1_dcm_old_name(~ismember(t1_dcm_old_name,{'.','..'}));
                t1_dcm_new_name=regexprep(t1_dcm_old_name,'\.(?!dcm)','_');
                file_t1_dcm_old_name=fullfile([path_this_id_dcm_this_run,'\dicom_raw'],t1_dcm_old_name);
                file_t1_dcm_new_name=fullfile([path_this_id_dcm_this_run,'\dicom_raw'],t1_dcm_new_name);
                for file=1:numel(file_t1_dcm_old_name)
                    movefile(file_t1_dcm_old_name{file},file_t1_dcm_new_name{file});
                end
                
            end
            %Check if the amount of unzipped slices are 176
            number_slice=length(cellstr(ls([path_this_id_dcm_this_run,'\dicom_raw'])))-2;%(minus  '.','..'.
            success_list{id,f}=char(strcat('Subject',num2str(id),...
                ' Run',num2str(f),' unzipped successfully,with slices amount:',...
                num2str(number_slice)))
            
            
            %% dcm2niix [Not using linux]
            % Create files for dcm2nii
            %(should implement dcm2nii
            % via window cmd code since this is a windoew OS
            % and dcm2niix.exe is also Window-type
            bash_dcm2nii=char(strcat("C:\Users\ychuang26\Desktop\mricrogl\dcm2niix -z n",...
                " -o ",'"',path_this_id_collect_nii,'"'," ",...
                '"',path_this_id_dcm_this_run,"\dicom_raw",'"'));
            dos(bash_dcm2nii);
            %Check if the amount moved is
            %correct
            number_mv=length(cellstr(ls(char(strcat(path_this_id_dcm_this_run,"\dicom_raw")))))-2;%(minus '.','..'.)
            success_mv_list{id,f}=char(strcat('Subject',num2str(id),...
                ' Run',num2str(f),' moved successfully,with slices amount:',...
                num2str(number_mv)));
            %Check if the amount tranformed is
            %correct
            number_dcm2nii=sum(~cellfun(@isempty,...
                regexp(cellstr(ls(path_this_id_collect_nii)),...
                char(strcat('BRAVO.*.nii')))));%if there is an nii. which match the run number, then output 1.
            success_dcm2nii_list{id,f}=char(strcat('Subject',num2str(id),...
                ' T1',' dcm2nii successfully,with nii amount:',...
                num2str(number_dcm2nii)));
        end
    end
    
end

cd(PATH_RAW_DATA);
save('zipping_record.mat','success_list');
save('moving_record_40_54.mat','success_mv_list');
save('dcm2nii_record_40_54.mat','success_dcm2nii_list');
save('error_list.mat','error_list_sub');

%% Rename the T1.nii
for id=2:length(subject_list)    
    %List out all old file names along with new names
    path_this_id=fullfile(PATH_RAW_DATA,subject_list{id});
    path_this_id_collect_nii=fullfile(path_this_id,'nii_raw');
    
    all_files=cellstr(ls(strcat(path_this_id,'\nii_raw')));
    file_interested=all_files(~cellfun(@isempty,regexp(all_files,'BRAVO.*.nii$','match')));
    if(length(file_interested)~=1)
        warning(strcat('Subject ',num2str(id),' does not have the T1 image.'));
    else        
        oldname=fullfile(path_this_id,'nii_raw',char(file_interested));
        newname=fullfile(path_this_id,'nii_raw','T1.nii');
        movefile(oldname,newname);
    end
end
