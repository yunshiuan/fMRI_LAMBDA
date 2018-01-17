%% Part 0- Clean up & Unzip data& Transfer all dcm. to nii. & Rename EPI images
addpath('D:\GoogleDrive\Lambda_code\m_file\LAMBDA\tool_code');
path='M:\Yun-Shiuan_LAMBDA\raw_data\';
cd(path);
subject_list=cellstr(ls('*')); % subject list
subject_list=subject_list(~cellfun(@isempty,regexp(subject_list,'(?<=df)\d+','match')));% include only subject name

%% Unzip all tgz to tar + transform dcm to nii 
success_list={};
success_mv_list={};
success_dcm2nii_list={};
error_list_sub={};
for i=2:length(subject_list)% Skip df1001 due to the missing T1 image
%     Should NOT use cd or relative directory in parfor
%     cd(strcat(path,subject_list{i},'\dicoms'));% Get the file dir. of EPI images
%     file_interested=cellstr(ls('*.fMRI'));
      %get the list of the interested folder names
      Allfiles=cellstr(ls(strcat(path,subject_list{i},'\dicoms')));
      file_interested =Allfiles(~cellfun(@isempty,regexp(Allfiles,'\.fMRI$','match')));
             if (length(file_interested)~=6) % Make sure there is 6 scans for each person
                warning(strcat('Subject ',num2str(i),' does not have 6 runs.'));
             % Check if nii_raw already exist
             elseif (sum(~cellfun(@isempty,regexp(cellstr(ls(strcat(path,subject_list{i},'\nii_raw'))),'fMRI.*.nii','match'))))==6
                warning(strcat('Subject ',num2str(i),' have already completed dcm2nii.'));
             else
                parfor (f=1:5,4)%Unzip all EPI images of this subject
                    [i,f]
                        path_interested_subject_run=...
                                strcat(path,subject_list{i},'\dicoms\',...
                                    file_interested{f});
                        file_zip=cellstr(ls(strcat(path_interested_subject_run)));
                        file_tgz=file_zip(~cellfun(@isempty,regexp(file_zip,'tgz$','match')));
                        file_tgz_full=strcat(path_interested_subject_run,'\',file_tgz);%Full path of the tgz file
                       
                                        % The 'untar' function provided by MATLAB is
                                        % freaking slow and error-prone! (It miss one file, remained it unzipped!!WTF)
                                        % untar(char(file_tgz_full),char(strcat(path_interested_subject_run,'\dicom')));%unzip all dicoms into the folder 'dicom'
                       
                        %Unzip: Instead, directly Call the bash function ' tar -zxvf {tarball.tgz'
                            % via bash terminal
                            %Change directory to linux style(('/cygdrive/m/Yun-Shiuan_LAMBDA/raw_data/df1003/dicoms/00005.fMRI'))
                            file_tgz_full_linux=regexprep(regexprep(file_tgz_full,'\','/'),'M:','/cygdrive/m');
                            folder_target_linux=regexprep(regexprep(char(strcat(path_interested_subject_run)),'\','/'),'M:','/cygdrive/m');
% 
%                             bash_unzip=char(strcat("C:\cygwin64\bin\bash.exe --login -c ",...
%                                 '"',...
%                                 "cd ",folder_target_linux,"&&",...
%                                 "tar -zxvf ",file_tgz_full_linux,'"'));
%                             dos(bash_unzip);
%                                 %Check if the amount of unzipped slices are 4180
%                                 number_slice=length(cellstr(ls(path_interested_subject_run)))-4;%(minus json, tgz, '.','..'.
%                                 success_list{i,f}=char(strcat('Subject',num2str(i),...
%                                      ' Run',num2str(f),' unzipped successfully,with slices amount:',...
%                                      num2str(number_slice)))
                         
                        %Prepare for: Convert dcm to nii (gather them in an independent folder nii_raw)
                            path_interested_subject=...
                                    strcat(path,subject_list{i});
                            path_interested_subject_linux=regexprep(regexprep(path_interested_subject,'\','/'),'M:','/cygdrive/m');
                            % Create files for dcm2nii
                                %(should implement dcm2vii
                                % via window cmd code since this is a windoew OS
                                % and dcm2niix.exe is also Window-type
                                bash_pre_dcm2nii=char(strcat("C:\cygwin64\bin\bash.exe --login -c ",...
                                    '"',...                                                    
                                    "mkdir -p ",folder_target_linux,"/dicom_raw","&&",...
                                    "mv ",folder_target_linux,"/*.dcm ",folder_target_linux,"/dicom_raw","&&",...
                                    "mkdir -p ",strcat(path_interested_subject_linux,'/nii_raw'),...% dcm2nii, gather the created nii into 'nii_raw' folder
                                    '"'));
                                dos(bash_pre_dcm2nii);
                        %dcm2niix
                              bash_dcm2nii=char(strcat("C:\Users\ychuang26\Desktop\mricrogl\dcm2niix -z n",...
                                  " -o ",'"',path_interested_subject,"\nii_raw",'"'," ",...
                                  '"',path_interested_subject_run,"\dicom_raw",'"'));
                              dos(bash_dcm2nii);
                               %Check if the amount moved is
                               %correct
                                number_mv=length(cellstr(ls(char(strcat(path_interested_subject_run,"\dicom_raw")))))-2;%(minus '.','..'.)
                                success_mv_list{i,f}=char(strcat('Subject',num2str(i),...
                                     ' Run',num2str(f),' moved successfully,with slices amount:',...
                                     num2str(number_mv)));
                              %Check if the amount tranformed is
                               %correct
                                file_index=regexp(file_interested{f},'(?<=000)\d+(?=\.fMRI)','match');
                                number_dcm2nii=sum(~cellfun(@isempty,...
                                    regexp(cellstr(ls(char(strcat(path_interested_subject,"\nii_raw")))),...
                                    char(strcat('fMRI_',file_index,'.nii')))));%if there is an nii. which match the run number, then output 1.
                                success_dcm2nii_list{i,f}=char(strcat('Subject',num2str(i),...
                                     ' Run',num2str(f),' dcm2nii successfully,with nii amount:',...
                                     num2str(number_dcm2nii)));
                   end

                end
end
cd(path);
save('zipping_record.mat','success_list');
save('moving_record_40_54.mat','success_mv_list');
save('dcm2nii_record_40_54.mat','success_dcm2nii_list');
save('error_list.mat','error_list_sub');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Rename the EPI.nii such that it correpsond to real order of task
%(note: df1014 & de1018 has abnormal amount of runs)
for i=2:length(subject_list)% Skip df1001 due to the missing T1 image
        
    %List out all old file names along with new names
          Allfiles=cellstr(ls(strcat(path,subject_list{i},'\nii_raw')));
          file_interested=Allfiles(~cellfun(@isempty,regexp(Allfiles,'fMRI_\d+.nii$','match')));
      if(length(file_interested)~=6)
              warning(strcat('Subject ',num2str(i),' does not have 6 runs.'));
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
                oldname=strcat(path,subject_list{i},'\nii_raw\',file_interested_table.filename_nii{f});
                newname=strcat(path,subject_list{i},'\nii_raw\',file_interested_table.new_name{f});
                movefile(oldname,newname);
            end           
      end
end
%% Move T1 image to the right place & Rename
success_T1_mv_list={};
error_list_T1={};
for i=2:length(subject_list)% Skip df1001 due to the missing T1 image
        try
    %List out all old file names along with new names
          Allfiles=cellstr(ls(strcat(path,subject_list{i},'\niftis\*MPnRAGE*\')));
          file_interested =Allfiles(~cellfun(@isempty,regexp(Allfiles,'_mc.nii$','match')));
          %Show the parent path to facilitate "movefile" later on
          file_intetested_parent_path=dir(strcat(path,subject_list{i},'\niftis\*MPnRAGE*\'));
          file_intetested_parent_path=file_intetested_parent_path(1).folder;
      if(length(file_interested)~=1)
              w=warning(strcat('Subject ',num2str(i),' does not have 1 MPnRAGE.'));
      else
          old_path=char(strcat(file_intetested_parent_path,'\',file_interested));
          new_path=char(strcat(path,subject_list{i},'\nii_raw\T1.nii'));
          success=movefile(old_path,new_path);
          
      end
        catch ME
            error_list_T1{i,1}=char(strcat('error:subject',num2str(i),ME.message))
        end
            success_T1_mv_list{i,1}=char(strcat('Subject',num2str(i),...
                                     ' T1 moved with status:',num2str(success))); 
         [i]
end
save('success_T1_mv.mat','success_T1_mv_list');

210361
