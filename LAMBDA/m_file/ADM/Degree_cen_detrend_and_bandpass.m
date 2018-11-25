%% Part 2: Degree Centrality: Detrend and Bandpass for Before "Partial out to remain only the Residuals"
addpath('/home/vimchiz/REST_V1.8_130615/') %Use the REST toolbox
path='/bml/Data/Bank6/ADM-YunShiuan/';
cd(strcat(path,'degree_centrality'));

id=csvread(strcat(path,'/young_idlist.csv')); %read only the young's IDs
young_name=strcat('S0',num2str(id));
young_name=cellstr(young_name); %convert char to cell to enable choosing string by index
young_name_40=cellstr(strvcat(young_name{~ismember(young_name,{'S031' 'S040' 'S044'}')})); %Exclude S031,S040,S044

%% Detrend 
for id=1:length(young_name_40);
%Should used unsmoothed preprocessed EPI images ("wradm*.nii")
% Copy each EPI file into the to-be-detrend folder on by one, each person's
% in each folder
path_old = strcat(path,'/IMG_nii/',young_name_40{id},'/degree_cen_normalise_3x3x3/wradm*');
path_new = strcat(path,'degree_centrality/Part_1_detrend_and_bandpass/',young_name_40{id},'/');
% copyfile (path_old,path_new); %Run it only if copying is needed


% rest_detrend, each person's detrended file will be writen under his/her own
% folder
D=dir(path_new);% get the 5 file names
D=struct2cell(D);
filename=D(1,:);
file_index=~cell2mat(D(4,:));%get the indexes of files, exluding folders
filename_valid=filename(file_index);

for f=1:length(filename_valid); %detrend 5 nii. files one by one
rest_detrend(strcat(path_new,filename_valid{f}),'_after_detrend');
end
end

%% Bandpass
for id=1:length(young_name_40);    
    %Get the 5 paths where detrended files locate
    path_detrend = strcat(path,'degree_centrality/Part_1_detrend_and_bandpass/',young_name_40{id},'/');
    D=dir(path_detrend);
    D=struct2cell(D);
    pathname=D(1,:);
    path_invalid=cell2mat(cellfun(@(x) ~(strcmp(x,'.')||strcmp(x,'..')),pathname,'un',0));% to exclude the '.' and '..' folders
    path_index=cell2mat(D(4,:)) & path_invalid ;%get the indexes of folders ,exluding files and those invalid path
    path_valid=pathname(path_index); %Get the 5 paths where detrended files locate
    full_file_index=strcat(path,'degree_centrality/Part_1_detrend_and_bandpass/',young_name_40{id},...
    '/',path_valid,'/detrend_4DVolume.nii');
    full_path_index=strcat(path,'degree_centrality/Part_1_detrend_and_bandpass/',young_name_40{id},...
    '/',path_valid)
    for f=1:5;
        % rest_bandpass
        rest_bandpass(full_file_index{f},2,0.08,0.009,'Yes',0,1);       
        % move the bandpassed file to "after_bandpass folder"
        movefile(strcat(full_path_index{f},'/detrend_4DVolume.nii_filtered/'),...
                 strcat(path,'degree_centrality/Part_1_detrend_and_bandpass/',young_name_40{id},...
                  '/wradm',num2str(f),'.nii_after_bandpass'))
    end
end
