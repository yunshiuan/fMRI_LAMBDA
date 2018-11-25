%% (Part1) Extract Categorical Betas of contrast(e.g.,Con_0001.nii) for R to do mixed-effects model analysis
%% Note: 
%(1)Should based on Categorical Approach First Level("First_design_matrix_with_contrastncovariate_categorical")
%(2)Select normalized 3x3x3 version(which is more similar to original voxel size, plus it could reduce computation demand)

%addpath('/usr/local/spm12/')
%spm_jobman('initcfg') -- before running the batch
%% Read in all Betas files from Categorical Approach folder
path='/bml/Data/Bank6/ADM-YunShiuan';
cd (strcat(path,'/First_design_matrix_with_contrastncovariate_categorical/'));
%Valid IDs
info=load(strcat(path,'/id_age_gender_sec_hed_sti.mat'));
index=logical(cell2mat(cellfun(@(x) ~(isequal(x,'S031')|isequal(x,'S040')|isequal(x,'S044')),info.a(:,1),'un',0)));%exclude i=9(S031)17(S040)21(S044) who have excessive motions
info.a=info.a(index,:);
young_name=info.a(:,1);

%Create the file names to be saved(cellstr:convert char to cell to enable choosing string by index)
one2nine=strcat(num2str([repmat(0,1,9)]'),num2str([1:9]'));
cond_15_files=cellstr(strcat('con_00',[one2nine;num2str([10:15]')],'.nii'));
collect_all_id_cond={};% the cell to collect all x,y,z,and cond betas result
for f=1:length(cond_15_files);
    for i=1:length(young_name);
    filelocate=char(strcat('./',young_name(i),'/normalise_3x3x3/'));
    V=spm_vol(strcat(filelocate,cond_15_files{f}));%Specify the nii. files interested
    Y=spm_read_vols(V);%Read in the Cond file as matrix(voxel-based coordinate albeit int mni space,
    %yet no need to transder, since one could input the results from R with voxel-based coordiante and save as nii. file. 
    %SPM will automatically convert it to mni spcae)
    c=find(~isnan(Y));;%Only read those which are not NaN and output their linear indexes('find' only select non-zero elements)
    Y_valid=Y(c);%Get the cond value from those valid coordiante
    [x y z]=ind2sub(size(Y),c); %Convert linear indexes to x,y,z indexes
    tmp=[[x,y,z],Y_valid];% Bind coordinate info with Cond.'s Betas
    %[Where could I find the coordinate in SPM Display??e.g., [27,22
    %1]=-0.3798](Although it is consistent with mricron results)
    collect_all_id_cond{i,f}=tmp;
    end
end
%% Write Betas with readable names
cd(strcat(path,'/First_design_matrix_with_contrastncovariate_categorical/15_cond_betas_normalise_3x3x3/'));
cond_15_name={'HHH','HHM','HHL','MHH','MHM','MHL','MMH','MMM','MML','MLH','MLM','MLL','LLH','LLM','LLL'};
for f=1:length(cond_15_files);
    for i=1:length(young_name);
csvwrite(char(strcat('cons_',cond_15_name{f},'_',young_name(i),'.csv')),collect_all_id_cond{i,f});
    end
end