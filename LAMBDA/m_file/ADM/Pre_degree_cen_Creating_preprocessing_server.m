%% Part 1: Create .mat for preprocessing (degree centrality)
path='/bml/Data/Bank6/ADM-YunShiuan/';% Select 3x3x3 (normarlized voxel) version
cd (strcat(path,'pre_mat_degree_cen_3x3x3/'));
load('template.mat');
m= matlabbatch;

%Checking:Only for the CDA data(since manual batch has been set)
%ck=1; %parameter for checking equivalence

% load('pre_23.mat');
% s023= matlabbatch;
% %CDA3= matlabbatch;
%load('.\CDA004\1st(pre).mat');
%CDA4= matlabbatch;
%load('.\CDA005\1st(pre).mat');
%CDA5= matlabbatch;
%manual={CDA3,CDA4,CDA5};% save manual batch(cell) to a cell
%% For loop for subjects
cd (strcat(path,'IMG_nii'));
for id=[23:95],
    %get the T1 file name for further substitution
    % d=ls(strcat('.\S00',num2str(id),'\T1'));% open the T1 folder
    %    dcell=cellstr(d);% convert file names (string array) to cell to enable regex
    %    index= cellfun('length',regexp(dcell,'^s\d')) == 1;%find string(the correct file name) in cell using regex
    %    t1index= dcell(index); %extract the T1 file name
    %get the T2 file name for further substitution
    %  d=ls(strcat('.\CDA00',num2str(id),'\T2'));% open the T2 folder
    %    dcell=cellstr(d);% convert file names (string array) to cell to enable regex
    %    index=cellfun('length',regexp(dcell,'^s\d')) == 1;%find string(the correct file name) in cell using regex
    %    t2index= dcell(index); %extract the T2 file name
  %% For loop for scan number
  for scan=1:5,
    %% Slice timing(m{1,1}.temporal.st)
     %get the EPI 4D file name of each subject for substitution
    %   d=ls(strcat('.\CDA00',num2str(id),'\F',num2str(scan)));% open the scanning session folder
    %    dcell=cellstr(d);% convert file names (string array) to cell to enable regex
    %    index=cellfun('length',regexp(dcell,'^s\d')) == 1;%find string(the correct file name) in cell using regex
    %    scanindex= dcell(index); %extract the file name of the scanning session
    %assign the scans data    
    for i=1:length(m{1,1}.spm.temporal.st.scans{1,1}),
    m{1,1}.spm.temporal.st.scans{1,scan}{i,1}=char(strcat(path,'/IMG_nii','/S0',...
        num2str(id),'/dm',num2str(scan),'.nii',',',num2str(i)));
    end
  end
    %% Realingment(m{2,1}.realign.estimeate)
      %data no need to be asigned, since it is dependent on the previous
      %slice-timed data
    %% Coregister: T2 to EPI (m{3.1}.spatial.coreg)
      %assign source image, i.e., T2
      m{1,3}.spm.spatial.coreg.estimate.source{1,1}=char(strcat(path,'/IMG_nii','/S0',...
        num2str(id),'/T2.nii',',1'));
    %% Coregister: T2,EPI to T1 (m{4.1}.spatial.coreg)
      %assign source image, i.e., T2
        m{1,4}.spm.spatial.coreg.estimate.source{1,1}=char(strcat(path,'/IMG_nii','/S0',...
        num2str(id),'/T2.nii',',1'));
      %assign reference image, i.e., T1
        m{1,4}.spm.spatial.coreg.estimate.ref{1,1}=char(strcat(path,'/IMG_nii','/S0',...
        num2str(id),'/T1.nii',',1'));
      %assign other images, i.e., EPI- dependency
  %% Normalise: Estimate and write (nornalise the T1)
   %Image to Align
    m{1,5}.spm.spatial.normalise.estwrite.subj.vol{1,1}=char(strcat(path,'/IMG_nii','/S0',...
        num2str(id),'/T1.nii',',1'));
    %Image to Write
    m{1,5}.spm.spatial.normalise.estwrite.subj.resample{1,1}=char(strcat(path,'/IMG_nii','/S0',...
        num2str(id),'/T1.nii',',1'));
    %tpm
    m{1,5}.spm.spatial.normalise.estwrite.eoptions.tpm{1,1}= char('/usr/local/spm12/tpm/TPM.nii');
    %affreg
    m{1,5}.spm.spatial.normalise.estwrite.eoptions.affreg=char('eastern');
    %% Normalise: Write(normalise the EPI)
    %all dependency
    %% Smooth
    %dependency
 %% Check equivalence (loop vs. manual)
%  matlabbatch=manual{1,ck};
%  %Check the whole batch
%   check(ck,1)=isequal(m,matlabbatch);
% %  %chech each cell
%   for i=2:length(m)+1    
%  check(ck,i)=isequal(m{1,i-1},matlabbatch{1,i-1});
%   end    
% %  %Slice timing
%   check(ck,i+1)=isequal(m{1,1}.spm.temporal.st.scans{1,scan},...
%   matlabbatch{1,1}.spm.temporal.st.scans{1,scan});
%   % Coregister: T2 to EPI 
%   check(ck,i+2)=isequal(m{1,3}.spm.spatial.coreg.estimate.source,...
%   matlabbatch{1,3}.spm.spatial.coreg.estimate.source);
% %  % Coregister: T2,EPI to T1 
% %    %source image, i.e., T2
%   check(ck,i+3)=isequal(m{1,4}.spm.spatial.coreg.estimate.source,...
%           matlabbatch{1,4}.spm.spatial.coreg.estimate.source);
% %    %reference image, i.e., T1
%   check(ck,i+4)=isequal(m{1,4}.spm.spatial.coreg.estimate.ref,...
%           matlabbatch{1,4}.spm.spatial.coreg.estimate.ref);
% %    % Normalise: Estimate and write (nornalise the T1)
% %        %Image to Align
%   check(ck,i+5)=isequal(m{1,5}.spm.spatial.normalise.estwrite.subj.vol,...
%          matlabbatch{1,5}.spm.spatial.normalise.estwrite.subj.vol);
% %         %Image to Write
%   check(ck,i+6)=isequal(m{1,5}.spm.spatial.normalise.estwrite.subj.vol,...
%          matlabbatch{1,5}.spm.spatial.normalise.estwrite.subj.vol);
% %  ck=ck+1;
 %% Save mat.  
 matlabbatch=m; % SHOULD set the name of batch as "matlabbatch" so that SPM GUI could read it
 save(strcat(path,'/pre_mat_degree_cen_3x3x3/','pre2_',num2str(id)),'matlabbatch');
end

