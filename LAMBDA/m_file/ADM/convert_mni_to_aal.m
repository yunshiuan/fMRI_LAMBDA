addpath('C:\Program Files\MATLAB\R2013a\spm12')
%addpath(strcat(path,'/Scripts(m)_files')) %To enable the usage of cor2mni
%function
cd('D:\desktop\googledrive\Brainimage\Scripts(m)_files')
V=spm_vol('aal.nii');%Specify the "aal.nii" file which is in "MRIcron"
Y=spm_read_vols(V);%size of Y: 181 217 181 (voxel-based coordinate)
%Convert the aal.nii to mni coordinate
T=V.mat; %abtain the transformation rule
    %cor2mni([94,153,66],T) %Convert to mni!
Y_mni=[];
for i=1:length(unique(Y));
c=find(Y==i);
[x y z]=ind2sub(size(Y),c);% find out the x,y,z cordinate values for the given indexes
tmp=[cor2mni([x,y,z],T),repmat(i,[size(c,1),1])];%convert to mni coordinate and then bind the column of aal index(i)
Y_mni=[Y_mni;tmp];
end
csvwrite('mni_to_aal.csv',Y_mni)