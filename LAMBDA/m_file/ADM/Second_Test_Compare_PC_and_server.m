path='/bml/Data/Bank6/ADM-YunShiuan';
cd (strcat(path,'/Second_ROI_defined_by_value2prob'));
Server=load('./prob/SPM.mat');
%% Server:The SPM.mat result of Prob.(With SVS)
%addpath('/usr/local/spm12/')
%spm_jobman('initcfg') -- before running the batch
path='/bml/Data/Bank6/ADM-YunShiuan';

%% My_PC:The SPM.mat result of Prob.(With SVS)
PC=load('./Test_compare_pc_and_server/My_PC/test_prob_withSVS/SPM.mat');
%Check Design Matrix
isequal(Server.SPM.xX.X,PC.SPM.xX.X)% Design Matrix Unequal-->Tiny diff at Security

%Check raw regressors' data (Note: SVS are from scaled value by 43n)
isequal(Server.SPM.xC(1),PC.SPM.xC(1))%Gender equal
isequal(Server.SPM.xC(2),PC.SPM.xC(2))%Hedonism equal

isequal(Server.SPM.xC(3),PC.SPM.xC(3))%Security unequal
isequal(Server.SPM.xC(3).rc,PC.SPM.xC(3).rc)%Raw data equal
isequal(Server.SPM.xC(3).c,PC.SPM.xC(3).c)% Centered data unqual-->Tiny diff at Security

%Other differences
isequal(PC.SPM.xVol.FWHM,Server.SPM.xVol.FWHM) %WTF is FWHM?
isequal(PC.SPM.xVol.R,Server.SPM.xVol.R) %WTF is R?
%SPM.xCon is different--PC.SPM.xCon.X0, PC.SPM.xCon.X1o, and eidf are different
isequal(PC.SPM.xCon,Server.SPM.xCon) 
pc_xcon=struct2cell(PC.SPM.xCon);
server_xcon=struct2cell(Server.SPM.xCon);
a=0;
for i=1:size(pc_xcon,1);
   for j=1:size(pc_xcon,2);
       for k=1:size(pc_xcon,3);
           a(i,j,k)=isequal(pc_xcon{i,j,k},server_xcon{i,j,k});
       end;
   end;
end;
isequal(PC.SPM.SPMid,Server.SPM.SPMid) %Version Different: Server:"v6015"; PC:"v6842"
%Recreate the PC-SPM.mat by Server----------------------------------




