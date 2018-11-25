%% (Part4)First_level:Check Design Matrix for Each person!
%% This Scripts To be ignored! ---> Skip directly to Part5!!

%addpath('/usr/local/spm12/')
%spm_jobman('initcfg') -- before running the batch
%mat = ls(strcat(path,'/pre_mat'))%Failed
path='/home/.bml/Data/Bank6/ADM-YunShiuan';
cd (strcat(path,'/first_design_matrix'));
%% Load Template(S023)

load('./S023/template_checking_design_matrix.mat');
template=SPM;
%% (1)Check if all batches (Design matrix: Prob/Mag/ProbxMag etc.) are identical

%% 1-1 Load Each Young Subject's Design Matrix
id=csvread(strcat(path,'/young_idlist.csv')); %read only the young's IDs
young_name=strcat('S0',num2str(id));
young_name=cellstr(young_name); %convert char to cell to enable choosing string by index
test={};
for id=1:length(young_name);
load(char(strcat('./',young_name(id),'/SPM.mat')));
test{1,id}=SPM;
end
%Check the Design Matrix of all batches(Prob/Mag/ProbxMag) with the template
check_dm_result=cellfun(@(x) isequal(x.xX.X,template.xX.X),test);
find(check_dm_result==0);

%The different ones:'S034'(id=12)'S041'(id=18)

%% 1-2 'S034'->diff found in 12th~15th var of Run 4
%Size of each Design Matrix: (218x5=1090) x 25; 218(208+10xblanks)
I=find(template.xX.X~=test{1,12}.xX.X);
I_c=ceil(I/(size(template.xX.X,1))); %column index
I_r=((I/(size(template.xX.X,1)))-floor(I/(size(template.xX.X,1))))*(size(template.xX.X,1)); %row index
I_run=ceil(round(I_r)./218); % run index
I_rinrun=I_r-(218*(I_run-1));%row index in the specific run
I_value_12=test{1,12}.xX.X(I);
I_value_template=template.xX.X(I);
I_rxc=[I_r,I_c,I_run,I_rinrun,I_value_12,I_value_template,I];

tabulate(I_rxc(:,2)); %'S034'->diff found in 13th~16th var (209x4=836)
tabulate(I_rxc(:,3));%'S034'--> diff found only in run 4th((208+1)*4=836)
%the 10 blanks in run4 of 'S034'&'S041' contain only 9 slides
I_rxc(1:20,:); %Run4 10th slide number differ
test{1,1}.xX.X(655:664,13);sum(ans==0);%10 slides
test{1,12}.xX.X(655:664,13);sum(ans==0); %9 slides (onset error!)
%should I concern?-->Ignore, will be canceled out during second level(summing up 15s)

%% 1-3 'S041'->diff found in 12th~15th var of Run 4
%Size of each Design Matrix: (218x5=1090) x 25; 218(208+10xblanks)
I=find(template.xX.X~=test{1,18}.xX.X);
I_c=ceil(I/(size(template.xX.X,1))); %column index
I_r=((I/(size(template.xX.X,1)))-floor(I/(size(template.xX.X,1))))*(size(template.xX.X,1)); %row index
I_run=ceil(round(I_r)./218); % run index
I_rinrun=I_r-(218*(I_run-1));%row index in the specific run
I_value_18=test{1,18}.xX.X(I);
I_value_template=template.xX.X(I);
I_rxc=[I_r,I_c,I_run,I_rinrun,I_value_18,I_value_template,I];

I_rxc(1:20,:); 
tabulate(I_rxc(:,2)); %'S034'->diff found in 13th~16th var (208x4=832)
tabulate(I_rxc(:,3));%'S034'--> diff found only in run 4th((208)*4=832)
figure();plot(I_rxc(1:208,5)); %var13 of 'S041'
figure();plot(I_rxc(1:208,6)); %var13 of template
figure();plot((I_rxc(1:208,5)-I_rxc(1:208,6))/I_rxc(1:208,6)); %differ up tp 15%
%should I concern?-->Ignore, will be canceled out during second level(summing up 15s)

%Check Onset Time-->Everyone is different
%check_onset_result=cellfun(@(x) isequal(x.Sess(1).U.ons,template.Sess(1).U.ons),test)

%% (2)Check if all batches (Design matrix: Prob/Mag/ProbxMag etc.) are corrected specified
check_dm_prob_result=[];
for i=1:5;
check_dm_prob_result=[check_dm_prob_result,cellfun(@(x) isequal(x.Sess(i).U.P(1).P,template.Sess(i).U.P(1).P),test)];
end %43*5(runs)=215
sum(check_dm_prob_result);

check_dm_mag_result=[];
for i=1:5;
check_dm_mag_result=[check_dm_mag_result,cellfun(@(x) isequal(x.Sess(i).U.P(2).P,template.Sess(i).U.P(2).P),test)];
end %43*5(runs)=215
sum(check_dm_mag_result);



