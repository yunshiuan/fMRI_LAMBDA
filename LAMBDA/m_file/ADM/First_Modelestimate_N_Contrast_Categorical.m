%% (Part5-b)First_level:Create Design Matrix with contrasts and covariates!(Categorical Design)
%% Should Ignore Part3 and Part 4, which did not included correct covariates
%% Tuned To normalise 3x3x3 in the current script
%% NOTE: DID explicit mask in the first level(categorical appraoch)--> Since Second Level will be done in R
%addpath('/usr/local/spm12/')
%spm_jobman('initcfg') -- before running the batch
%mat = ls(strcat(path,'/pre_mat'))%Failed
path='/bml/Data/Bank6/ADM-YunShiuan';
cd (strcat(path,'/First_design_matrix_with_contrastncovariate_categorical/'));

%% Load template
scan=num2cell(1:218);
scan=cellfun(@(x) num2str(x),scan,'un',0);
%% Derive all needed information(Choice Condition;Outcome_accepted_condition;Outcome_rejected_condition;Outcome_null_condition)
id=csvread(strcat(path,'/young_idlist.csv')); %read only the young's IDs
young_name=strcat('S0',num2str(id));
young_name=cellstr(young_name); %convert char to cell to enable choosing string by index

for i=[2:20,22:length(young_name)]; % Should skip i=21(S044) to avoid bug
    load('./template.mat');
    template=matlabbatch;
    %% Read the trials information of each person
        A2={0,0,0,0,0};% save info to "A2"
for j=1:5;
    fid = fopen(char(strcat(path,'/Run_info(csv)/',young_name(i),'/r',num2str(j),'.csv')));  %open file
 %read all contents into data as a char array 
    A = textscan(fid,'%s','Delimiter', ',');
    A2{1,j}=reshape(A{1},8,90)';
end
%% Derive Choice Condition:Onset, Prob., Mag., Current total points

% Ignore those continuous variables in the categorical design, used 15
% categories instead

% prob_1=str2double(A2{1,1}(1:2:length(A2{1,1}),3));
% prob_2=str2double(A2{1,2}(1:2:length(A2{1,1}),3));
% prob_3=str2double(A2{1,3}(1:2:length(A2{1,1}),3));
% prob_4=str2double(A2{1,4}(1:2:length(A2{1,1}),3));
% prob_5=str2double(A2{1,5}(1:2:length(A2{1,1}),3));
% prob_mean=mean(mean([prob_1,prob_2,prob_3,prob_4,prob_5]));
% prob_center=cellfun(@(x) (x-prob_mean),{prob_1,prob_2,prob_3,prob_4,prob_5},'un',0);
% 
% mag_1=str2double(A2{1,1}(1:2:length(A2{1,1}),4));
% mag_2=str2double(A2{1,2}(1:2:length(A2{1,1}),4));
% mag_3=str2double(A2{1,3}(1:2:length(A2{1,1}),4));
% mag_4=str2double(A2{1,4}(1:2:length(A2{1,1}),4));
% mag_5=str2double(A2{1,5}(1:2:length(A2{1,1}),4));
% mag_mean=mean(mean([mag_1,mag_2,mag_3,mag_4,mag_5]));
% mag_center=cellfun(@(x) (x-mag_mean),{mag_1,mag_2,mag_3,mag_4,mag_5},'un',0);
% 
% probxmag=cell2mat(prob_center).*cell2mat(mag_center);
% probxmag=mat2cell(probxmag,45,[1 1 1 1 1]);


% 15 Condition Parameter(HHH=[1 0 0 1 0 0 ...])
cond_list={};
HHH={};HHM={};HHL={};
MHH={};MHM={};MHL={};
MMH={};MMM={};MML={};
MLH={};MLM={};MLL={};
LLH={};LLM={};LLL={};
for s =1:5;
cond=cell2mat(A2{1,s}(1:2:length(A2{1,1}),2));%15 conditions of the run
cond_list{1,s}=cond;%collect 5 runs
HHH{1,s}=(cellfun('length',regexp(cellstr(cond),'(?<=Prob)HH_valueH'))>0);
HHM{1,s}=(cellfun('length',regexp(cellstr(cond),'(?<=Prob)HH_valueM'))>0);
HHL{1,s}=(cellfun('length',regexp(cellstr(cond),'(?<=Prob)HH_valueL'))>0);
MHH{1,s}=(cellfun('length',regexp(cellstr(cond),'(?<=Prob)MH_valueH'))>0);
MHM{1,s}=(cellfun('length',regexp(cellstr(cond),'(?<=Prob)MH_valueM'))>0);
MHL{1,s}=(cellfun('length',regexp(cellstr(cond),'(?<=Prob)MH_valueL'))>0);
MMH{1,s}=(cellfun('length',regexp(cellstr(cond),'(?<=Prob)MM_valueH'))>0);
MMM{1,s}=(cellfun('length',regexp(cellstr(cond),'(?<=Prob)MM_valueM'))>0);
MML{1,s}=(cellfun('length',regexp(cellstr(cond),'(?<=Prob)MM_valueL'))>0);
MLH{1,s}=(cellfun('length',regexp(cellstr(cond),'(?<=Prob)ML_valueH'))>0);
MLM{1,s}=(cellfun('length',regexp(cellstr(cond),'(?<=Prob)ML_valueM'))>0);
MLL{1,s}=(cellfun('length',regexp(cellstr(cond),'(?<=Prob)ML_valueL'))>0);
LLH{1,s}=(cellfun('length',regexp(cellstr(cond),'(?<=Prob)LL_valueH'))>0);
LLM{1,s}=(cellfun('length',regexp(cellstr(cond),'(?<=Prob)LL_valueM'))>0);
LLL{1,s}=(cellfun('length',regexp(cellstr(cond),'(?<=Prob)LL_valueL'))>0);
end
cond_15_value={HHH,HHM,HHL,MHH,MHM,MHL,MMH,MMM,MML,MLH,MLM,MLL,LLH,LLM,LLL};
cond_15_name={'HHH','HHM','HHL','MHH','MHM','MHL','MMH','MMM','MML','MLH','MLM','MLL','LLH','LLM','LLL'};
%Check if correctly specify:[cellstr(cond_15{1,s}) ,cellstr(num2str(cellfun('length',regexp(cellstr(cond_15{1,s}),'(?<=Prob)HH_valueH'))>0))]
%Check if total counts=225 
%sum(cell2mat(cellfun(@(x) sum(cell2mat(x(:))), cond_15_value,'un',0)));

% Ignore those which regard all choice condition as the same
% condition(mudulated by p/m), regard 15 conditions each with unique onset
% time instead
onset_choice_1=str2double(A2{1,1}(1:2:length(A2{1,1}),1));
onset_choice_2=str2double(A2{1,2}(1:2:length(A2{1,1}),1));
onset_choice_3=str2double(A2{1,3}(1:2:length(A2{1,1}),1));
onset_choice_4=str2double(A2{1,4}(1:2:length(A2{1,1}),1));
onset_choice_5=str2double(A2{1,5}(1:2:length(A2{1,1}),1));
onset_choice=cellfun(@(x) (x./1000),{onset_choice_1,onset_choice_2,onset_choice_3,onset_choice_4,onset_choice_5},'un',0);

%Deriving 15 conditions each with unique onset time instead
onset_cond_name=strcat('onset_',cond_15_name);
onset_cond_15={};
for cond=1:15;
tmp=cell2mat(onset_choice).*cell2mat(cond_15_value{cond});
tmp=cellfun(@(x) x(x~=0),num2cell(tmp,1),'un',0);
onset_cond_15{1,cond}=tmp;
end

%total=accumulated points after this trial(the outcome phase)
total_1=str2double(A2{1,1}(2:2:length(A2{1,1}),8));
total_2=str2double(A2{1,2}(2:2:length(A2{1,1}),8));
total_3=str2double(A2{1,3}(2:2:length(A2{1,1}),8));
total_4=str2double(A2{1,4}(2:2:length(A2{1,1}),8));
total_5=str2double(A2{1,5}(2:2:length(A2{1,1}),8));
total={total_1,total_2,total_3,total_4,total_5};

%% Derive Outcome_accept Condition:Onset, Gain, Loss, Total points
%all responses
response_1=str2double(A2{1,1}(1:2:length(A2{1,1}),6));
response_2=str2double(A2{1,2}(1:2:length(A2{1,1}),6));
response_3=str2double(A2{1,3}(1:2:length(A2{1,1}),6));
response_4=str2double(A2{1,4}(1:2:length(A2{1,1}),6));
response_5=str2double(A2{1,5}(1:2:length(A2{1,1}),6));
response={response_1,response_2,response_3,response_4,response_5};
%all outcomes
outcome_1=[str2double(A2{1,1}(2:2:length(A2{1,1}),5)) response_1];
outcome_2=[str2double(A2{1,2}(2:2:length(A2{1,1}),5)) response_2];
outcome_3=[str2double(A2{1,3}(2:2:length(A2{1,1}),5)) response_3];
outcome_4=[str2double(A2{1,4}(2:2:length(A2{1,1}),5)) response_4];
outcome_5=[str2double(A2{1,5}(2:2:length(A2{1,1}),5)) response_5];
outcome={outcome_1,outcome_2,outcome_3,outcome_4,outcome_5};

%index for accpetance
accept_index=cellfun(@(x) (x==4),response,'un',0);
accept_onset_index=cellfun(@(x) logical(reshape(permute([zeros(length(x),1),x],[2,1]),length([x,zeros(length(x),1)])*2,1)),accept_index,'un',0);
%onset under acceptance
onset_outcome_accept={};
for j=1:5;
  onset_outcome_accept{j}=str2double(A2{1,j}(accept_onset_index{j},1));
end
onset_outcome_accept=cellfun(@(x) (x./1000),onset_outcome_accept,'un',0);

%Gain under acceptance
gain=cellfun(@(x) (x(:,2)==4&x(:,1)>0).*x(:,1)+(x(:,2)==3).*0,outcome,'un',0);
for j=1:5;
  gain{j}=gain{j}(accept_index{j},:);
end
%Loss under acceptance
loss=cellfun(@(x) (x(:,2)==4&x(:,1)<0).*-x(:,1)+(x(:,2)==3).*0,outcome,'un',0);
for j=1:5;
  loss{j}=loss{j}(accept_index{j},:);
end
%Total Points under accepttance
total_accept={};
for j=1:5;
tmp=total{j}.*accept_index{j};
total_accept{j}=tmp(tmp~=0);
end
%% Derive Outcome_reject Condition:Onset, Ungain, Unloss, Total points
%index for rejection
reject_index=cellfun(@(x) (x==3),response,'un',0);
reject_onset_index=cellfun(@(x) logical(reshape(permute([zeros(length(x),1),x],[2,1]),length([x,zeros(length(x),1)])*2,1)),reject_index,'un',0);

%onset under rejection
onset_outcome_reject={};
for j=1:5;
  onset_outcome_reject{j}=str2double(A2{1,j}(reject_onset_index{j},1));
end
onset_outcome_reject=cellfun(@(x) (x./1000),onset_outcome_reject,'un',0);

%ungain under rejection
ungain=cellfun(@(x) (x(:,2)==4).*0+(x(:,2)==3&x(:,1)>0).*x(:,1),outcome,'un',0);
for j=1:5;
  ungain{j}=ungain{j}(reject_index{j},:);
end
%unloss under rejection
unloss=cellfun(@(x) (x(:,2)==4).*0+(x(:,2)==3&x(:,1)<0).*x(:,1).*-1,outcome,'un',0);
for j=1:5;
  unloss{j}=unloss{j}(reject_index{j},:);
end
%Total Points under accepttance
total_reject={};
for j=1:5;
tmp=total{j}.*reject_index{j};
total_reject{j}=tmp(tmp~=0);
end
%% Derive Outcome_null Condition:Onset
%index for null
null_index=cellfun(@(x) (isnan(x)),response,'un',0);
null_onset_index=cellfun(@(x) logical(reshape(permute([zeros(length(x),1),x],[2,1]),length([x,zeros(length(x),1)])*2,1)),null_index,'un',0);
%onset under null
onset_outcome_null={};
for j=1:5;
  onset_outcome_null{j}=str2double(A2{1,j}(null_onset_index{j},1));
end
onset_outcome_null=cellfun(@(x) (x./1000),onset_outcome_null,'un',0);

%% Substitute Choice info(Onset/Prob/Mag/Total) to the template
%% Substitute Outcome info(Onset/Gain/Loss/Ungain/Unloss) to the template
%% Substitute multiple regressors(head motions) to the template
%% Specify contrast manager
template{1}.spm.stats.fmri_spec.dir = {char(strcat('/bml/Data/Bank6/ADM-YunShiuan/First_design_matrix_with_contrastncovariate_categorical/',young_name(i),'/normalise_3x3x3'))};
for s=1:5;
% Substitute Choice info(Onset/15Cond/Total) to the template
template{1}.spm.stats.fmri_spec.sess(s).scans = (strcat('/bml/Data/Bank6/ADM-YunShiuan/IMG_nii/',young_name(i),'/normalise_3x3x3/swadm',num2str(s),'.nii,',scan))';
    for cond=1:15; %15 conditions
template{1}.spm.stats.fmri_spec.sess(s).cond(cond).duration = 0;
template{1}.spm.stats.fmri_spec.sess(s).cond(cond).tmod = 0;
template{1}.spm.stats.fmri_spec.sess(s).cond(cond).pmod = struct('name', {}, 'param', {}, 'poly', {});
template{1}.spm.stats.fmri_spec.sess(s).cond(cond).orth = 0;
template{1}.spm.stats.fmri_spec.sess(s).cond(cond).onset = [onset_cond_15{cond}{s}];
template{1}.spm.stats.fmri_spec.sess(s).cond(cond).name = cond_15_name{cond};
    end
% Substitute Outcome_accept info(gain/loss) to the template
template{1}.spm.stats.fmri_spec.sess(s).cond(16).onset= [onset_outcome_accept{1,s}];
template{1}.spm.stats.fmri_spec.sess(s).cond(16).name='accept';
template{1}.spm.stats.fmri_spec.sess(s).cond(16).pmod(1).name = 'gain';
template{1}.spm.stats.fmri_spec.sess(s).cond(16).pmod(1).param = [gain{1,s}];
template{1}.spm.stats.fmri_spec.sess(s).cond(16).pmod(1).poly = 1;
template{1}.spm.stats.fmri_spec.sess(s).cond(16).pmod(2).name = 'loss';
template{1}.spm.stats.fmri_spec.sess(s).cond(16).pmod(2).param = [loss{1,s}];
template{1}.spm.stats.fmri_spec.sess(s).cond(16).pmod(2).poly = 1;
template{1}.spm.stats.fmri_spec.sess(s).cond(16).pmod(3).name = 'p_total_accept';
template{1}.spm.stats.fmri_spec.sess(s).cond(16).pmod(3).param = [total_accept{1,s}];
template{1}.spm.stats.fmri_spec.sess(s).cond(16).pmod(3).poly = 1;


% Substitute Outcome_reject info(ungain/unloss) to the template
template{1}.spm.stats.fmri_spec.sess(s).cond(17).onset= [onset_outcome_reject{1,s}];
template{1}.spm.stats.fmri_spec.sess(s).cond(17).name='reject';
template{1}.spm.stats.fmri_spec.sess(s).cond(17).pmod(1).name = 'ungain';
template{1}.spm.stats.fmri_spec.sess(s).cond(17).pmod(1).param = [ungain{1,s}];
template{1}.spm.stats.fmri_spec.sess(s).cond(17).pmod(1).poly = 1;
template{1}.spm.stats.fmri_spec.sess(s).cond(17).pmod(2).name = 'unloss';
template{1}.spm.stats.fmri_spec.sess(s).cond(17).pmod(2).param = [unloss{1,s}];
template{1}.spm.stats.fmri_spec.sess(s).cond(17).pmod(2).poly = 1;
template{1}.spm.stats.fmri_spec.sess(s).cond(17).pmod(3).name = 'p_total_reject';
template{1}.spm.stats.fmri_spec.sess(s).cond(17).pmod(3).param = [total_reject{1,s}];
template{1}.spm.stats.fmri_spec.sess(s).cond(17).pmod(3).poly = 1;
% Substitute Outcome_null info to the template
if length(onset_outcome_null{s})>0;
    template{1}.spm.stats.fmri_spec.sess(s).cond(18).name='outcome_null';
    template{1}.spm.stats.fmri_spec.sess(s).cond(18).onset= [onset_outcome_null{1,s}];
else
    template{1}.spm.stats.fmri_spec.sess(s).cond(18)=[];
end

% Substitute multiple regressors(head motions) to the template
template{1}.spm.stats.fmri_spec.sess(s).multi_reg{1}=...
char(strcat('/bml/Data/Bank6/ADM-YunShiuan/IMG_nii/',young_name(i),'/normalise_3x3x3/rp_adm',num2str(s),'.txt'));
end %End the loop of sessions
%Explicit mask
template{1}.spm.stats.fmri_spec.mthresh=-Inf;
template{1}.spm.stats.fmri_spec.mask = {'/bml/Data/Bank6/ADM-YunShiuan/smooth_mask_grey_white.nii,1'};

%Contrast manager
    for cond=1:15; %15 conditions
template{3}.spm.stats.con.consess{cond}.tcon.name=cond_15_name{cond};
template{3}.spm.stats.con.consess{cond}.tcon.weights=[repmat(0,1,cond-1),1];% SPM will automatically add the remained 0s
    end
index_end_cond=length(cond_15_name);
template{3}.spm.stats.con.consess{index_end_cond+1}.tcon.name='Accept';
template{3}.spm.stats.con.consess{index_end_cond+1}.tcon.weights=[repmat(0,1,index_end_cond),1];
template{3}.spm.stats.con.consess{index_end_cond+2}.tcon.name='Gain';
template{3}.spm.stats.con.consess{index_end_cond+2}.tcon.weights=[repmat(0,1,index_end_cond+1),1];
template{3}.spm.stats.con.consess{index_end_cond+3}.tcon.name='Loss';
template{3}.spm.stats.con.consess{index_end_cond+3}.tcon.weights=[repmat(0,1,index_end_cond+2),1];
template{3}.spm.stats.con.consess{index_end_cond+4}.tcon.name='Total_accept';
template{3}.spm.stats.con.consess{index_end_cond+4}.tcon.weights=[repmat(0,1,index_end_cond+3),1];

template{3}.spm.stats.con.consess{index_end_cond+5}.tcon.name='Reject';
template{3}.spm.stats.con.consess{index_end_cond+5}.tcon.weights=[repmat(0,1,index_end_cond+4),1];
template{3}.spm.stats.con.consess{index_end_cond+6}.tcon.name='Ungain';
template{3}.spm.stats.con.consess{index_end_cond+6}.tcon.weights=[repmat(0,1,index_end_cond+5),1];
template{3}.spm.stats.con.consess{index_end_cond+7}.tcon.name='Unloss';
template{3}.spm.stats.con.consess{index_end_cond+7}.tcon.weights=[repmat(0,1,index_end_cond+6),1];
template{3}.spm.stats.con.consess{index_end_cond+8}.tcon.name='Total_reject';
template{3}.spm.stats.con.consess{index_end_cond+8}.tcon.weights=[repmat(0,1,index_end_cond+7),1];
%Save .mat
matlabbatch=template;
%Save matlabbatch.mat (the file which could be open by GUI)
save(char(strcat(path,'/First_design_matrix_with_contrastncovariate_categorical/',young_name(i),'/normalise_3x3x3/design_matrix_with_contrast_',num2str(id(i)))),'matlabbatch');
 %Create "SPM.mat" (the design matrix which could be checked)
spm_jobman('run', matlabbatch);
end

%Stop at S044(i=21)--since no ungain in Run1 and 4-> He is excluded due to
%excess of head motion anyways... Just Skip him.

spm_jobman('initcfg')
load('../First_design_matrix_with_contrastncovariate/S023/normalise_3x3x3/SPM.mat')
figure()
plot(SPM.xX.X) % total score increase
plot(SPM.xX.X(1:440,[2,3,4]))%p,m pxm
plot(SPM.xX.X(1:220,[6,7,8]))%accept
plot(SPM.xX.X(1:220,[9,10,11]))%reject
plot(SPM.xX.X(:,[1,18,36,53,70]))%choice baseline
plot(SPM.xX.X(:,[6,24,41,58,75]))%accept baseline
plot(SPM.xX.X(:,[9,27,44,61,78]))%reject baseline

plot(SPM.xX.X(:,[12:17,30:35,47:52,64:69,81:86]))%motion
plot(SPM.xX.X(:,[87:91]))