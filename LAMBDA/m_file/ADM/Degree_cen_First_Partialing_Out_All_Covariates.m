%% Part 4: Degree Centrality: Right after Detrend and Bandpas, Get Brain-Relevant Covariates for partialing out in the next step
%Covariates for "Partial out to remain only the Residuals" . They include:
%lateral ventricle activations/ white matter/ head motions, and also their
%first temporal derivatives (one-image-backward difference). (Generated in the previous step)

%In addtion, all task-related covariates should be imported(e.g., prob,
%gain,...), if dealing with task EPI data.(I will deal with it in this script)

%addpath('/usr/local/spm12/')
%spm_jobman('initcfg') -- before running the batch
path='/bml/Data/Bank6/ADM-YunShiuan';
cd(path);
addpath(strcat(path,'/Scripts(m)_files'));% to enable using the function "cor2mni"

id=csvread(strcat(path,'/young_idlist.csv')); %read only the young's IDs
young_name=strcat('S0',num2str(id));
young_name=cellstr(young_name); %convert char to cell to enable choosing string by index
young_name_40=cellstr(strvcat(young_name{~ismember(young_name,{'S031' 'S040' 'S044'}')})); %Exclude S031,S040,S044

%% Retrieve the saved covariates (white matter/CSF/motion/1st derivatives)
cd(strcat(path,'/degree_centrality/Part_2_partial_out_to_remain_only_residuals/covarites/derived_covariate_mat'))
load('collect_6_head_motion.mat');load('collect_6_head_motion_and_white_csf_whole_signal_1_drv.mat');
load('collect_white_csf_whole_signal.mat');

%% Load template
%load('./S023/template_modelestimate_contrastmanager2.mat');
%template=matlabbatch;
scan=num2cell(1:218);
scan=cellfun(@(x) num2str(x),scan,'un',0);

%% Derive all needed information(Choice Condition;Outcome_accepted_condition;Outcome_rejected_condition;Outcome_null_condition)
cd(strcat(path,'/degree_centrality/Part_2_partial_out_to_remain_only_residuals/residual_result_after_partialing'))


for i=[2:length(young_name_40)]; % Only do for the 40 participants
    load('template.mat');
    template=matlabbatch;
    %% Read the trials information of each person
        A2={0,0,0,0,0};% save info to "A2"
for j=1:5;
    fid = fopen(char(strcat(path,'/Run_info(csv)/',young_name_40(i),'/r',num2str(j),'.csv')));  %open file
 %read all contents into data as a char array 
    A = textscan(fid,'%s','Delimiter', ',');
    A2{1,j}=reshape(A{1},8,90)';
end
%% Derive Choice Condition:Onset, Prob., Mag., Current total points

onset_choice_1=str2double(A2{1,1}(1:2:length(A2{1,1}),1));
onset_choice_2=str2double(A2{1,2}(1:2:length(A2{1,1}),1));
onset_choice_3=str2double(A2{1,3}(1:2:length(A2{1,1}),1));
onset_choice_4=str2double(A2{1,4}(1:2:length(A2{1,1}),1));
onset_choice_5=str2double(A2{1,5}(1:2:length(A2{1,1}),1));
onset_choice=cellfun(@(x) (x./1000),{onset_choice_1,onset_choice_2,onset_choice_3,onset_choice_4,onset_choice_5},'un',0);

prob_1=str2double(A2{1,1}(1:2:length(A2{1,1}),3));
prob_2=str2double(A2{1,2}(1:2:length(A2{1,1}),3));
prob_3=str2double(A2{1,3}(1:2:length(A2{1,1}),3));
prob_4=str2double(A2{1,4}(1:2:length(A2{1,1}),3));
prob_5=str2double(A2{1,5}(1:2:length(A2{1,1}),3));
prob_mean=mean(mean([prob_1,prob_2,prob_3,prob_4,prob_5]));
prob_center=cellfun(@(x) (x-prob_mean),{prob_1,prob_2,prob_3,prob_4,prob_5},'un',0);

mag_1=str2double(A2{1,1}(1:2:length(A2{1,1}),4));
mag_2=str2double(A2{1,2}(1:2:length(A2{1,1}),4));
mag_3=str2double(A2{1,3}(1:2:length(A2{1,1}),4));
mag_4=str2double(A2{1,4}(1:2:length(A2{1,1}),4));
mag_5=str2double(A2{1,5}(1:2:length(A2{1,1}),4));
mag_mean=mean(mean([mag_1,mag_2,mag_3,mag_4,mag_5]));
mag_center=cellfun(@(x) (x-mag_mean),{mag_1,mag_2,mag_3,mag_4,mag_5},'un',0);

probxmag=cell2mat(prob_center).*cell2mat(mag_center);
probxmag=mat2cell(probxmag,45,[1 1 1 1 1]);
%total=accumulated - gain/loss in the current trial
total_1=str2double(A2{1,1}(2:2:length(A2{1,1}),8))-str2double(A2{1,1}(2:2:length(A2{1,1}),5));
total_2=str2double(A2{1,2}(2:2:length(A2{1,1}),8))-str2double(A2{1,2}(2:2:length(A2{1,1}),5));
total_3=str2double(A2{1,3}(2:2:length(A2{1,1}),8))-str2double(A2{1,3}(2:2:length(A2{1,1}),5));
total_4=str2double(A2{1,4}(2:2:length(A2{1,1}),8))-str2double(A2{1,4}(2:2:length(A2{1,1}),5));
total_5=str2double(A2{1,5}(2:2:length(A2{1,1}),8))-str2double(A2{1,5}(2:2:length(A2{1,1}),5));
total={total_1,total_2,total_3,total_4,total_5};

%% Derive Outcome_accept Condition:Onset, Gain, Loss, 
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
%% Derive Outcome_reject Condition:Onset, Ungain, Unloss, 
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
result_locate=char(strcat(path,'/degree_centrality/Part_2_partial_out_to_remain_only_residuals/',...
    'residual_result_after_partialing/',young_name_40(i))); %where SPM.mat is going to be written
mkdir(result_locate);%Only use it if needed
template{1}.spm.stats.fmri_spec.dir = {result_locate};
for s=1:5;
% Substitute Choice info(Onset/Prob/Mag/Total) to the template
template{1}.spm.stats.fmri_spec.sess(s).scans = ...
    strcat(path,'/degree_centrality/Part_1_raw_data_and_detrend_and_bandpass/',young_name_40{i},...
           '/wradm',num2str(s),'.nii_after_bandpass/Filtered_4DVolume.nii,',scan)';
template{1}.spm.stats.fmri_spec.sess(s).cond(1).onset = [onset_choice{1,s}];
template{1}.spm.stats.fmri_spec.sess(s).cond(1).pmod(1).name = 'prob';
template{1}.spm.stats.fmri_spec.sess(s).cond(1).pmod(1).param = [prob_center{1,s}];
template{1}.spm.stats.fmri_spec.sess(s).cond(1).pmod(2).name = 'mag';
template{1}.spm.stats.fmri_spec.sess(s).cond(1).pmod(2).param = [mag_center{1,s}];
template{1}.spm.stats.fmri_spec.sess(s).cond(1).pmod(3).name = 'probxmag';
template{1}.spm.stats.fmri_spec.sess(s).cond(1).pmod(3).param = [probxmag{1,s}];
template{1}.spm.stats.fmri_spec.sess(s).cond(1).pmod(4).name = 'p_total';
template{1}.spm.stats.fmri_spec.sess(s).cond(1).pmod(4).param = [total{1,s}];
% Substitute Outcome_accept info(gain/loss) to the template
template{1}.spm.stats.fmri_spec.sess(s).cond(2).onset= [onset_outcome_accept{1,s}];
template{1}.spm.stats.fmri_spec.sess(s).cond(2).pmod(1).name = 'gain';
template{1}.spm.stats.fmri_spec.sess(s).cond(2).pmod(1).param = [gain{1,s}];
template{1}.spm.stats.fmri_spec.sess(s).cond(2).pmod(2).name = 'loss';
template{1}.spm.stats.fmri_spec.sess(s).cond(2).pmod(2).param = [loss{1,s}];
% Substitute Outcome_reject info(ungain/unloss) to the template
template{1}.spm.stats.fmri_spec.sess(s).cond(3).onset= [onset_outcome_reject{1,s}];
template{1}.spm.stats.fmri_spec.sess(s).cond(3).pmod(1).name = 'ungain';
template{1}.spm.stats.fmri_spec.sess(s).cond(3).pmod(1).param = [ungain{1,s}];
template{1}.spm.stats.fmri_spec.sess(s).cond(3).pmod(2).name = 'unloss';
template{1}.spm.stats.fmri_spec.sess(s).cond(3).pmod(2).param = [unloss{1,s}];
% Substitute Outcome_null info to the template
if length(onset_outcome_null{s})>0;
    template{1}.spm.stats.fmri_spec.sess(s).cond(4).name='outcome_null';
    template{1}.spm.stats.fmri_spec.sess(s).cond(4).onset= [onset_outcome_null{1,s}];
else
    template{1}.spm.stats.fmri_spec.sess(s).cond(4)=[];
end
%% Add Regressors(White Matter/Ventricle CSF/1st div of White,CSF,and 6 head motions)
template{1}.spm.stats.fmri_spec.sess(s).regress(1).name='white';
template{1}.spm.stats.fmri_spec.sess(s).regress(1).val=collect_white_deep_signal{i,s};
template{1}.spm.stats.fmri_spec.sess(s).regress(2).name='ventricle';
template{1}.spm.stats.fmri_spec.sess(s).regress(2).val=collect_csf_deep_signal{i,s};
template{1}.spm.stats.fmri_spec.sess(s).regress(3).name='whole';
template{1}.spm.stats.fmri_spec.sess(s).regress(3).val=collect_whole_signal{i,s};
template{1}.spm.stats.fmri_spec.sess(s).regress(4).name='white_div';
template{1}.spm.stats.fmri_spec.sess(s).regress(4).val=collect_white_deep_signal_1_drv{i,s};
template{1}.spm.stats.fmri_spec.sess(s).regress(5).name='ventricle_div';
template{1}.spm.stats.fmri_spec.sess(s).regress(5).val=collect_csf_deep_signal_1_drv{i,s};
template{1}.spm.stats.fmri_spec.sess(s).regress(6).name='whole_div';
template{1}.spm.stats.fmri_spec.sess(s).regress(6).val=collect_whole_signal_1_drv{i,s};


template{1}.spm.stats.fmri_spec.sess(s).regress(7).name='M1_div';
template{1}.spm.stats.fmri_spec.sess(s).regress(7).val=collect_6_motion_1_drv{i,s}{1};
template{1}.spm.stats.fmri_spec.sess(s).regress(8).name='M2_div';
template{1}.spm.stats.fmri_spec.sess(s).regress(8).val=collect_6_motion_1_drv{i,s}{2};
template{1}.spm.stats.fmri_spec.sess(s).regress(9).name='M3_div';
template{1}.spm.stats.fmri_spec.sess(s).regress(9).val=collect_6_motion_1_drv{i,s}{3};
template{1}.spm.stats.fmri_spec.sess(s).regress(10).name='M4_div';
template{1}.spm.stats.fmri_spec.sess(s).regress(10).val=collect_6_motion_1_drv{i,s}{4};
template{1}.spm.stats.fmri_spec.sess(s).regress(11).name='M5_div';
template{1}.spm.stats.fmri_spec.sess(s).regress(11).val=collect_6_motion_1_drv{i,s}{5};
template{1}.spm.stats.fmri_spec.sess(s).regress(12).name='M6_div';
template{1}.spm.stats.fmri_spec.sess(s).regress(12).val=collect_6_motion_1_drv{i,s}{6};

% Substitute multiple regressors(6 head motions) to the template
template{1}.spm.stats.fmri_spec.sess(s).multi_reg{1}=...
    char(strcat(path,'/degree_centrality/Part_2_partial_out_to_remain_only_residuals/covarites/head_motion_raw/',...
    young_name_40(i),'/rp_adm',num2str(s),'.txt'));
end

template{1}.spm.stats.fmri_spec.mthresh=-Inf;
template{2}.spm.stats.fmri_est.write_residuals = 1;
%Contrast manager
template{3}.spm.stats.con.consess{1}.tcon.name='Prob';
template{3}.spm.stats.con.consess{1}.tcon.weights=[0,1];% SPM will automatically add the remained 0s
template{3}.spm.stats.con.consess{2}.tcon.name='Mag';
template{3}.spm.stats.con.consess{2}.tcon.weights=[repmat(0,1,2),1];
template{3}.spm.stats.con.consess{3}.tcon.name='ProbxMag';
template{3}.spm.stats.con.consess{3}.tcon.weights=[repmat(0,1,3),1];
template{3}.spm.stats.con.consess{4}.tcon.name='Accept';
template{3}.spm.stats.con.consess{4}.tcon.weights=[repmat(0,1,5),1];
template{3}.spm.stats.con.consess{5}.tcon.name='Gain';
template{3}.spm.stats.con.consess{5}.tcon.weights=[repmat(0,1,6),1];
template{3}.spm.stats.con.consess{6}.tcon.name='Loss';
template{3}.spm.stats.con.consess{6}.tcon.weights=[repmat(0,1,7),1];
template{3}.spm.stats.con.consess{7}.tcon.name='Reject';
template{3}.spm.stats.con.consess{7}.tcon.weights=[repmat(0,1,8),1];
template{3}.spm.stats.con.consess{8}.tcon.name='Ungain';
template{3}.spm.stats.con.consess{8}.tcon.weights=[repmat(0,1,9),1];
template{3}.spm.stats.con.consess{9}.tcon.name='Unloss';
template{3}.spm.stats.con.consess{9}.tcon.weights=[repmat(0,1,10),1];
template{3}.spm.stats.con.consess{10}.tcon.name='Choice';
template{3}.spm.stats.con.consess{10}.tcon.weights=[1];
%Save .mat
matlabbatch=template;
%Save matlabbatch.mat (the file which could be open by GUI)
save(char(strcat(path,'/degree_centrality/Part_2_partial_out_to_remain_only_residuals/residual_result_after_partialing/',...
    young_name_40(i),'/design_matrix_with_contrast_',young_name_40(i))),'matlabbatch');
 %Create "SPM.mat" (the design matrix which could be checked)
spm_jobman('run', matlabbatch);
end

