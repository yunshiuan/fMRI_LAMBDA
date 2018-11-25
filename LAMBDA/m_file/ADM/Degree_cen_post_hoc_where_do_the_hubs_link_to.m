%% Part 8: Degree Centrality: Find out where do the hubs where "(DC~SVS*+sex)" link to
%Select the without cerebelum mode
% Two Mode to Choose- Standardized(Z_matrix) or Unstandardized(D_matrix)
%addpath('/usr/local/spm12/')
%spm_jobman('initcfg') -- before running the batch
%addpath('/usr/local/spm12/')
%spm_jobman('initcfg') -- before running the batch
path='/bml/Data/Bank6/ADM-YunShiuan';
cd(path);
addpath(strcat(path,'/Scripts(m)_files'));% to enable using the function "cor2mni"

id=csvread(strcat(path,'/young_idlist.csv')); %read only the young's IDs
young_name=strcat('S0',num2str(id));
young_name=cellstr(young_name); %convert char to cell to enable choosing string by index
young_name_40=cellstr(strvcat(young_name{~ismember(young_name,{'S031' 'S040' 'S044'}')})); %Exclude S031,S040,S044
%% Load in the hubs center location and relevant data
load('/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Dosenbach_Science_160ROIs_Center.mat',...
'Dosenbach_Science_160ROIs_Center');
ROI_160_center=Dosenbach_Science_160ROIs_Center;

cd('/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_4_Regression_degree_cen_by_SVS_by_R/lm_result_unstandardized/without_cerebelum');
load('collect_lm_160_ROI_by_2_SVS.mat'); % Load in (DC~SVS*) hubs
load('overall_t_p_160_ROI.mat'); % Load in (overall*) hubs
load('cerebelum_aal_name.mat'); % Load in cerebelum ROI name to exclude them
aal_valid_index=~ismember(aal_name,cerebelum_aal_name);%142 of the ROIs are not cerebelum regions
aal_valid_name=aal_name(aal_valid_index);% Those 142 ROI aal names

ROI_valid_center=ROI_160_center(aal_valid_index,:); %Those 142 ROI center location
%Should transpose to 142x1 first
overall_t_map=overall_t_map';
overall_p_map=overall_p_map';

%% Derive the (DC~SVS*) & (overall*) hubs locations
%% Security %  2 of ROI: 'Supp_Motor_Area_L''Supp_Motor_Area_R'
overall_and_security_hub_index = (security_p<0.1)&(overall_p_map<0.1)&...
                                 (security_t>0)&(overall_t_map>0);
overall_and_security_hub_name = aal_valid_name(overall_and_security_hub_index);
overall_and_security_hub_center=ROI_valid_center(overall_and_security_hub_index,:);
%
%% Hedonism %  4 of ROI: 'Frontal_Sup_L''Postcentral_L''Parietal_Inf_L''Angular_R'
% use one-tail p<0.05 criteria (Sign should be positive)
overall_and_hedonism_hub_index = (hedonism_p<0.1)&(overall_p_map<0.1)&...
                                 (hedonism_t>0)&(overall_t_map>0);%  7 of ROI
overall_and_hedonism_hub_name = aal_valid_name(overall_and_hedonism_hub_index);
overall_and_hedonism_hub_center=ROI_valid_center(overall_and_hedonism_hub_index,:);
%
%% Securiy & Hedonism % None of ROI
% use one-tail p<0.05 criteria (Sign should be positive)
overall_and_security_hedonism_hub_index=overall_and_security_hub_index & overall_and_hedonism_hub_index;

%% Derive the (DC~SVS*) & (overall n.s.) hubs locations
%% Security only (overall n.s.) % 9 of ROI: 
%% 'Cingulum_Mid_R''Precentral_R''Supp_Motor_Area_R''Insula_R''Putamen_R''Insula_L''Putamen_R'
%% 'Putamen_L''Temporal_Mid_R'

% use one-tail p<0.05 criteria (Sign should be positive); 
% overall_t shouldbe n.s. or negative
only_security_hub_index = (security_p<0.1)&(overall_p_map>0.1|overall_t_map<0)&...
                                 (security_t>0);
only_security_hub_name = aal_valid_name(only_security_hub_index);
only_security_hub_center=ROI_valid_center(only_security_hub_index,:);
%
%% Hedonism only (overall n.s.)  % 9 of ROI:'Frontal_Inf_Oper_R''Insula_R''Putamen_R''Rolandic_Oper_L''Insula_R''Insula_L'
%% 'Rolandic_Oper_R' 'Temporal_Sup_R''Temporal_Mid_L'
% use one-tail p<0.05 criteria (Sign should be positive); 
% overall_t shouldbe n.s. or negative
only_hedonism_hub_index = (hedonism_p<0.1)&(overall_p_map>0.1|overall_t_map<0)&...
                                 (hedonism_t>0);%  7 of ROI
only_hedonism_hub_name = aal_valid_name(only_hedonism_hub_index);
only_hedonism_hub_center=ROI_valid_center(only_hedonism_hub_index,:);
tidy_only_hedonism_hub=table(hedonism_t(only_hedonism_hub_index),hedonism_p(only_hedonism_hub_index),...
    only_hedonism_hub_name,only_hedonism_hub_center);
tidy_only_hedonism_hub.Properties.VariableNames={'t','p_unc','aal_name','mni'};
tidy_only_hedonism_hub=sortrows(tidy_only_hedonism_hub,'p_unc','ascend');
%
%% Securiy & Hedonism % 2 of ROI 'Insula_R''Putamen_R''Insula_L'
% use one-tail p<0.05 criteria (Sign should be positive)
only_security_hedonism_hub_index=only_security_hub_index & only_hedonism_hub_index;
only_security_hedonism_hub_name = aal_valid_name(only_security_hedonism_hub_index);
only_security_hedonism_hub_center=ROI_valid_center(only_security_hedonism_hub_index,:);
%

%% Logic of FDR(Storey)
% (Remained to be read)http://www.nonlinear.com/support/progenesis/comet/faq/v2.0/pq-values.aspx
% https://www.bioconductor.org/packages/devel/bioc/vignettes/qvalue/inst/doc/qvalue.pdf
%% (None pass FDR correction)FDR & Reduce Comparisons(Focus only on overall significant positive hubs)
overall_p_positive_map_index=overall_t_map>0;%Only focus on 81 positive hubs
overall_p_positive_map=overall_p_map(overall_p_positive_map_index);

overall_p_positive_map_fdr=mafdr(overall_p_positive_map,'BHFDR',0);
% focus on 64 overall FDR significant positive hubs
overall_p_positive_fdr_index=logical([]);
order=1;
for i=1:size(overall_p_positive_map_index,1)
if overall_p_positive_map_index(i)==1
    overall_p_positive_fdr_index(i,1)=(overall_p_positive_map_fdr(order)<0.05);
    order=order+1;
else
    overall_p_positive_fdr_index(i,1)=logical(0);
end
end
% overall_p_positive_fdr_index=(find(overall_p_positive_map_index).*(overall_p_positive_map_fdr<0.05))>0; 
overall_p_positive_fdr_name=aal_valid_name(overall_p_positive_fdr_index);
tidy_overall_p_positive_fdr=...
    table(round(overall_t_map(overall_p_positive_fdr_index),2),...
    round(overall_p_map(overall_p_positive_fdr_index),3),...
    round(overall_p_positive_map_fdr(overall_p_positive_map_fdr<0.05),3),...
    overall_p_positive_fdr_name,ROI_valid_center(overall_p_positive_fdr_index,:),...
    round(security_t(overall_p_positive_fdr_index),2),...
    round(security_p(overall_p_positive_fdr_index),2),...
    round(hedonism_t(overall_p_positive_fdr_index),2),...
    round(hedonism_p(overall_p_positive_fdr_index),2));
tidy_overall_p_positive_fdr.Properties.VariableNames=...
    {'t','p_unc','p_fdr','aal_name','mni','sec_t','sec_p_unc','hed_t','hed_p_unc'};
tidy_overall_p_positive_fdr=sortrows(tidy_overall_p_positive_fdr,'t','descend');
% DC~SVS (Only focus on 67 of overall positive significant hubs)
% Security
security_p_in_positive_overall_hub=security_p(overall_p_positive_fdr_index);
security_p_in_positive_overall_hub_fdr=mafdr(security_p_in_positive_overall_hub,'BHFDR',0);
% Hedonism
hedonism_p_in_positive_overall_hub=hedonism_p(overall_p_positive_fdr_index);
hedonism_p_in_positive_overall_hub_fdr=mafdr(hedonism_p_in_positive_overall_hub,'BHFDR',0);
%% (Unused)FDR of overall SVS effects
% overall value effect (All ns. QAQ)
% security_p_fdr=mafdr(security_p,'BHFDR',0); %Use Storey (2002) instead of 'BH'
% hedonism_p_fdr=mafdr(hedonism_p,'BHFDR',0); %Use Storey (2002) instead of 'BH'
security_p_fdr=mafdr(security_p,'BHFDR',1); %Use 'BH'
hedonism_p_fdr=mafdr(hedonism_p,'BHFDR',1); %Use 'BH'

tidy_security_p_fdr=table(...
    round(security_t,2),...
    round(security_p,2),...
    round(security_p_fdr,2),...
    aal_valid_name,...
    ROI_valid_center,...
    round(overall_t_map,2),...
    round(overall_p_map,2));
tidy_security_p_fdr.Properties.VariableNames={'t','p_unc','p_fdr','aal_name','mni','overall_t','overall_p_unc'};
tidy_security_p_fdr=sortrows(tidy_security_p_fdr,'p_unc','ascend');

tidy_hedonism_p_fdr=table(...
    round(hedonism_t,2),...
    round(hedonism_p,2),...
    round(hedonism_p_fdr,2),...
    aal_valid_name,...
    ROI_valid_center,...
    round(overall_t_map,2),...
    round(overall_p_map,2));
tidy_hedonism_p_fdr.Properties.VariableNames={'t','p_unc','p_fdr','aal_name','mni','overall_t','overall_p_unc'};
tidy_hedonism_p_fdr=sortrows(tidy_hedonism_p_fdr,'p_unc','ascend');

% One tail: FDR is not suitable for one-tail p value
%Positive tail
% security_p_one_tail=(security_t>0).*(security_p/2)+(security_t<0).*(1-(security_p/2));
% security_p_fdr_one_tail=mafdr(security_p_one_tail,'BHFDR',0);
% hedonism_p_one_tail=(hedonism_t>0).*(hedonism_p/2)+(hedonism_t<0).*(1-(hedonism_p/2));
% hedonism_p_fdr_one_tail=mafdr(hedonism_p_one_tail,'BHFDR',0);
%Negative tail
% security_p_one_tail=(security_t<0).*(security_p/2)+(security_t>0).*(1-(security_p/2));
% security_p_fdr_one_tail=mafdr(security_p_one_tail,'BHFDR',0);
% hedonism_p_one_tail=(hedonism_t>0).*(hedonism_p/2)+(hedonism_t<0).*(1-(hedonism_p/2));
% hedonism_p_fdr_one_tail=mafdr(hedonism_p_one_tail,'BHFDR',0)
%% (Used)FDR & Reduce Comparisons (Focus only on subcortical regions)
% get the indexes of subcortical regions
subcortical_index=~cellfun('isempty',regexpi(aal_valid_name,...
    '(insula)|(amygdala)|(thalamus)|(putamen)|(caudate)|(pallidum)|(hippocampus)'));
subcortical_name=aal_valid_name(subcortical_index);
subcortical_center=ROI_valid_center(subcortical_index,:);
% Security
% DC~SVS (Only focus on subcortical reigion)
security_p_in_subcortical=security_p(subcortical_index);
% [security_p_in_subcortical_fdr security_q security_pi]=mafdr(security_p_in_subcortical,'BHFDR',0);
security_p_in_subcortical_fdr=mafdr(security_p_in_subcortical,'BHFDR',1);

security_t_in_subcortical=security_t(subcortical_index);% All positive effects
%Tidy Output
[num2cell(security_t_in_subcortical) num2cell(security_p_in_subcortical),...
    num2cell(security_p_in_subcortical_fdr),...
    cellstr(subcortical_name) mat2cell(subcortical_center,...
    [repmat(1,size(security_p_in_subcortical_fdr,1),1)],[1 1 1])]

% Get ihe indices of subcortical Security Effect hubs
security_p_in_subcortical_fdr_index=logical([]);
order=1;
for i=1:size(subcortical_index,1)
if subcortical_index(i)==1
   security_p_in_subcortical_fdr_index(i,1)=(security_p_in_subcortical_fdr(order)<0.05);
   order=order+1;
else
    security_p_in_subcortical_fdr_index(i,1)=logical(0);
end
end


% Hedonism
% DC~SVS (Only focus on subcortical reigion)
hedonism_p_in_subcortical=hedonism_p(subcortical_index);
% [hedonism_p_in_subcortical_fdr security_q security_pi]=mafdr(hedonism_p_in_subcortical,'BHFDR',0);
hedonism_p_in_subcortical_fdr=mafdr(hedonism_p_in_subcortical,'BHFDR',1);

hedonism_t_in_subcortical=hedonism_t(subcortical_index);
%Tidy Output
[num2cell(hedonism_t_in_subcortical) num2cell(hedonism_p_in_subcortical),...
    num2cell(hedonism_p_in_subcortical_fdr),...
    cellstr(subcortical_name) mat2cell(subcortical_center,...
    [repmat(1,size(hedonism_p_in_subcortical_fdr,1),1)],[1 1 1])]
% Get ihe indices of subcortical Hedonism Effect hubs
hedonism_p_in_subcortical_fdr_index=logical([]);
order=1;
for i=1:size(subcortical_index,1)
if subcortical_index(i)==1
   hedonism_p_in_subcortical_fdr_index(i,1)=(hedonism_p_in_subcortical_fdr(order)<0.05);
   order=order+1;
else
    hedonism_p_in_subcortical_fdr_index(i,1)=logical(0);
end
end
%% (Used)FDR & Reduce Comparisons (Focus only on SMA)
% get the indexes of sma regions
sma_index=~cellfun('isempty',regexpi(aal_valid_name,...
    'Supp_Motor_Area'));
sma_name=aal_valid_name(sma_index);
sma_center=ROI_valid_center(sma_index,:);

% Security
% DC~SVS (Only focus on sma reigion)
security_p_in_sma=security_p(sma_index);
% [security_p_in_sma_fdr security_q security_pi]=mafdr(security_p_in_sma,'BHFDR',0);
security_p_in_sma_fdr=mafdr(security_p_in_sma,'BHFDR',1);

security_t_in_sma=security_t(sma_index);% All positive effects
%Tidy Output
[num2cell(security_t_in_sma) num2cell(security_p_in_sma),...
    num2cell(security_p_in_sma_fdr),...
    cellstr(sma_name) mat2cell(sma_center,...
    [repmat(1,size(security_p_in_sma_fdr,1),1)],[1 1 1])]

% Get ihe indices of sma Security Effect hubs
security_p_in_sma_fdr_index=logical([]);
order=1;
for i=1:size(sma_index,1)
if sma_index(i)==1
   security_p_in_sma_fdr_index(i,1)=(security_p_in_sma_fdr(order)<0.05);
   order=order+1;
else
    security_p_in_sma_fdr_index(i,1)=logical(0);
end
end

% Hedonism
% DC~SVS (Only focus on sma reigion)
hedonism_p_in_sma=hedonism_p(sma_index);
% [hedonism_p_in_sma_fdr hedonism_q hedonism_pi]=mafdr(hedonism_p_in_sma,'BHFDR',0);
hedonism_p_in_sma_fdr=mafdr(hedonism_p_in_sma,'BHFDR',1);
hedonism_t_in_sma=hedonism_t(sma_index);% All positive effects
%Tidy Output
[num2cell(hedonism_t_in_sma) num2cell(hedonism_p_in_sma),...
    num2cell(hedonism_p_in_sma_fdr),...
    cellstr(sma_name) mat2cell(sma_center,...
    [repmat(1,size(hedonism_p_in_sma_fdr,1),1)],[1 1 1])]

% Get ihe indices of sma hedonism Effect hubs
hedonism_p_in_sma_fdr_index=logical([]);
order=1;
for i=1:size(sma_index,1)
if sma_index(i)==1
   hedonism_p_in_sma_fdr_index(i,1)=(hedonism_p_in_sma_fdr(order)<0.05);
   order=order+1;
else
    hedonism_p_in_sma_fdr_index(i,1)=logical(0);
end
end

%% Post-hoc (Focus only on subcortical regions)
%% See which links from a given seed could be prediced by SVS
%% Load R matrix(excluding cerebelum) and SVS info
cd('/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_3_R_and_Z_Matrix_of_160_ROI/2_R_D_and_Z_matrix/exclude_cerebelum')
load('collect_R_matrix_ROI_160_id_40.mat'); % R matrix(142x142) each per person
cd(path);
load('id_age_gender_sec_hed_sti.mat') %Load in SVS information
id=a(:,1);
id_valid_index=~ismember(id,{'S031' 'S040' 'S044'}');
age=str2double(a(id_valid_index,2));
gender=str2double(a(id_valid_index,3));
hedonism=str2double(a(id_valid_index,5));
security=str2double(a(id_valid_index,4));
%% Use Security subcortical-hubs as seeds
%% Security (R matrix of: Seeds ~ 142 ROI)
m=cell2mat(collect_R_matrix);%collpase the cell to (142 ROI) x (142x40=5680) matrix
% seed_index=find(overall_and_security_hub_index);%get the index of each interested center(seeds)
seed_index=subcortical_index;%get the index of each interested center(seeds)

seed_R_matrix=m(seed_index,:);
seed_R_cell=mat2cell(seed_R_matrix,repmat(1,1,size(seed_R_matrix,1)))%split the matrix by seed (turn to 15x [1x5680])
seed_R_cell=cellfun(@(x)reshape(x,[size(m,1), size(collect_R_matrix,2)]),seed_R_cell,'un',0);%turn into 15 (142x40) matrices

% Derive the R matrix which record the links between seeds and ROIs
%corrcoef(seed_R_cell{1}(1,:),security')
% cellfun(@(x) corrcoef(x,security'),mat2cell(seed_R_cell{1},repmat(1,1,142)),'un',0)
seed_security_160_ROI_R_matrix_r=[];
seed_security_160_ROI_R_matrix_p=[];
subindex=@(M,r,c) M{r,c};
n_seed=size(seed_R_cell,1);
for SEED=1:size(1:n_seed,2) %15 seeds
    for ROI_index=1:size(seed_R_cell{1},1) %142 ROI
    ROI=subindex(mat2cell(seed_R_cell{SEED},repmat(1,1,size(seed_R_cell{1},1))),ROI_index,1);
    [r p]=corrcoef(ROI,security');  
    r=r(1,2);
    p=p(1,2);
    seed_security_160_ROI_R_matrix_r(SEED,ROI_index)=r;
    seed_security_160_ROI_R_matrix_p(SEED,ROI_index)=p;
    end
end

% [seed_security_160_ROI_R_matrix_r seed_security_160_ROI_R_matrix_p]=...
%     cellfun(@(SEED) cellfun(@(ROI)...
%     corrcoef(ROI,security'),...
%     mat2cell(seed_R_cell{SEED},repmat(1,1,size(seed_R_cell{1},1))),'un',0),...
%     mat2cell(1:n_seed,1,repmat(1,1,size(1:n_seed,2))),'un',0);
%% Tidy Output of Security (R matrix of: Seeds ~ 142 ROI) (15 seeds)
collect_tidy_seed_security_160_ROI_R_matrix={};
for SEED=1:n_seed;
%Tidy Output (FDR,aal,mni etc, one table per seed)
[seed_security_160_ROI_R_matrix_p_fdr seed_security_q seed_security_pi]=...
    mafdr(seed_security_160_ROI_R_matrix_p(SEED,:),'BHFDR',0);
tidy_result_link_level=[cellstr(repmat(subcortical_name(SEED),1,size(aal_valid_name,1)))',...
    num2cell(seed_security_160_ROI_R_matrix_r(SEED,:))',...
    num2cell(seed_security_160_ROI_R_matrix_p(SEED,:))',...
    num2cell(seed_security_160_ROI_R_matrix_p_fdr)',...
    cellstr(aal_valid_name),...
    mat2cell(ROI_valid_center,[repmat(1,size(aal_valid_name,1),1)],[1 1 1])];
tidy_table_link_level=cell2table(tidy_result_link_level,...
    'VariableNames',{'seed_name','r','p_unc','p_fdr','aal_name','x','y','z'});
collect_tidy_seed_security_160_ROI_R_matrix(SEED,1)={sortrows(tidy_table_link_level,'p_fdr','ascend')};
end

% append seed info.(seed aal_name, p_FDR of seed, how many significant links from the seed, mni)
tidy_result_seed_level=...
    [num2cell(security_t_in_subcortical) num2cell(security_p_in_subcortical),...
    num2cell(security_p_in_subcortical_fdr),...
    cellstr(subcortical_name) mat2cell(subcortical_center,...
    [repmat(1,size(security_p_in_subcortical_fdr,1),1)],[1 1 1])];
tidy_table_seed_level=cell2table(tidy_result_seed_level,...
    'VariableNames',{'t','p_unc','p_fdr','aal_name','x','y','z'});
collect_tidy_seed_security_160_ROI_R_matrix=...
    [tidy_table_seed_level,...
    collect_tidy_seed_security_160_ROI_R_matrix,...
    ];
%Rename the links table
collect_tidy_seed_security_160_ROI_R_matrix.Properties.VariableNames{end}='links_table';
collect_tidy_seed_security_160_ROI_R_matrix=...
    sortrows(collect_tidy_seed_security_160_ROI_R_matrix,'p_unc','ascend');
%% Use Hedonism subcortical-hubs as seeds
%% Hedonism (R matrix of: Seeds ~ 142 ROI)
m=cell2mat(collect_R_matrix);%collpase the cell to (142 ROI) x (142x40=5680) matrix
% seed_index=find(overall_and_hedonism_hub_index);%get the index of each interested center(seeds)
seed_index=subcortical_index;%get the index of each interested center(seeds)

seed_R_matrix=m(seed_index,:);
seed_R_cell=mat2cell(seed_R_matrix,repmat(1,1,size(seed_R_matrix,1)));%split the matrix by seed (turn to 15x [1x5680])
seed_R_cell=cellfun(@(x)reshape(x,[size(m,1), size(collect_R_matrix,2)]),seed_R_cell,'un',0);%turn into 15 (142x40) matrices

% Derive the R matrix which record the links between seeds and ROIs
%corrcoef(seed_R_cell{1}(1,:),hedonism')
% cellfun(@(x) corrcoef(x,hedonism'),mat2cell(seed_R_cell{1},repmat(1,1,142)),'un',0)
seed_hedonism_160_ROI_R_matrix_r=[];
seed_hedonism_160_ROI_R_matrix_p=[];
subindex=@(M,r,c) M{r,c};
n_seed=size(seed_R_cell,1);
for SEED=1:size(1:n_seed,2) %15 seeds
    for ROI_index=1:size(seed_R_cell{1},1) %142 ROI
    ROI=subindex(mat2cell(seed_R_cell{SEED},repmat(1,1,size(seed_R_cell{1},1))),ROI_index,1);
    [r p]=corrcoef(ROI,hedonism');  
    r=r(1,2);
    p=p(1,2);
    seed_hedonism_160_ROI_R_matrix_r(SEED,ROI_index)=r;
    seed_hedonism_160_ROI_R_matrix_p(SEED,ROI_index)=p;
    end
end

% [seed_hedonism_160_ROI_R_matrix_r seed_hedonism_160_ROI_R_matrix_p]=...
%     cellfun(@(SEED) cellfun(@(ROI)...
%     corrcoef(ROI,hedonism'),...
%     mat2cell(seed_R_cell{SEED},repmat(1,1,size(seed_R_cell{1},1))),'un',0),...
%     mat2cell(1:n_seed,1,repmat(1,1,size(1:n_seed,2))),'un',0);
%% Tidy Output of Hedonism, (R matrix of: Seeds ~ 142 ROI) (15 seeds)
collect_tidy_seed_hedonism_160_ROI_R_matrix={};
for SEED=1:n_seed;
%Tidy Output (FDR,aal,mni etc, one table per seed)
[seed_hedonism_160_ROI_R_matrix_p_fdr seed_hedonism_q seed_hedonism_pi]=...
    mafdr(seed_hedonism_160_ROI_R_matrix_p(SEED,:),'BHFDR',0);
tidy_result_link_level=[cellstr(repmat(subcortical_name(SEED),1,size(aal_valid_name,1)))',...
    num2cell(seed_hedonism_160_ROI_R_matrix_r(SEED,:))',...
    num2cell(seed_hedonism_160_ROI_R_matrix_p(SEED,:))',...
    num2cell(seed_hedonism_160_ROI_R_matrix_p_fdr)',...
    cellstr(aal_valid_name),...
    mat2cell(ROI_valid_center,[repmat(1,size(aal_valid_name,1),1)],[1 1 1])];
tidy_table_link_level=cell2table(tidy_result_link_level,...
    'VariableNames',{'seed_name','r','p_unc','p_fdr','aal_name','x','y','z'});
collect_tidy_seed_hedonism_160_ROI_R_matrix(SEED,1)={sortrows(tidy_table_link_level,'p_fdr','ascend')};
end

% append seed info.(seed aal_name, p_FDR of seed, how many significant links from the seed, mni)
tidy_result_seed_level=...
    [num2cell(hedonism_t_in_subcortical) num2cell(hedonism_p_in_subcortical),...
    num2cell(hedonism_p_in_subcortical_fdr),...
    cellstr(subcortical_name) mat2cell(subcortical_center,...
    [repmat(1,size(hedonism_p_in_subcortical_fdr,1),1)],[1 1 1])];
tidy_table_seed_level=cell2table(tidy_result_seed_level,...
    'VariableNames',{'t','p_unc','p_fdr','aal_name','x','y','z'});
collect_tidy_seed_hedonism_160_ROI_R_matrix=...
    [tidy_table_seed_level,...
    collect_tidy_seed_hedonism_160_ROI_R_matrix,...
    ];
%Rename the links table
collect_tidy_seed_hedonism_160_ROI_R_matrix.Properties.VariableNames{end}='links_table';
collect_tidy_seed_hedonism_160_ROI_R_matrix=...
    sortrows(collect_tidy_seed_hedonism_160_ROI_R_matrix,'p_unc','ascend');
%% Focus only on Striatum
%% Security
t_sec=collect_tidy_seed_security_160_ROI_R_matrix...
    (~cellfun('isempty',regexpi(collect_tidy_seed_security_160_ROI_R_matrix.aal_name,...
    '(putamen)|(caudate)')),:);
t_sec.p_fdr=mafdr(t_sec.p_unc,'BHFDR',1);
%%H edonism
t_hed=collect_tidy_seed_hedonism_160_ROI_R_matrix...
    (~cellfun('isempty',regexpi(collect_tidy_seed_hedonism_160_ROI_R_matrix.aal_name,...
    '(putamen)|(caudate)')),:);
t_hed.p_fdr=mafdr(t_hed.p_unc,'BHFDR',1);
%% Post-hoc (Focus only on SMA regions)
%% See which links from a given seed could be prediced by SVS
%% Use Security sma-hubs as seeds
%% Security (R matrix of: Seeds ~ 142 ROI)
m=cell2mat(collect_R_matrix);%collpase the cell to (142 ROI) x (142x40=5680) matrix
% seed_index=find(overall_and_security_hub_index);%get the index of each interested center(seeds)
seed_index=sma_index;%get the index of each interested center(seeds)

seed_R_matrix=m(seed_index,:);
seed_R_cell=mat2cell(seed_R_matrix,repmat(1,1,size(seed_R_matrix,1)))%split the matrix by seed (turn to 15x [1x5680])
seed_R_cell=cellfun(@(x)reshape(x,[size(m,1), size(collect_R_matrix,2)]),seed_R_cell,'un',0);%turn into 15 (142x40) matrices

% Derive the R matrix which record the links between seeds and ROIs
%corrcoef(seed_R_cell{1}(1,:),security')
% cellfun(@(x) corrcoef(x,security'),mat2cell(seed_R_cell{1},repmat(1,1,142)),'un',0)
seed_security_160_ROI_R_matrix_r=[];
seed_security_160_ROI_R_matrix_p=[];
subindex=@(M,r,c) M{r,c};
n_seed=size(seed_R_cell,1);
for SEED=1:size(1:n_seed,2) %15 seeds
    for ROI_index=1:size(seed_R_cell{1},1) %142 ROI
    ROI=subindex(mat2cell(seed_R_cell{SEED},repmat(1,1,size(seed_R_cell{1},1))),ROI_index,1);
    [r p]=corrcoef(ROI,security');  
    r=r(1,2);
    p=p(1,2);
    seed_security_160_ROI_R_matrix_r(SEED,ROI_index)=r;
    seed_security_160_ROI_R_matrix_p(SEED,ROI_index)=p;
    end
end

% [seed_security_160_ROI_R_matrix_r seed_security_160_ROI_R_matrix_p]=...
%     cellfun(@(SEED) cellfun(@(ROI)...
%     corrcoef(ROI,security'),...
%     mat2cell(seed_R_cell{SEED},repmat(1,1,size(seed_R_cell{1},1))),'un',0),...
%     mat2cell(1:n_seed,1,repmat(1,1,size(1:n_seed,2))),'un',0);
%% Tidy Output of Security (R matrix of: Seeds ~ 142 ROI) (15 seeds)
collect_tidy_seed_security_160_ROI_R_matrix={};
for SEED=1:n_seed;
%Tidy Output (FDR,aal,mni etc, one table per seed)
[seed_security_160_ROI_R_matrix_p_fdr seed_security_q seed_security_pi]=...
    mafdr(seed_security_160_ROI_R_matrix_p(SEED,:),'BHFDR',0);
tidy_result_link_level=[cellstr(repmat(sma_name(SEED),1,size(aal_valid_name,1)))',...
    num2cell(seed_security_160_ROI_R_matrix_r(SEED,:))',...
    num2cell(seed_security_160_ROI_R_matrix_p(SEED,:))',...
    num2cell(seed_security_160_ROI_R_matrix_p_fdr)',...
    cellstr(aal_valid_name),...
    mat2cell(ROI_valid_center,[repmat(1,size(aal_valid_name,1),1)],[1 1 1])];
tidy_table_link_level=cell2table(tidy_result_link_level,...
    'VariableNames',{'seed_name','r','p_unc','p_fdr','aal_name','x','y','z'});
collect_tidy_seed_security_160_ROI_R_matrix(SEED,1)={sortrows(tidy_table_link_level,'p_fdr','ascend')};
end

% append seed info.(seed aal_name, p_FDR of seed, how many significant links from the seed, mni)
tidy_result_seed_level=...
    [num2cell(security_t_in_sma) num2cell(security_p_in_sma),...
    num2cell(security_p_in_sma_fdr),...
    cellstr(sma_name) mat2cell(sma_center,...
    [repmat(1,size(security_p_in_sma_fdr,1),1)],[1 1 1])];
tidy_table_seed_level=cell2table(tidy_result_seed_level,...
    'VariableNames',{'t','p_unc','p_fdr','aal_name','x','y','z'});
collect_tidy_seed_security_160_ROI_R_matrix=...
    [tidy_table_seed_level,...
    collect_tidy_seed_security_160_ROI_R_matrix,...
    ];
%Rename the links table
collect_tidy_seed_security_160_ROI_R_matrix.Properties.VariableNames{end}='links_table';
collect_tidy_seed_security_160_ROI_R_matrix=...
    sortrows(collect_tidy_seed_security_160_ROI_R_matrix,'p_unc','ascend');
%% Use Hedonism sma-hubs as seeds
%% Hedonism (R matrix of: Seeds ~ 142 ROI)
m=cell2mat(collect_R_matrix);%collpase the cell to (142 ROI) x (142x40=5680) matrix
% seed_index=find(overall_and_hedonism_hub_index);%get the index of each interested center(seeds)
seed_index=sma_index;%get the index of each interested center(seeds)

seed_R_matrix=m(seed_index,:);
seed_R_cell=mat2cell(seed_R_matrix,repmat(1,1,size(seed_R_matrix,1)));%split the matrix by seed (turn to 15x [1x5680])
seed_R_cell=cellfun(@(x)reshape(x,[size(m,1), size(collect_R_matrix,2)]),seed_R_cell,'un',0);%turn into 15 (142x40) matrices

% Derive the R matrix which record the links between seeds and ROIs
%corrcoef(seed_R_cell{1}(1,:),hedonism')
% cellfun(@(x) corrcoef(x,hedonism'),mat2cell(seed_R_cell{1},repmat(1,1,142)),'un',0)
seed_hedonism_160_ROI_R_matrix_r=[];
seed_hedonism_160_ROI_R_matrix_p=[];
subindex=@(M,r,c) M{r,c};
n_seed=size(seed_R_cell,1);
for SEED=1:size(1:n_seed,2) %15 seeds
    for ROI_index=1:size(seed_R_cell{1},1) %142 ROI
    ROI=subindex(mat2cell(seed_R_cell{SEED},repmat(1,1,size(seed_R_cell{1},1))),ROI_index,1);
    [r p]=corrcoef(ROI,hedonism');  
    r=r(1,2);
    p=p(1,2);
    seed_hedonism_160_ROI_R_matrix_r(SEED,ROI_index)=r;
    seed_hedonism_160_ROI_R_matrix_p(SEED,ROI_index)=p;
    end
end

% [seed_hedonism_160_ROI_R_matrix_r seed_hedonism_160_ROI_R_matrix_p]=...
%     cellfun(@(SEED) cellfun(@(ROI)...
%     corrcoef(ROI,hedonism'),...
%     mat2cell(seed_R_cell{SEED},repmat(1,1,size(seed_R_cell{1},1))),'un',0),...
%     mat2cell(1:n_seed,1,repmat(1,1,size(1:n_seed,2))),'un',0);
%% Tidy Output of Hedonis, (R matrix of: Seeds ~ 142 ROI) (15 seeds)
collect_tidy_seed_hedonism_160_ROI_R_matrix={};
for SEED=1:n_seed;
%Tidy Output (FDR,aal,mni etc, one table per seed)
[seed_hedonism_160_ROI_R_matrix_p_fdr seed_hedonism_q seed_hedonism_pi]=...
    mafdr(seed_hedonism_160_ROI_R_matrix_p(SEED,:),'BHFDR',0);
tidy_result_link_level=[cellstr(repmat(sma_name(SEED),1,size(aal_valid_name,1)))',...
    num2cell(seed_hedonism_160_ROI_R_matrix_r(SEED,:))',...
    num2cell(seed_hedonism_160_ROI_R_matrix_p(SEED,:))',...
    num2cell(seed_hedonism_160_ROI_R_matrix_p_fdr)',...
    cellstr(aal_valid_name),...
    mat2cell(ROI_valid_center,[repmat(1,size(aal_valid_name,1),1)],[1 1 1])];
tidy_table_link_level=cell2table(tidy_result_link_level,...
    'VariableNames',{'seed_name','r','p_unc','p_fdr','aal_name','x','y','z'});
collect_tidy_seed_hedonism_160_ROI_R_matrix(SEED,1)={sortrows(tidy_table_link_level,'p_fdr','ascend')};
end

% append seed info.(seed aal_name, p_FDR of seed, how many significant links from the seed, mni)
tidy_result_seed_level=...
    [num2cell(hedonism_t_in_sma) num2cell(hedonism_p_in_sma),...
    num2cell(hedonism_p_in_sma_fdr),...
    cellstr(sma_name) mat2cell(sma_center,...
    [repmat(1,size(hedonism_p_in_sma_fdr,1),1)],[1 1 1])];
tidy_table_seed_level=cell2table(tidy_result_seed_level,...
    'VariableNames',{'t','p_unc','p_fdr','aal_name','x','y','z'});
collect_tidy_seed_hedonism_160_ROI_R_matrix=...
    [tidy_table_seed_level,...
    collect_tidy_seed_hedonism_160_ROI_R_matrix,...
    ];
%Rename the links table
collect_tidy_seed_hedonism_160_ROI_R_matrix.Properties.VariableNames{end}='links_table';
collect_tidy_seed_hedonism_160_ROI_R_matrix=...
    sortrows(collect_tidy_seed_hedonism_160_ROI_R_matrix,'p_unc','ascend');