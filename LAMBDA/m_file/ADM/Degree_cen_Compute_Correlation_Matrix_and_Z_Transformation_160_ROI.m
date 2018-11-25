%% Part 6-A: Degree Centrality (160 ROI approach): Compute Correlation Matrix(R),degree of each ROI(D) and perform Z transformation(Z)
%% Two Mode to choose: Exclude Crebelum or NOT
path='/bml/Data/Bank6/ADM-YunShiuan';
cd(path);
addpath(strcat(path,'/Scripts(m)_files'));% to enable using the function "cor2mni"

id=csvread(strcat(path,'/young_idlist.csv')); %read only the young's IDs
young_name=strcat('S0',num2str(id));
young_name=cellstr(young_name); %convert char to cell to enable choosing string by index
young_name_40=cellstr(strvcat(young_name{~ismember(young_name,{'S031' 'S040' 'S044'}')})); %Exclude S031,S040,S044

%% Load in the .mat of Residuals in the 160 ROI of 40 people
cd('/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_3_R_and_Z_Matrix_of_160_ROI/1_raw_residual_160_ROI');
load('beta_160_ROI.mat','beta_160_ROI');

cd('/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_3_R_and_Z_Matrix_of_160_ROI/2_R_D_and_Z_matrix');
%% Exclude cerebelum (use it only if needed)
load('/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_4_Regression_degree_cen_by_SVS_by_R/lm_result/without_cerebelum/cerebelum_aal_name.mat')
cd('/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_4_Regression_degree_cen_by_SVS_by_R/lm_result')
load('collect_lm_160_ROI_by_2_SVS.mat','aal_name');
cerebelum_index=ismember(aal_name,cerebelum_aal_name);
cd('/bml/Data/Bank6/ADM-YunShiuan/degree_centrality/Part_3_R_and_Z_Matrix_of_160_ROI/2_R_D_and_Z_matrix/exclude_cerebelum')
beta_160_ROI=beta_160_ROI(:,~cerebelum_index);

%Conver cell to 1090 0bservations x 1160 ROI matrix (142 ROI if excluding cerebelum)
collect_R_matrix={};
for id=1:size(beta_160_ROI,1);
df=cell2mat(beta_160_ROI(id,:)); % 1090 0bservations x 160 ROI matrix
collect_R_matrix{1,id}=corrcoef(df);
id
end
save('collect_R_matrix_ROI_160_id_40','collect_R_matrix');
% Check how many positive and negative correlation for each person
% pn_ratio={};
% for id=1:size(beta_160_ROI,1);
% pn_ratio{1,id}=sum(sum(collect_R_matrix{id}>0 & collect_R_matrix{id}~=1))/sum(sum(collect_R_matrix{id}<0));
% end
% mean(cell2mat(pn_ratio))% positive-negative ratio mean=0.95
%% Compute D matrix ('degree' of each ROI)
%criteria=0.25
collect_D_matrix={};
for id=1:size(beta_160_ROI,1);
d=sum(collect_R_matrix{1,id}>0.25 & collect_R_matrix{1,id}~=1,1);%Sum for each column(each ROI has a D value)
collect_D_matrix{id,1}=d;clea
id
end
collect_D_matrix=cell2mat(collect_D_matrix); % 40id x 160 ROI matrix
save('collect_D_matrix_ROI_160_id_40','collect_D_matrix');
%% Derive Z matrix (strandardize across 160 ROI within the person)
collect_Z_matrix=zscore(collect_D_matrix,1,2); %Standardise across row(160 ROI)
%mean(collect_Z_matrix,2)
%std(collect_Z_matrix,0,2)
save('collect_Z_matrix_ROI_160_id_40','collect_Z_matrix');


