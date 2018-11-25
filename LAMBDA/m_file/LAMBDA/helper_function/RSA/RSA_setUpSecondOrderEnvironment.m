function [preprocessed_parameters] = RSA_setUpSecondOrderEnvironment(raw_parameters)
%[preprocessed_parameters] = RSA_setUpSecondOrderEnvironment(raw_parameters)
%
% Set up environment for MDS and relatedness test.
% This script contain reletively stable parameter and preprocess.
% %Files
% raw_parameters.file_user_options;
% raw_parameters.file_rdm_result;
% raw_parameters.file_valid_run_child;
% raw_parameters.file_valid_run_adult;
% raw_parameters.file_demographic_child;
% raw_parameters.file_demographic_adult;
% %Paths
% raw_parameters.path_2nd_result;
% 
% %Parameters
% raw_parameters.version_date;
% raw_parameters.pattern_reduced_conceptual_RDMs;
% raw_parameters.list_age_abrv;
% raw_parameters.list_age_group_name;
% raw_parameters.acc_include;
% raw_parameters.list_interested_brain_region;

%% Part 3&4 - iv-c: RSA : 2nd order analysis - Set up environments 
% trial 5 in the README.txt.
% Part3 & 4: 2nd order analysis - MDS and relatedness test
%
% iv: RDA_ID_discrete 18 ID:
% - RSA ID : 1~18 (3 formats x 3 distance levels x 2 directions)
% - The motivation is to increase the statistical efficiency (Conditions in the Trial 4 is too noisy),
% while still keep that same-notation trials have the same statistical efficiency
% as the cross-notation trials (as Trial 4 did).
% [LL: 12 trials = 3 distance levels x 2 direction x 2 repetition]
% [FF: 12 trials = 3 distance levels x 2 direction x 2 repetition]
% [cross-notation: 12 trials = 3 distance levels x 2 direction x 2 [repetition]
% (Having more repetition in a condition entails a smaller degree of random noise.)

% The motivation for having 18 trials is to ensure beta estimate of each condition is subjected to same amount of noise.
% (Having more repetition in a condition entails a smaller degree of random noise.)

% c: Include adult neural data (both collapse split by age group)
% Split into 2nd grader, 5th grader, and the adult

% [Run info has already been averaged out. Therefore, the number of runs no longer matters.]

%% Set up
% Constants----------------------------------------------

% The conceptual RDMS to be included in the 2nd analysis
% (MDS, model selection)
% pattern_reduced_conceptual_RDMs=['^g$|^n$|null|^ad$|^sd$|',...
%                                  '(?<=ad|sd)3$|(?<=ad|sd)5$|(?<=ad|sd)7$'];

%Get the vraibles from the fields
%Files
file_user_option=raw_parameters.file_user_options;
file_rdm_result=raw_parameters.file_rdm_result;
file_valid_run_child=raw_parameters.file_valid_run_child;
file_valid_run_adult=raw_parameters.file_valid_run_adult;
file_demographic_child=raw_parameters.file_demographic_child;
file_demographic_adult=raw_parameters.file_demographic_adult;
%Paths
path_2nd_result=raw_parameters.path_2nd_result;

%Parameters
version_date=raw_parameters.version_date;
pattern_reduced_conceptual_RDMs=raw_parameters.pattern_reduced_conceptual_RDMs;
list_age_abrv=raw_parameters.list_age_abrv;
list_age_group_name=raw_parameters.list_age_group_name;
boolean_acc_include=raw_parameters.acc_include;
list_interested_brain_region=raw_parameters.list_interested_brain_region;

% Load in needed data-----------------------------------------------
%load in user option
user_option=load(file_user_option);
user_option_field_name=regexp(fieldnames(user_option),'user.*','match');
user_option=user_option.(char(user_option_field_name{1}));

%load in neural RDMs
load(file_rdm_result,...
     'neural_RDMs',...
     'neural_RDMs_avg_2_grade','neural_RDMs_avg_5_grade',...
     'neural_RDMs_avg_child','neural_RDMs_avg_adult','neural_RDMs_avg_all',...
     'conceptual_RDMs',...
     'user_options_RDM_construction')
pattern_neural_RDMs_avg=strcat('neural_RDMs_avg_','(',...
    strjoin({list_age_group_name{:}}','|'),...
    ')');
neural_RDMs_avg_interested=regexp(who',pattern_neural_RDMs_avg,'match');
list_neural_avg_RDMs=neural_RDMs_avg_interested(~cellfun(@isempty,...
    neural_RDMs_avg_interested));
list_neural_avg_RDMs=...
    eval(strcat('{',strjoin([list_neural_avg_RDMs{:}],','),'}'));
% list_neural_avg_RDMs={neural_RDMs_avg_2_grade,neural_RDMs_avg_5_grade,neural_RDMs_avg_adult};

% Only filter in the interested brain region-----------------------------
%The list_neural_avg_RDMs
for age_group_index=1:numel(list_age_group_name)
    raw_brain_region_names={list_neural_avg_RDMs{age_group_index}.name};
    interested_brain_region_index=ismember(raw_brain_region_names,list_interested_brain_region);
    list_neural_avg_RDMs{age_group_index}=list_neural_avg_RDMs{age_group_index}(interested_brain_region_index);
end
%The neural_RDMs
pattern_interested_brain_region=strjoin(list_interested_brain_region,'|');
pattern_interested_brain_region=regexprep(pattern_interested_brain_region,'_','-');
neural_RDMs_name_reshape=reshape({neural_RDMs.name},size(neural_RDMs));
neural_RDMs_interested_brain_region_index=...
    ~cellfun('isempty',regexp(neural_RDMs_name_reshape,pattern_interested_brain_region,'match'));
n_voi=numel(list_interested_brain_region);
n_sub=size(neural_RDMs,2);
neural_RDMs_interested_brain_region=neural_RDMs(neural_RDMs_interested_brain_region_index);
neural_RDMs_interested_brain_region=reshape(neural_RDMs_interested_brain_region,n_voi,n_sub);

% Options for 2nd order analysis-----------------------------------------
clear user_options_2nd_order
% Edit the color bar color due to Matlab version issue
%(R2015b change default color maps. Should now be colormap(jet(256)))
% user_options_visualize_RDM.colourScheme=colormap(jet(256)); [Already edify the source codes]
user_options_2nd_order=user_option;
user_options_2nd_order.analysisName = version_date;
user_options_2nd_order.rootPath=path_2nd_result;                   


% Recolor the conceptual models
model_colors=cell(size(conceptual_RDMs));
for model_index=1:numel(conceptual_RDMs)
    model_name=conceptual_RDMs(model_index).name;
    %FCF and null model (pink)
    if ~isempty(regexp(model_name,'^f$|null','match'))
        model_colors{model_index}=[1,0.7,1];
        
    %GCF (dark green)
    elseif ~isempty(regexp(model_name,'^g$','match'))
        model_colors{model_index}=[0,0.5,0];
    %GCF weighted by LL-distance (light green)
    elseif ~isempty(regexp(model_name,'^gl','match'))
        model_colors{model_index}=[0.5,1,0.5];
    %GCF weighted by distance (green)
    elseif ~isempty(regexp(model_name,'^g((a)|(s))','match'))
        model_colors{model_index}=[0,1,0];
    %GCF weighted by difficulty (yellowish green)
    elseif ~isempty(regexp(model_name,'^gac','match'))
        model_colors{model_index}=[0.5,0.7,0];
        
    %NCF (dark blue)
    elseif ~isempty(regexp(model_name,'^n$','match'))
        model_colors{model_index}=[0,0,0.5];
    %NCF weighted by LL-distance (light blue)
    elseif ~isempty(regexp(model_name,'^nl','match'))
        model_colors{model_index}=[0.5,0.9,1];
    %NCF weighted by distance (blue)
    elseif ~isempty(regexp(model_name,'^n((a)|(s))','match'))
        model_colors{model_index}=[0,0,1];
    %NCF weighted by difficulty (purplish blue)
    elseif ~isempty(regexp(model_name,'^nac','match'))
        model_colors{model_index}=[0.4,0,1];
        
    %Pure distance models (brown)
    elseif~isempty(regexp(model_name,'^((a)|(s))$','match')) 
        model_colors{model_index}=[1,0.7,0];
    %Pure difficulty models (dark red)
    elseif~isempty(regexp(model_name,'^ac\w$','match')) 
        model_colors{model_index}=[0.8,0,0];
    %Pure log distance model (dark brown)
    elseif~isempty(regexp(model_name,'^(a|s)gv\d$','match'))
        model_colors{model_index}=[0.8,0.5,0];
    end
end

for model_index=1:numel(conceptual_RDMs)  
    conceptual_RDMs(model_index).color=...
        model_colors{model_index};
end

% Derive the conceptual RDMs to be included for each age group------------
list_reduced_conceptual_RDMs_names=...
    {'reduced_conceptual_RDMs_2_grade',...
     'reduced_conceptual_RDMs_5_grade',...
     'reduced_conceptual_RDMs_adult'};
full_conceptual_RDMS_names={conceptual_RDMs.name}';

for group_index=1:numel(list_age_group_name)
    
    if(boolean_acc_include)
        %(NOT care about ACC now) Include the age-specific difficulty models (pure and hybrid)
        pattern_reduced_for_age_group=[pattern_reduced_conceptual_RDMs,...
            '|ac',list_age_abrv{group_index},'$'];
        
        %The version with ACC-hybrid models
        %         pattern_reduced_for_age_group=[pattern_reduced_conceptual_RDMs,...
        %             '|ac',LIST_AGE_ABRV{group_index},'$',...
        %             '|ac',LIST_AGE_ABRV{group_index},'(?=3|5|7)'];
    else
        pattern_reduced_for_age_group=pattern_reduced_conceptual_RDMs;
    end
    
    reduced_conceptual_RMDS_index=...
        ~cellfun('isempty',...
        regexp(full_conceptual_RDMS_names,pattern_reduced_for_age_group,'match'));
    reduced_conceptual_RDMs=conceptual_RDMs(1,reduced_conceptual_RMDS_index);
    
    list_reduced_conceptual_RDMs{group_index}=reduced_conceptual_RDMs;
%     assignin('base',list_reduced_conceptual_RDMs_names{group_index},reduced_conceptual_RDMs);
end

% list_reduced_conceptual_RDMs={reduced_conceptual_RDMs_2_grade,...
%                               reduced_conceptual_RDMs_5_grade,...
%                               reduced_conceptual_RDMs_adult};

% %Check the names of the reduced set of models
% {reduced_conceptual_RDMs.name}'

% Read in subject and run inclusion index info ------------------------------------
%Child
run_inclusion_child_index=read_mixed_csv_to_table(file_valid_run_child);
subject_list_child=unique(run_inclusion_child_index.sub_id);% Derive subjects with valid runs
%Adult
run_inclusion_adult_index=read_mixed_csv_to_table(file_valid_run_adult);
subject_list_adult=unique(run_inclusion_adult_index.sub_id);% Derive subjects with valid runs
% subject_list_adult=subject_list_adult(~ismember(subject_list_adult,{'XFC305'}));
% %305 have missing run (Already resolved)

% Read in demographic data (for age group split)------------------------
%Child
table_demographic_child=read_mixed_csv_to_table(file_demographic_child);
table_demographic_child=table_demographic_child(:,{'sub_id','grade'});
table_demographic_child.sub_id=lower(table_demographic_child.sub_id);
table_demographic_child=join(table(subject_list_child,'VariableNames',{'sub_id'}),...
    table_demographic_child);
%Adult
table_demographic_adult=read_mixed_csv_to_table(file_demographic_adult);
table_demographic_adult=table_demographic_adult(:,{'sub_id','grade'});
table_demographic_adult=join(table(subject_list_adult,'VariableNames',{'sub_id'}),...
    table_demographic_adult);
%Join the Two demographic Table (Child+Adult)
table_demographic_all=vertcat(table_demographic_child,table_demographic_adult);
index_all=true(numel(table_demographic_all.grade),1);
index_2_grade=strcmp(table_demographic_all.grade,'2');
index_5_grade=strcmp(table_demographic_all.grade,'5');
index_child=logical(index_2_grade+index_5_grade);%All children includind 2 and 5 graders
index_adult=strcmp(table_demographic_all.grade,'Adult');
subject_list=table_demographic_all.sub_id;

% list_group_index={index_2_grade,index_5_grade,index_adult}; % Should be in line with 'LIST_GROUP_NAME'

pattern_group_index=strcat('index_','(',...
    strjoin({list_age_group_name{:}}','|'),...
    ')');
group_index_interested=regexp(who',pattern_group_index,'match');
list_group_index=group_index_interested(~cellfun(@isempty,...
    group_index_interested));
list_group_index=...
    eval(strcat('{',strjoin([list_group_index{:}],','),'}'));

%Save to output
preprocessed_parameters.list_neural_avg_RDMs=list_neural_avg_RDMs;
preprocessed_parameters.list_reduced_conceptual_RDMs=list_reduced_conceptual_RDMs;
preprocessed_parameters.list_group_index=list_group_index;
preprocessed_parameters.subject_list=subject_list;
preprocessed_parameters.list_group_index=list_group_index;
preprocessed_parameters.list_age_group_name=list_age_group_name;
preprocessed_parameters.user_options_2nd_order=user_options_2nd_order;
preprocessed_parameters.neural_RDMs=neural_RDMs_interested_brain_region;
end
