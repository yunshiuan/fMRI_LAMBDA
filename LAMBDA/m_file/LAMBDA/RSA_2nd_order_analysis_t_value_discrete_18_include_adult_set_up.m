%% Part 3&4 - iv-c: RSA : 2nd order analysis - Set up environments for 
% Set up environment for MDS and relatedness test.
% This script contain reletively stable parameter and preprocess.
%
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
% PATTERN_REDUCED_CONCEPTUAL_RDMS=['^g$|^n$|null|^ad$|^sd$|',...
%                                  '(?<=ad|sd)3$|(?<=ad|sd)5$|(?<=ad|sd)7$'];

PATTERN_REDUCED_CONCEPTUAL_RDMS=PATTERN_BASED_CONCEPTUAL_RDMS;
% Load in needed data
load(FILE_1ST_USER_OPTION,'user_options_visualize_RDM');
load(FILE_RDM_RESULT,...
     'neural_RDMs',...
     'neural_RDMs_avg_2_grade','neural_RDMs_avg_5_grade',...
     'neural_RDMs_avg_child','neural_RDMs_avg_adult','neural_RDMs_avg_all',...
     'conceptual_RDMs',...
     'user_options_RDM_construction')
list_neural_avg_RDMs={neural_RDMs_avg_2_grade,neural_RDMs_avg_5_grade,neural_RDMs_avg_adult};

% Options for 2nd order analysis
clear user_options_2nd_order
% Edit the color bar color due to Matlab version issue
%(R2015b change default color maps. Should now be colormap(jet(256)))
% user_options_visualize_RDM.colourScheme=colormap(jet(256)); [Already edify the source codes]
user_options_2nd_order=user_options_visualize_RDM;
user_options_2nd_order.analysisName = VERSION_DATE;
user_options_2nd_order.rootPath=PATH_2ND_RESULT;                   


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
    elseif~isempty(regexp(model_name,'^agv\d$','match'))
        model_colors{model_index}=[0.8,0.5,0];
    end
end

for model_index=1:numel(conceptual_RDMs)  
    conceptual_RDMs(model_index).color=...
        model_colors{model_index};
end

% Derive the conceptual RDMs to be included for each age group
list_reduced_conceptual_RDMs_names=...
    {'reduced_conceptual_RDMs_2_grade',...
     'reduced_conceptual_RDMs_5_grade',...
     'reduced_conceptual_RDMs_adult'};
full_conceptual_RDMS_names={conceptual_RDMs.name}';

for group_index=1:numel(LIST_AGE_GROUP_NAME)
    
    if(BOOLEAN_ACC_INCLUDE)
        %(NOT care about it now) Include the age-specific difficulty models (pure and hybrid)
        pattern_reduced_for_age_group=[PATTERN_REDUCED_CONCEPTUAL_RDMS,...
            '|ac',LIST_AGE_ABRV{group_index},'$'];
        
        %The version with ACC-hybride models
        %         pattern_reduced_for_age_group=[PATTERN_REDUCED_CONCEPTUAL_RDMS,...
        %             '|ac',LIST_AGE_ABRV{group_index},'$',...
        %             '|ac',LIST_AGE_ABRV{group_index},'(?=3|5|7)'];
    else
        pattern_reduced_for_age_group=PATTERN_REDUCED_CONCEPTUAL_RDMS;
    end
    
    reduced_conceptual_RMDS_index=...
        ~cellfun('isempty',...
        regexp(full_conceptual_RDMS_names,pattern_reduced_for_age_group,'match'));
    reduced_conceptual_RDMs=conceptual_RDMs(1,reduced_conceptual_RMDS_index);
    
    assignin('base',list_reduced_conceptual_RDMs_names{group_index},reduced_conceptual_RDMs);
end
list_reduced_conceptual_RDMs={reduced_conceptual_RDMs_2_grade,...
                              reduced_conceptual_RDMs_5_grade,...
                              reduced_conceptual_RDMs_adult};

% %Check the names of the reduced set of models
% {reduced_conceptual_RDMs.name}'

% Read in subject and run inclusion index info ------------------------------------
%Child
run_inclusion_child_index=read_mixed_csv_to_table(FILE_VALID_RUN_CHILD);
subject_list_child=unique(run_inclusion_child_index.sub_id);% Derive subjects with valid runs
%Adult
run_inclusion_adult_index=read_mixed_csv_to_table(FILE_VALID_RUN_ADULT);
subject_list_adult=unique(run_inclusion_adult_index.sub_id);% Derive subjects with valid runs
% subject_list_adult=subject_list_adult(~ismember(subject_list_adult,{'XFC305'}));
% %305 have missing run (Already resolved)

% Read in demographic data (for age group split)------------------------
%Child
table_demographic_child=read_mixed_csv_to_table(FILE_DEMOGRAPHIC_CHILD);
table_demographic_child=table_demographic_child(:,{'sub_id','grade'});
table_demographic_child.sub_id=lower(table_demographic_child.sub_id);
table_demographic_child=join(table(subject_list_child,'VariableNames',{'sub_id'}),...
    table_demographic_child);
%Adult
table_demographic_adult=read_mixed_csv_to_table(FILE_DEMOGRAPHIC_ADULT);
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
list_group_index={index_2_grade,index_5_grade,index_adult}; % Should be in line with 'LIST_GROUP_NAME'
% LIST_GROUP_NAME={'2_grade','5_grade','adult'};