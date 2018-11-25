%% Part 2 - iv-c: RSA : 1st order analysis 
% trial 6 in the README.txt.
% Part2: 1st order analysis (visualize RDM, MDS, clustering)
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

% [Run info has already been averaged out. Therefore, the number of runs no longer matters.]

% c: Include adult neural data (both collapse split by age group)
% Split into 2nd grader, 5th grader, and the adult

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Update subject list for children (i.e., inclusive_runs_indexes_new_June_10.csv) (2018/6/10)
%Update subject list for children (i.e., inclusive_runs_indexes_new_June_11.csv) (2018/6/10)

%% Set up
PATH_ROOT='D:\Yun-Shiuan_LAMBDA';
PATH_TOOL_CODE='D:\GoogleDrive\Lambda_code\m_file\LAMBDA\tool_code';
PATH_HELPER_FUNCTION='D:\GoogleDrive\Lambda_code\m_file\LAMBDA\helper_function';
addpath(genpath(fullfile(PATH_ROOT,'rsatoolbox')));
addpath(PATH_TOOL_CODE);% enable read_mixed_csv()
addpath(genpath(PATH_HELPER_FUNCTION));% enable my rsa helper functions

% Constants---------------------------------------------
%(the input of this script)
PATH_RDM_RESULT=fullfile(PATH_ROOT,'RSA','trial_22','Part1_neural_and_conceptual_RDMs');
%(the output of this script)
PATH_1ST_RESULT=fullfile(PATH_ROOT,'RSA','trial_22','Part2_1st_order_analysis');

%files
%(the input of this script)
FILE_RDM_RESULT=fullfile(PATH_RDM_RESULT,'RDMs','RDMs_and_options_neural_and_conceptual_trial_22_10-Jun-2018.mat');
FILE_VALID_RUN_CHILD=fullfile(PATH_ROOT,'Run_inclusion_info','inclusive_runs_indexes_new_June_11.csv');
FILE_VALID_RUN_ADULT=fullfile(PATH_ROOT,'Adult','Run_inclusion_info','inclusive_runs_indexes.csv');
FILE_DEMOGRAPHIC_CHILD=fullfile(PATH_ROOT,'demographic','tidy_demographic.csv');
FILE_DEMOGRAPHIC_ADULT=fullfile(PATH_ROOT,'Adult','demographic','tidy_demographic.csv');
% Parameters----------------------
LIST_VOI_NAME={"R_V1","L_V1",...
               "R_IPS","L_IPS","R_DLPFC"}';
LIST_GROUP_NAME={'all_sub','2_grade','5_grade','adult','child'};
LIST_MDS_CRITERION={'stress','metricstress'}; % The options for MDS
NUM_TRIALS_EACH_FORMAT=6; % 12 trials per format but only 6 are unique
NUM_TRIALS=18;
NUM_FORMAT=3;
NUM_DISTANCE_LEVEL=3;
CONDITION_LABELS=cellstr(string(1:NUM_TRIALS));
VERSION_NUM='22';
BOOLEAN_COLOR_GRADIENT_1st_order=true; % Could be adjusted if needed.
BOOLEAN_PLOT_RDM=true;
BOOLEAN_1ST_ANALYSYS=true;

% Load in needed data--------------------------------------------------
load(FILE_RDM_RESULT,...
     'neural_RDMs',...
     'neural_RDMs_avg_2_grade','neural_RDMs_avg_5_grade',...
     'neural_RDMs_avg_child','neural_RDMs_avg_adult','neural_RDMs_avg_all',...
     'conceptual_RDMs',...
     'user_options_RDM_construction')
list_neural_RDMs={neural_RDMs_avg_all,neural_RDMs_avg_2_grade,...
                  neural_RDMs_avg_5_grade,neural_RDMs_avg_adult,neural_RDMs_avg_child};
% Read in run inclusion index info ------------------------------------
%Child
run_inclusion_child_index=read_mixed_csv_to_table(FILE_VALID_RUN_CHILD);
subject_list_child=unique(run_inclusion_child_index.sub_id);% Derive subjects with valid runs
%Adult
run_inclusion_adult_index=read_mixed_csv_to_table(FILE_VALID_RUN_ADULT);
subject_list_adult=unique(run_inclusion_adult_index.sub_id);% Derive subjects with valid runs
% subject_list_adult=subject_list_adult(~ismember(subject_list_adult,{'XFC305'}));
% %305 have missing run (Already solved)

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
index_2_grade=strcmp(table_demographic_all.grade,'2');
index_5_grade=strcmp(table_demographic_all.grade,'5');
index_child=logical(index_2_grade+index_5_grade);%All children includind 2 and 5 graders
index_adult=strcmp(table_demographic_all.grade,'Adult');
subject_list=table_demographic_all.sub_id;

% Options for visualization----------------------------------------------
clear user_options_visualize_RDM
% Edit the color bar color due to Matlab version issue
%(R2015b change default color maps. Should now be colormap(jet(256)))
% user_options_visualize_RDM.colourScheme=colormap(jet(256)); [Already edify the source codes]
user_options_visualize_RDM=user_options_RDM_construction;
user_options_visualize_RDM.saveFiguresPDF=true;
user_options_visualize_RDM.saveFigurePDF=true;
user_options_visualize_RDM.tightInset=true;
user_options_visualize_RDM.analysisName = ['RSA_trial_',VERSION_NUM,'_',date];
user_options_visualize_RDM.conditionLabels=cellstr(string(1:NUM_TRIALS));
user_options_visualize_RDM.rootPath=PATH_RDM_RESULT;
user_options_visualize_RDM.displayFigures=false;
user_options_visualize_RDM.rubberbands=false;
% Color grandient for conditions to reflect the value of absolute distance
%Bright->Dark: Small Distance -> Large distance
%  If not needed, could be set as
if (BOOLEAN_COLOR_GRADIENT_1st_order)
    color_FF=brewermap(NUM_DISTANCE_LEVEL*2,'Reds');
    color_FF=color_FF([2,4,5],:);
    color_CN=brewermap(NUM_DISTANCE_LEVEL*2,'Greens');
    color_CN=color_CN([2,4,5],:);
    color_LL=brewermap(NUM_DISTANCE_LEVEL*2,'Blues');
    color_LL=color_LL([2,4,5],:);
    
    color_FF=repelem(color_FF,2,1); 
    color_CN=repelem(color_CN,2,1);
    color_LL=repelem(color_LL,2,1);
    
else % No gradient
    color_FF=repmat([1 0 0],NUM_TRIALS_EACH_FORMAT,1);
    color_CN=repmat([0 1 0],NUM_TRIALS_EACH_FORMAT,1);
    color_LL=repmat([0 0 1],NUM_TRIALS_EACH_FORMAT,1);
end
user_options_visualize_RDM.conditionColours=...
    [color_FF;...
    color_CN;...
    color_LL];
%% (Step 1) Visulize the Neural and Conceptual RDMs: Each subject and average across subjects
% Save the RDM matrix plot into the RDMs folder instead of the 1st order
% one
if(BOOLEAN_PLOT_RDM)
    % Plot Neural RDMs--------------------------------
    % Avg subjects-------------------
    % all subjects
    figureRDMs(concatenateRDMs(neural_RDMs_avg_all), user_options_visualize_RDM,...
        struct('fileName', 'neural_RDMs_avg_all'));
    % 2nd grade
    figureRDMs(concatenateRDMs(neural_RDMs_avg_2_grade), user_options_visualize_RDM,...
        struct('fileName', 'neural_RDMs_avg_2_grade'));
    % 5th grade
    figureRDMs(concatenateRDMs(neural_RDMs_avg_5_grade), user_options_visualize_RDM,...
        struct('fileName', 'neural_RDMs_avg_5_grade' ));
    % All childred
    figureRDMs(concatenateRDMs(neural_RDMs_avg_child), user_options_visualize_RDM,...
        struct('fileName', 'neural_RDMs_avg_child' ));
    % Adult
    figureRDMs(concatenateRDMs(neural_RDMs_avg_adult), user_options_visualize_RDM,...
        struct('fileName', 'neural_RDMs_avg_adult' ));
    % For each VOI-------------------------------------
    % all subjects----------------
    for v=1:size(LIST_VOI_NAME,1)
        figureRDMs(concatenateRDMs(neural_RDMs(v,:)), user_options_visualize_RDM,...
            struct('fileName', ['neural_RDMs_',char(LIST_VOI_NAME{v}),'_all']));
    end
    % 2nd grade----------------
    for v=1:size(LIST_VOI_NAME,1)
        figureRDMs(concatenateRDMs(neural_RDMs(v,index_2_grade)), user_options_visualize_RDM,...
            struct('fileName', ['neural_RDMs_',char(LIST_VOI_NAME{v}),'_2_grade']));
    end
    % 5th grade----------------
    for v=1:size(LIST_VOI_NAME,1)
        figureRDMs(concatenateRDMs(neural_RDMs(v,index_5_grade)), user_options_visualize_RDM,...
            struct('fileName', ['neural_RDMs_',char(LIST_VOI_NAME{v}),'_5_grade']));
    end
    % adult grade----------------
    for v=1:size(LIST_VOI_NAME,1)
        figureRDMs(concatenateRDMs(neural_RDMs(v,index_adult)), user_options_visualize_RDM,...
            struct('fileName', ['neural_RDMs_',char(LIST_VOI_NAME{v}),'_adult']));
    end
    % child----------------
    for v=1:size(LIST_VOI_NAME,1)
        figureRDMs(concatenateRDMs(neural_RDMs(v,index_child)), user_options_visualize_RDM,...
            struct('fileName', ['neural_RDMs_',char(LIST_VOI_NAME{v}),'_child']));
    end
   
%     % Plot Conceptual RDMs--------------------------------
%     num_plot_on_page=30;
%     num_pages=ceil(numel(conceptual_RDMs)/num_plot_on_page);
%     for page_index=1:num_pages
%         plot_start=(page_index-1)*num_plot_on_page+1;
%         plot_end=(page_index*num_plot_on_page);
%         if plot_end>numel(conceptual_RDMs)
%             plot_end=numel(conceptual_RDMs);
%         end
%         plot_range=plot_start:plot_end;
%         figureRDMs(concatenateRDMs(conceptual_RDMs(plot_range)), user_options_visualize_RDM,...
%             struct('fileName', ['conceptual_RDMs_',num2str(page_index)], 'figureNumber', 1));
%     end
end
%% (Step 2) Determine dendrograms for the clustering of the conditions
%%%%%%%%%%%%%WEIRD! The condition color and condition numbers do not match!
if(BOOLEAN_1ST_ANALYSYS)
    user_options_visualize_RDM.saveFiguresPDF=0;
    user_options_visualize_RDM.saveFigurePDF=0;
    user_options_visualize_RDM.rootPath=PATH_1ST_RESULT; %Save the results starting from here to the 1st order part
%     % (all/2nd grade/5th grade/child/adult) x  VOIs----------------
%     for group_index=1:numel(LIST_GROUP_NAME)
%         dendrogramConditions(list_neural_RDMs{group_index}, user_options_visualize_RDM,...
%             struct(...
%             'useAlternativeConditionLabels', true,...
%             'alternativeConditionLabels', {CONDITION_LABELS}));
%         % The pdf file name is set as default, so I need to modify them into
%         % desired names.
%         for v=1:numel(LIST_VOI_NAME)
%             old_name = [user_options_visualize_RDM.analysisName,'_Dendrogram_',...
%                 regexprep(char(LIST_VOI_NAME{v}),'_','-'),'.pdf'];
%             path_old_name = fullfile(user_options_visualize_RDM.rootPath,old_name);
%             new_name= [user_options_visualize_RDM.analysisName,'_Dendrogram_',char(LIST_VOI_NAME{v}),...
%                 '_',LIST_GROUP_NAME{group_index},'.pdf'];
%             path_new_name = fullfile(user_options_visualize_RDM.rootPath,new_name);
%             movefile(path_old_name,path_new_name);
%         end
%     end
    
    %% (Step 3) MDS of conditions
    % (all/2nd grade/5th grade/adult/child)x  VOIs----------------
    user_options_visualize_RDM.rootPath=PATH_1ST_RESULT; %Save the results starting from here to the 1st order part
    user_options_visualize_RDM.dotSize=8;
    user_options_visualize_RDM.saveFiguresPDF=0;
    user_options_visualize_RDM.saveFigurePDF=0;

    % Save the MDS results in the structure
    % list_MDS_results.(name_MDS_type).(name_group)=MDS_results;
    
    for MDS_criterion=1:numel(LIST_MDS_CRITERION) %Test different MDS criteria
        name_MDS_type=LIST_MDS_CRITERION{MDS_criterion};
        user_options_visualize_RDM.criterion=name_MDS_type;
        path_MDS_type=fullfile(user_options_visualize_RDM.rootPath,name_MDS_type);
        mkdir(path_MDS_type)
        
        for group_index=1:numel(LIST_GROUP_NAME)
            [MDS_results]=MDSConditionsReturnResults(list_neural_RDMs{group_index},...
                                user_options_visualize_RDM,...
                                struct('alternativeConditionLabels', {CONDITION_LABELS}));
            name_group=LIST_GROUP_NAME{group_index};
            
            %Save the MDS results in the structure           
            valid_name_group=matlab.lang.makeValidName(name_group);
            valid_name_MDS_type=matlab.lang.makeValidName(name_MDS_type);

            list_MDS_results.(valid_name_MDS_type).(valid_name_group)=MDS_results;
            
%             % The pdf file name is set as default, so I need to modify them into
%             % desired names.
%             for v=1:numel(LIST_VOI_NAME)
%                 name_VOI=char(LIST_VOI_NAME{v});
%                 old_name = [user_options_visualize_RDM.analysisName,'_FirstOrderMDS_',...
%                     regexprep(name_VOI,'_','-'),'.pdf'];
%                 path_old_name = fullfile(user_options_visualize_RDM.rootPath,old_name);
%                 new_name= [user_options_visualize_RDM.analysisName,'_FirstOrderMDS_',name_VOI,...
%                     '_',name_group,'.pdf'];
%           
%                 path_new_name = fullfile(path_MDS_type,new_name);                
%                 movefile(path_old_name,path_new_name);
%             end
        end
    end
end
%Save the options for visualization
%make sure the path to save is set properly
user_options_visualize_RDM.rootPath=PATH_1ST_RESULT; %Save the results starting from here to the 1st order part
mkdir(PATH_1ST_RESULT);
file_options_save=...
    fullfile(user_options_visualize_RDM.rootPath,['user_options_visualize_RDM_',VERSION_NUM,'_',date,'.mat']);

save(file_options_save,'user_options_visualize_RDM')

%Save the MDS results for visualization in R
file_MDS_result_save=...
    fullfile(user_options_visualize_RDM.rootPath,['MDS_Results_',VERSION_NUM,'_',date,'.mat']);
save(file_MDS_result_save,'list_MDS_results');
