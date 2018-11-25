%% Part 3 - iv-b: RSA : Prepare data for RSA usage (before construct RDM)
% trial 5 in the README.txt.
% Part3: 2nd order analysis.
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

% b: Split neural RDM by age: 2nd grader and 5th grader. 
% Run 1st and 2nd order RSA on each age group separately.

% [Run info has already been averaged out. Therefore, the number of runs no longer matters.]
addpath(genpath('D:\Yun-Shiuan_LAMBDA\rsatoolbox\'));

%% Set up
% Constants
% path
path_1st_result='D:\Yun-Shiuan_LAMBDA\RSA\trial_6_id_18_split_age\Part2_1st_oder_analysis';
path_2nd_result='D:\Yun-Shiuan_LAMBDA\RSA\trial_6_id_18_split_age\Part3_2nd_order_analysis';
% file
file_RDM_result='D:\Yun-Shiuan_LAMBDA\RSA\trial_6_id_18_split_age\Part1_neural_and_conceptual_RDMs\RDMs\RDMs_and_options_neural_and_conceptual_trial_6_09-Mar-2018.mat';
file_1st_user_option=fullfile(path_1st_result,'user_options_visualize_RDM_6_09-Mar-2018.mat');
file_demographic='D:\Yun-Shiuan_LAMBDA\demographic\tidy_demographic.csv';
file_valid_run='D:\Yun-Shiuan_LAMBDA\Run_inclusion_info\inclusive_runs_indexes.csv';

% Parameters
list_VOI_name={"R_V1","L_V1",...
               "R_IPS","L_IPS"}';
amount_each_format=6;
version_num='6';
version_date= ['RSA_trial_',version_num,'_',date];
list_group_name={'all_sub','2_grade','5_grade'};
% Load in needed data
load(file_RDM_result,...
     'neural_RDMs',...
     'neural_RDMs_avg_2_grade','neural_RDMs_avg_5_grade',...
     'neural_RDMs_avg_all','conceptual_RDMs',...
     'user_options_RDM_construction')
load(file_1st_user_option,'user_options_visualize_RDM');
list_neural_RDMs={neural_RDMs_avg_all,neural_RDMs_avg_2_grade,neural_RDMs_avg_5_grade};

% Options for2nd order analysis
clear user_options_2nd_order
% Edit the color bar color due to Matlab version issue
%(R2015b change default color maps. Should now be colormap(jet(256)))
% user_options_visualize_RDM.colourScheme=colormap(jet(256)); [Already edify the source codes]
user_options_2nd_order=user_options_visualize_RDM;

user_options_2nd_order.analysisName = version_date;
user_options_2nd_order.rootPath=path_2nd_result;                   
user_options_2nd_order.RDMcorrelationType='Kendall_taua';
user_options_2nd_order.distanceMeasure='Kendall_taua';
reduced_conceptual_RDMs=conceptual_RDMs(1,...
    [1,1+3,1+5,1+7,...
    11,11+3,11+5,11+7,...
    21,22,23,...
    23+3,23+5,23+7,...
    33,33+3,33+5,33+7]);

% Recolor the conceptual models
model_colors=cell(numel(conceptual_RDMs));
for model_index=1:numel(conceptual_RDMs)
    model_name=conceptual_RDMs(model_index).name;
    if ~isempty(regexp(model_name,'FCF','match'))
        model_colors{model_index}=[1,0.7,1];
    elseif ~isempty(regexp(model_name,'GCF','match'))
        model_colors{model_index}=[0,1,0];
    elseif ~isempty(regexp(model_name,'NCF','match'))
        model_colors{model_index}=[0,0,1];
    else
        model_colors{model_index}=[1,0.7,0];
    end
end

for model_index=1:numel(conceptual_RDMs)  
    conceptual_RDMs(model_index).color=...
        model_colors{model_index};
end
% Read in run inclusion index info (note that only subject matters for
% t-map based analysis)
run_inclusion_index=cellfun(@(x) regexprep(x,'"',''),...
    read_mixed_csv(file_valid_run,','),'un',0);
run_inclusion_index=table(run_inclusion_index(2:end,2),run_inclusion_index(2:end,3),...
    'VariableNames',{'sub_id','run_num'});
% Derive subjects with valid runs
subject_list=unique(run_inclusion_index.sub_id);
% Read in demographic data (for age group split)
table_demographic=read_mixed_csv_to_table(file_demographic);
table_demographic=table_demographic(:,{'sub_id','grade'});
table_demographic.sub_id=lower(table_demographic.sub_id);
table_demographic=join(table(subject_list,'VariableNames',{'sub_id'}),...
    table_demographic);
index_all=true(numel(table_demographic.grade),1);
index_2_grade=strcmp(table_demographic.grade,'2');
index_5_grade=strcmp(table_demographic.grade,'5');
list_group_index={index_all,index_2_grade,index_5_grade};
%% (Step 1) Visualizatopm: 2nd order simmilarity matrix and 2nd order MDS
% Display a 2nd order simmilarity matrix between neural RDMs and conceptual
% RDMs             
for group_index=1:numel(list_neural_RDMs)
    pairwiseCorrelateRDMs({list_neural_RDMs{group_index}, conceptual_RDMs},...
        user_options_2nd_order);
    % The pdf file name is et as default, so I need to modify them into
    % desired names.
    old_name = [user_options_visualize_RDM.analysisName,'_secondOrderSM.pdf'];
    path_old_name = fullfile(user_options_2nd_order.rootPath,old_name);
    new_name= [user_options_visualize_RDM.analysisName,'_secondOrderSM_',...
        list_group_name{group_index},'.pdf'];
    path_new_name = fullfile(user_options_2nd_order.rootPath,new_name);
    movefile(path_old_name,path_new_name);
end

% MDS: Plot all RDMs on a MDS plot to visualise pairwise
% distances.-----------------------------------------------------------
% Reduce the conceptual RDMs to be plotted.
% Without reducing: Points in the configuration have co-located.  Try a different
% starting point, or use a different criterion.
for group_index=1:numel(list_neural_RDMs)
    MDSRDMs({list_neural_RDMs{group_index}, reduced_conceptual_RDMs}, user_options_2nd_order,...
        struct('titleString', 'MDS of neural and conceptual'));
    % The pdf file name is et as default, so I need to modify them into
    % desired names.
    old_name = [user_options_visualize_RDM.analysisName,'_secondOrderMDS_.pdf'];
    path_old_name = fullfile(user_options_2nd_order.rootPath,old_name);
    new_name= [user_options_visualize_RDM.analysisName,'_secondOrderMDS_',...
        list_group_name{group_index},'.pdf'];
    path_new_name = fullfile(user_options_2nd_order.rootPath,new_name);
    movefile(path_old_name,path_new_name);
end
% If one want to save manually instead of using the automatic saving function of the
% MDSRDMs()
% saveas(gcf,fullfile(path_2nd_result,[user_options_2nd_order.analysisName,'_SecondOrderMDS_.pdf']))


%% (Step 2)stats_p_r = compareRefRDM2candRDMs(refRDM, candRDMs[, userOptions]) - statistical inference:
clear models
% test the relatedness and compare the candidate RDMs
user_options_2nd_order.RDMcorrelationType='Kendall_taua';
user_options_2nd_order.RDMrelatednessTest = 'subjectRFXsignedRank';
user_options_2nd_order.RDMrelatednessThreshold = 0.05;
user_options_2nd_order.figureIndex = [1:2];
user_options_2nd_order.RDMrelatednessMultipleTesting = 'FDR';
user_options_2nd_order.candRDMdifferencesTest = 'subjectRFXsignedRank';
user_options_2nd_order.candRDMdifferencesThreshold = 0.05;
user_options_2nd_order.candRDMdifferencesMultipleTesting = 'none';
%Remain to be
%debugged: Undefined function 'mod' for input arguments of type 'matlab.ui.Figure'.----------------------------------------------------------------
% R_V1
for group_index=1:numel(list_neural_RDMs)
    % The reference model (the neural RDMs from every subjects)
    neural_RDMs_avg=list_neural_RDMs{group_index};
    % Candidate Models: Neural RDMS (averaged across subjects)
    for conceptual_index=1:numel(reduced_conceptual_RDMs)
        models{conceptual_index}=reduced_conceptual_RDMs(conceptual_index);
    end
    % Candidate Models: Neural RDMS (averaged across subjects)
    for neural_index=1:numel(neural_RDMs_avg)
        models{conceptual_index+neural_index}=neural_RDMs_avg(neural_index);
    end
    for v=1:numel(list_VOI_name)
        % Before the bugs are fixed, use try-catch to prevent the
        % program from collapsing.
        try
            stats_p_r=compareRefRDM2candRDMs(neural_RDMs(v,list_group_index{group_index}), models, user_options_2nd_order);
        catch ME
            ME.message
        end
        file_figure = fullfile(path_2nd_result,...
            [version_date,'_',char(list_VOI_name{v}),'_',list_group_name{group_index},...
            '_relatedness_test.pdf']);
        print(file_figure,'-dpdf','-fillpage');
        
    end
end
