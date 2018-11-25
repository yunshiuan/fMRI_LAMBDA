function [subj] = DCM_BMS_create_model_space(DCM_file_paths)
% Create a model space for Bayesian Model Selection (BMS) based on the
% estimated DCMs results.
% function model_space] = DCM_BMS_create_model_space(file_paths)
%
% inputs:

% DCM_file_paths -
% A cell of full path of each DCM file.
% Each path refers to an estimated DCMs of every subject x
% run x model.
% That is, the "file_paths(sub_id,run_index,model_index)" refers to 
% the estimated DCM file of the specicif model of a sub_id's run_num.
% Note that  all subjects should have same number of
% runs and models for the sake of folloing BMS 
% (this is the way SPM's DCM implements BMS).
% 
%outputs:
% subj-
% The model_space is the input for subsequent BMS.
%
% Example:
% [model_space] = DCM_BMS_create_model_space(file_paths)
%________________________________________________________________________
% Yun-Shiuan Chuang - March 5, 2018
% Contact: yunshiuan.chuang@gmail.com

% Constants---------------------------------------------
% Numbder of subjects
num_sub = size(DCM_file_paths,1);
% Number of runs
num_run = size(DCM_file_paths,2);
% Numder of models
num_model = size(DCM_file_paths,3);

% Initialize the empty model_space-----------------------
subj(1,num_sub)=struct();

% Insert information into the model_space----------------------------
for sub_id=1:num_sub

    for run_index=1:num_run
        
        for model_index=1:num_model
            file_DCM=DCM_file_paths{sub_id,run_index,model_index};
            load(file_DCM);
            % fname: The whole path of the estimated DCM file (without the '.mat' suffix)
            subj(sub_id).sess(run_index).model(model_index).fname = ...
                char(regexprep(DCM_file_paths(sub_id,run_index,model_index),'.mat$',''));
            % F:  free energy approximation to the model evidenct
            subj(sub_id).sess(run_index).model(model_index).F = DCM.F;
            % Ep: posterior mean parameters
            subj(sub_id).sess(run_index).model(model_index).Ep = DCM.Ep;
            % Cp: posterior covariance of parameters
            subj(sub_id).sess(run_index).model(model_index).Cp = DCM.Cp;
            fprintf(['model : ',num2str(model_index),...
                     '; run : ',num2str(run_index),...
                     '; subject : ',num2str(sub_id),'\n'])
        end
    end
end
end

