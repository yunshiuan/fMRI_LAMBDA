function [matlabbatch] = DCM_BMS_specify_analysis( BMS_path,file_model_space,...
             inference_method,family_file,bma,verify)
% Specify the type of Bayesian Model Selection(BMS).
% function [ matlabbatch ] = DCM_BMS_specify_analysis( BMS_path,file_model_space,...
%              inference_method,family_file,bma,verify)
% Inputs:
% BMS_path - 
% The working directory of the BMS.
% file_model_space -
% A string which represents the full path of the model space.
% inference_method -
% The method of BMS analysis. Should be either 'RFX' of 'FFX'.
%
% Optional inputs:
% family_file -
% The mat file which defines the famaly which each model belongs to.
% (defaults to not performing family-wise analysis: use empty string '')
% bma -
% A boolean value which decides whether Bayesian Model Averagin is
% performed. (defaults to 0)
% verify - 
% A boolean value which decides whether the validity of the analysis is being
% verified automatically. (defaults to 1)
% ________________________________________________________________________
% Yun-Shiuan Chuang - March 5, 2018
% Contact: yunshiuan.chuang@gmail.com

%% Default Argument check----------------------------------------------------------
if nargin < 3
    help DCM_BMS_specify_analysis; error('Necessary parameters are missing');
end
if  isempty(BMS_path)|| isempty(file_model_space)|| isempty(inference_method)
    help DCM_BMS_specify_analysis; error('Necessary parameters are missing');
end
if ~exist('family_file','var'); family_file = ''; end
if ~exist('bma','var'); bma = 0; end
if ~exist('verify','var'); verify = 1; end
if bma
    matlabbatch{1}.spm.dcm.bms.inference.bma.bma_yes.bma_famwin = 'famwin'; 
else
    matlabbatch{1}.spm.dcm.bms.inference.bma.bma_no = bma;
end
matlabbatch{1}.spm.dcm.bms.inference.dir = {BMS_path};
matlabbatch{1}.spm.dcm.bms.inference.sess_dcm = cell(1, 0);
matlabbatch{1}.spm.dcm.bms.inference.model_sp = {file_model_space};
matlabbatch{1}.spm.dcm.bms.inference.load_f = {''};
matlabbatch{1}.spm.dcm.bms.inference.method = inference_method;
matlabbatch{1}.spm.dcm.bms.inference.family_level.family_file = {family_file};
matlabbatch{1}.spm.dcm.bms.inference.verify_id = verify;
end

