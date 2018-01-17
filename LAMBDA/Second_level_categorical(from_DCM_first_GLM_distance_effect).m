%% Part 0-a-i - Second Level (Categorical): Collapse all age and set age as a covariate in the second levle GLM
%% NOTE1: Since the first level GLM is for DCM (distance effect), so notation is ignore.
%That is, all trials are only categorized into 3 distance condition.
%% NOTE2: This second level GLM is for reviewing only. This help to check
% if the first level GLMs of DCM is correct.
addpath(char("D:\GoogleDrive\Lambda_code\m_file\LAMBDA\tool_code"));% enable read_mixed_csv()
path='D:\Yun-Shiuan_LAMBDA\DCM\Part0_First_GLM_for_DCM'; % To get the subject list with valid runs (see below)
cd(path);
subject_list=cellstr(ls('*')); % subject list (which has valid first level contrast, i.e., not being excluded)
subject_list=subject_list(~cellfun(@isempty,regexp(subject_list,'(?<=df)\d+','match')));% include only subject name

%% Read in Gender and Age as covariates
demographic=cellfun(@(x) regexprep(x,'"',''),...
                read_mixed_csv("D:\Yun-Shiuan_LAMBDA\demographic\tidy_demographic.csv",','),'un',0);
table_demopraphic=cell2table(demographic(2:end,2:end),...
                         'VariableNames',cellstr(demographic(1,2:end)));
table_demopraphic=sortrows(table_demopraphic,'sub_id','ascend'); % Need to sort bi subID to make sure the order is ascending
valid_Gender_Dummy=str2double(table2array(table_demopraphic(ismember(lower(table_demopraphic.sub_id),subject_list),'gender_dummy')));
valid_Age_Precise=str2double(table2array(table_demopraphic(ismember(lower(table_demopraphic.sub_id),subject_list),'scan_age_precise')));

%% Second Level Matrix: T_Parameter(from first level)~ B0 + Age + Gender
%Note that T contrast are named from 002~014 (001 is an F contrast)
T_Parameters_name={'B0_Near' 'B0_Medium' 'B0_Far'... 
                   't_M-N' 't_N-M' ...
                   't_F-N' 't_N-F' ...
                   't_F-M' 't_M-F'
                   };
numeric_Contrast_Index=pad(strtrim(cellstr((num2str([2:12]')))),4,'left','0');              
explicit_Mask={'D:\Yun-Shiuan_LAMBDA\template\for_this_study\binary_overlapped_WM_GM_common_mask_n40_with_signal.nii,1'};

for con=1:size(T_Parameters_name,2) %NOTE: con=1 refers to spmT_0002
    clear matlabbatch
    second_Result_Dir = strcat('D:\Yun-Shiuan_LAMBDA\DCM\Second_GLM_Reivew_Only\',T_Parameters_name{con});
    mkdir(second_Result_Dir);
    matlabbatch{1}.spm.stats.factorial_design.dir = {second_Result_Dir};
    
    first_T_Contrast=strcat('D:\Yun-Shiuan_LAMBDA\DCM\First_GLM_for_DCM\',...
                             subject_list,'\spmT_',numeric_Contrast_Index(con),'.nii,1');
    matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = first_T_Contrast;
    % Covariate
    % Age
    matlabbatch{1}.spm.stats.factorial_design.cov(1).c = valid_Age_Precise;
    matlabbatch{1}.spm.stats.factorial_design.cov(1).cname = 'Age';
    matlabbatch{1}.spm.stats.factorial_design.cov(1).iCFI = 1;
    matlabbatch{1}.spm.stats.factorial_design.cov(1).iCC = 1;
    % Gender
    matlabbatch{1}.spm.stats.factorial_design.cov(2).c = valid_Gender_Dummy;
    matlabbatch{1}.spm.stats.factorial_design.cov(2).cname = 'Gender';
    matlabbatch{1}.spm.stats.factorial_design.cov(2).iCFI = 1;
    matlabbatch{1}.spm.stats.factorial_design.cov(2).iCC = 1;
    % Other parameters
    matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.em = explicit_Mask;
    matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
    
    % Estimation------------------------------
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
    % Contrast Manager------------------------
    matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = ['pos_' T_Parameters_name{con}];
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = ['neg_' T_Parameters_name{con}];
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = -1;
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
    
    matlabbatch{3}.spm.stats.con.consess{3}.tcon.name = 'pos_Age';
    matlabbatch{3}.spm.stats.con.consess{3}.tcon.weights = [0 1];
    matlabbatch{3}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
    
    matlabbatch{3}.spm.stats.con.consess{4}.tcon.name = 'neg_Age';
    matlabbatch{3}.spm.stats.con.consess{4}.tcon.weights = [0 -1];
    matlabbatch{3}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
    
    matlabbatch{3}.spm.stats.con.consess{5}.tcon.name = 'pos_Gender';
    matlabbatch{3}.spm.stats.con.consess{5}.tcon.weights = [0 0 1];
    matlabbatch{3}.spm.stats.con.consess{5}.tcon.sessrep = 'none';
    
    matlabbatch{3}.spm.stats.con.consess{6}.tcon.name = 'neg_Gender';
    matlabbatch{3}.spm.stats.con.consess{6}.tcon.weights = [0 0 -1];
    matlabbatch{3}.spm.stats.con.consess{6}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.delete = 1;
    save(strcat(second_Result_Dir,'\second_glm_batch.mat'),'matlabbatch')
end

%% Run Second Level
for con=1:size(T_Parameters_name,2) 
    second_Result_Dir = strcat('D:\Yun-Shiuan_LAMBDA\DCM\Second_GLM_Reivew_Only\',T_Parameters_name{con});
    load(strcat(second_Result_Dir,'\second_glm_batch.mat'));
    spm_jobman('run',matlabbatch);
end

%% Output Result Table
for first_con=1:size(T_Parameters_name,2) 
    matlabbatch='';
    spm_File=strcat('D:\Yun-Shiuan_LAMBDA\DCM\Second_GLM_Reivew_Only\',T_Parameters_name{first_con},'\SPM.mat');
    matlabbatch{1}.spm.stats.results.spmmat = {spm_File};

    for second_con=1:6
        threshold=(second_con<=2)*(0.005)+(second_con>=3)*(0.05); 
        % Set a higher threshold for age and gender, 
        % as they're exploratoty.
        matlabbatch{1}.spm.stats.results.conspec(second_con).titlestr = '';
        matlabbatch{1}.spm.stats.results.conspec(second_con).contrasts = second_con;
        matlabbatch{1}.spm.stats.results.conspec(second_con).threshdesc = 'none';
        matlabbatch{1}.spm.stats.results.conspec(second_con).thresh = threshold;
        matlabbatch{1}.spm.stats.results.conspec(second_con).extent = 15;
        matlabbatch{1}.spm.stats.results.conspec(second_con).conjunction = 1;
        matlabbatch{1}.spm.stats.results.conspec(second_con).mask.none = 1;
    end
    
    matlabbatch{1}.spm.stats.results.units = 1;
    matlabbatch{1}.spm.stats.results.export{1}.csv = true;

    spm_jobman('run',matlabbatch);
    
    
    % Rename and move the csv to result_table folder
    result_table_Dir=strcat('D:\Yun-Shiuan_LAMBDA\DCM\Second_GLM_Reivew_Only\',T_Parameters_name{first_con},'\Result_table');
    mkdir(result_table_Dir);
    csv_Old_Name=ls(strcat('D:\Yun-Shiuan_LAMBDA\DCM\Second_GLM_Reivew_Only\',T_Parameters_name{first_con},'\spm*.csv'));
    csv_New_Name={['pos_' T_Parameters_name{first_con},'_p005_k15'] ['neg_' T_Parameters_name{first_con},'_p005_k15'],...
                  'pos_Age_p05_k15' 'neg_Age_p05_k15',...
                  'pos_Gender_p05_k15' 'neg_Gender_p05_k15'};    
    csv_Old_Full_Name=strcat('D:\Yun-Shiuan_LAMBDA\DCM\Second_GLM_Reivew_Only\',T_Parameters_name{first_con},'\',csv_Old_Name);
    csv_New_Full_Name=strcat(result_table_Dir,'\',csv_New_Name','.csv');

    for csv=1:6
        movefile(csv_Old_Full_Name(csv,:),char(csv_New_Full_Name(csv,:)));
    end
end

