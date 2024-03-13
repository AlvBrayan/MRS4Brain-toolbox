%%  check_prexisting_results.m 

% Copyright All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, MRS4Brain research group @ CIBM MRI EPFL AIT, 2024
% See the LICENSE.TXT file for more details.

% Guillaume Briand, CIBM - MRS4Brain group, 2023
% 
% USAGE : MRSI Class public method 
% msg = obj.check_prexisting_results(brain_mask,save_figure)
% 
% DESCRIPTION :
% Check if there are already quantified data or RAW files in the studied
% result folder
%
% INPUTS :
% obj               = MRSI Class object with properties and methods
% MRI_brain_mask    = MRI Brain mask computed by registration with atlas
%
% OUTPUT :
% msg               = Error message
% processed         = Status of processing of data (0:not processed, 1:quantified, 2:RAW files but not quantified)
function [msg,processed] = check_prexisting_results(obj,MRI_brain_mask)
% check_prexisting_results allows the program to check if there are already
% quantifed data or RAW files in the studied result folder
msg = {''};

MatSize = obj.acq_params.matrix_sz;
% Check the number of slices that have been processed
slice_directories = dir(fullfile(obj.results_folder,'Slice_N*'));
slice_directories = {slice_directories.name};

quantify_paths = cell(length(slice_directories),1);
qp = 1;
for k = 1:length(slice_directories)
    % Check if table folder is filled with fitted and quantified data
    if isfolder(fullfile(obj.results_folder,slice_directories{k},'Data','tables'))
        if length(dir(fullfile(obj.results_folder,slice_directories{k},'Data','tables'))) > 2
            quantify_paths{qp} = fullfile(obj.results_folder,slice_directories{k});
            qp = qp + 1;
        end
    end
end
qp = qp - 1;
quantify_paths = quantify_paths(~cellfun('isempty',quantify_paths));

% Case of no quantified data
if or(qp < length(slice_directories),qp == 0)
    qp2 = 1;
    for k = 1:length(slice_directories)
        data_folder = fullfile(obj.results_folder,slice_directories{k},'Data');
        if isfolder(data_folder)
            if ~isempty(dir(fullfile(data_folder,'*.RAW')))
                quantify_paths{qp2} = fullfile(obj.results_folder,slice_directories{k});
                qp2 = qp2 + 1;
            end
        end
    end
    qp2 = qp2 - 1;
    if or(qp2 < length(slice_directories),qp2 == 0)
        processed = 0;
    % RAW files exist but not fitted/quantified
    else
        % brain_mask_filename = ['brain_mask_' num2str(MatSize(1)) '_' num2str(MatSize(2)) '.mat'];
        % if isfile(fullfile(obj.results_folder,'Registration',brain_mask_filename))
        %     final_mask_temp = load(fullfile(obj.results_folder,'Registration',brain_mask_filename));
        %     obj.Final_mask = final_mask_temp.Brain_mask;
        % else
            msg = obj.Brain_mask_comp(MRI_brain_mask,false); % call the Brain_mask_comp function
        % end
        processed = 2;
    end
% Case of fitted/quantified data
else
    processed = 1;
    obj.Final_met_map = struct;
    % Brain mask
    % brain_mask_filename = ['brain_mask_' num2str(MatSize(1)) '_' num2str(MatSize(2)) '.mat'];
    % if isfile(fullfile(obj.results_folder,'Registration',brain_mask_filename))
    %     final_mask_temp = load(fullfile(obj.results_folder,'Registration',brain_mask_filename));
    %     obj.Final_mask = final_mask_temp.Brain_mask;
    % else
        msg = obj.Brain_mask_comp(MRI_brain_mask,false); % call the Brain_mask_comp function
    % end
    % Fill the Final_met_map structure with the metabolite concentrations
    for pa = 1:qp
        met_map_filename = fullfile(quantify_paths{pa},'Results','metab_map.mat');
        % If the strcture already exist
        if isfile(met_map_filename)
            met_map_temp = load(met_map_filename); % Loading
            % Check if the there are the absolute and relative concentrations
            if and(isfield(met_map_temp.met_map,'rel_conc'), ...
                isfield(met_map_temp.met_map,'abs_conc'))
                slice = str2double(regexp(quantify_paths{pa},'\d*','Match'));
                if pa > 1
                    obj.Final_met_map(slice) = met_map_temp.met_map;
                else
                    obj.Final_met_map = met_map_temp.met_map;
                    obj.Met_names = obj.Final_met_map.met_names;
                end
            else
                obj.read_coord_tables(quantify_paths{pa})
            end
        % If the structure do not exist
        else
            obj.read_coord_tables(quantify_paths{pa})
        end
    end
end
end