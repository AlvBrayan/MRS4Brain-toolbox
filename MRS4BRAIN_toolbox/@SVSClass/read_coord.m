%% read_coord.m 

% Copyright All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, MRS4Brain research group @ CIBM MRI EPFL AIT, 2024
% See the LICENSE.TXT file for more details.

% Guillaume Briand, CIBM - MRS4Brain group, 2023
% 
% USAGE : SVS Class public method 
% out = obj.read_coord()
% 
% DESCRIPTION :
% Read coord files to extract LCModel fitting data points, fit points and
% baseline
%
% INPUTS :
% obj       = SVS Class object with properties and methods
%
% OUTPUT :
% msg       = Error message
function msg = read_coord(obj)
msg = {''};

idx_met = find([obj.SVS_struct.met_bool]);
folder_name = fullfile(obj.result_dir,obj.foldername,'quantified');

% Find the number of points calculated by lcmodel
study_met = obj.SVS_struct(idx_met(1)).sum_processed_study;
filenameRAW = study_met.filename(1:end-4);
lcmodel_folder = fullfile(folder_name,filenameRAW);
coord_filename = fullfile(lcmodel_folder,[filenameRAW,'.coord']);
assignin("base","coord_filename",coord_filename)
fid = fopen(coord_filename);
if isfile(coord_filename)
    % Number of metabolite calculated
    [S, ~] = textscan(fid, '%s','Delimiter', '\n');
    fclose(fid);
    S = S{1,1};
    idx_ppmaxis = find(contains(S, 'points on ppm-axis = NY'));
    idx_data_points = find(contains(S, 'NY phased data points follow'));
    idx_fit_points = find(contains(S, 'NY points of the fit to the data follow'));
    idx_baseline_points = find(contains(S, 'NY background values follow'));
    fid = fopen(coord_filename);
    diff_points = idx_data_points - idx_ppmaxis;
    C = textscan(fid, '%f', 10*diff_points, 'headerlines', idx_ppmaxis); % ppm scale from the COORD file
    data = C{1,1};
end

% Merged results
N_points = length(data);
ppm_scale = data;
data_pts_tot = zeros(N_points,length(idx_met));
fit_pts_tot = zeros(N_points,length(idx_met));
baseline_pts_tot = zeros(N_points,length(idx_met));
for m = 1:length(idx_met)
    study_met = obj.SVS_struct(idx_met(m)).sum_processed_study;
    filenameRAW = study_met.filename(1:end-4);
    lcmodel_folder = fullfile(folder_name,filenameRAW);
    coord_filename = fullfile(lcmodel_folder,[filenameRAW,'.coord']);

    % Results initialisation
    data_corr = zeros(N_points,1);
    fit_data = zeros(N_points,1);
    baseline_fit = zeros(N_points,1);

    fid = fopen(coord_filename);
    if isfile(coord_filename)
        % Getting the corrected data points from coord file
        [data, ~] = textscan(fid, '%f', 10*diff_points, 'headerlines', idx_data_points);
        data_corr = data{1,1};
        fclose(fid);
        % Getting the fitted data points from coord file
        fid = fopen(coord_filename);
        [data, ~] = textscan(fid, '%f', 10*diff_points, 'headerlines', idx_fit_points);
        fit_data = data{1,1};
        fclose(fid);
        % Getting the baseline points from coord file
        fid = fopen(coord_filename);
        [data, ~] = textscan(fid, '%f', 10*diff_points, 'headerlines', idx_baseline_points);
        baseline_fit = data{1,1};
        fclose(fid);
    end
    quantif_fit = struct;
    quantif_fit.ppm_pts = ppm_scale;
    quantif_fit.data_pts = data_corr;
    data_pts_tot(:,m) = data_corr;
    quantif_fit.fit_pts = fit_data;
    fit_pts_tot(:,m) = fit_data;
    quantif_fit.baseline_pts = baseline_fit;
    baseline_pts_tot(:,m) = baseline_fit;
    obj.SVS_struct(idx_met(m)).quantif_fit = quantif_fit;
end
    obj.ppm_scale = ppm_scale;
    obj.data_pts_tot = data_pts_tot;
    obj.fit_pts_tot = fit_pts_tot;
    obj.baseline_pts_tot = baseline_pts_tot;
end