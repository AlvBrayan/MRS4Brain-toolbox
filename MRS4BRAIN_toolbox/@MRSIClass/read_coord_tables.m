%%  read_coord_tables.m 

% Copyright All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, MRS4Brain research group @ CIBM MRI EPFL AIT, 2024
% See the LICENSE.TXT file for more details.

% Brayan Alves, CIBM - MRS4Brain group, 2022
% Guillaume Briand, CIBM - MRS4Brain group, 2023
% 
% USAGE : MRSI Class public method 
% msg = obj.read_coord_tables(quantify_path)
% 
% DESCRIPTION :
% Read tables and coord files to extract metabolites and fitting properties
%
% INPUTS :
% obj           = MRSI Class object with properties and methods
% quantify_path = lcmodel coordfiles folder directory
%
% OUTPUT :
% msg           = Error message
function msg = read_coord_tables(obj,quantify_path)
msg = {''};
N_met = 24;

if nargin < 2
    slice_directories = dir(fullfile(obj.results_folder,'Slice_N*'));
    slice_directories = {slice_directories.name};
%     slice_directories = {'Slice_N1'};
    quantify_paths = cell(length(slice_directories),1);
    qp = 1;
    for k = 1:length(slice_directories)
        if isfolder(fullfile(obj.results_folder,slice_directories{k},'Data','tables'))
            quantify_paths{qp} = fullfile(obj.results_folder,slice_directories{k});
            qp = qp + 1;
        end
    end
    qp = qp - 1;
    quantify_paths = quantify_paths(~cellfun('isempty',quantify_paths));
    if qp == 0
        msg = {'No existing quantified data','Please check your process and your quantifications'};
        return
    end
else
    qp = 1;
    quantify_paths{1} = quantify_path;
end

% Voxel matrix
Matrix_Size = obj.acq_params.matrix_sz;
Nx = Matrix_Size(1); Ny = Matrix_Size(2);

for pa = 1:qp
    [~,folder_name,~] = fileparts(quantify_paths{pa});
    slice = str2double(regexp(folder_name,'\d*','Match'));
    directory_data = fullfile(quantify_paths{pa},'Data');
    directory_table = fullfile(directory_data,'tables');
    directory_results = fullfile(quantify_paths{pa},'Results');
    if ~isfolder(directory_results)
        mkdir(directory_results)
    end
    src_dir = dir([directory_data '\*.RAW']);
    name_list = {src_dir.name};
    check_fullname = name_list{1};
    check_name_list = strsplit(check_fullname,'@');
    check_name = check_name_list{1};
    if ~contains(check_name,folder_name)
        folder_name = check_name;
    end
    % Intitialization, find N_met
    [x,y] = find(squeeze(obj.Final_mask(slice,:,:)) == 1);
    % folder_name = '1avg_1st'; % tune if not the same as the folder_name
    fid = -1;
    h = 1;
    while fid == -1
        filename = [folder_name '@' num2str(x(h)) '_' num2str(y(h)) '.table'];
        fid = fopen(fullfile(directory_table,filename));
        if fid~=-1
            [C, ~] = textscan(fid, '%f %s %f %s', 'headerlines', 5);
            if isempty(C{1,4})
                fid=-1;
            end
        end
        h = h + 1;
        if h>length(x)
            fid = 0;
        end
    end
    h=h-1;
    if isfile(fullfile(directory_table,filename))
        % Number of metabolite calculated
        fid = fopen(fullfile(directory_table,filename));
        [C, ~] = textscan(fid, '%f %s %f %s', 'headerlines', 5);
        N_met = length(C{1,1});
        name_conc_tot = C{1,4};
        fclose(fid);
    end
    coord_filename = [folder_name '@' num2str(x(h)) '_' num2str(y(h)) '.coord'];
    fid = fopen(fullfile(directory_data,coord_filename));
    if isfile(fullfile(directory_data,coord_filename))
        % Number of metabolite calculated
        [S, ~] = textscan(fid, '%s','Delimiter', '\n');
        fclose(fid);
        S = S{1,1};
        idx_ppmaxis = find(contains(S, 'points on ppm-axis = NY'));
        idx_SNR_FWHM = contains(S, 'S/N');
        idx_data_points = find(contains(S, 'NY phased data points follow'));
        idx_baseline_points = find(contains(S, 'NY background values follow'));
        diff_points = idx_data_points - idx_ppmaxis;
        S2 = cellfun(@str2double,cellfun(@strsplit, ...
            S(idx_ppmaxis+1:idx_ppmaxis+diff_points-1),'UniformOutput',false),'UniformOutput',false);
        data = cat(2,S2{:}).';
    end
    % Matrix initialization
    abs_conc = zeros(N_met,Nx,Ny);
    rel_conc = zeros(N_met,Nx,Ny);
    CRLB = zeros(N_met,Nx,Ny);
    % Merged results
    N_points = length(data);
    ppm_scale = data;
    data_pts_tot = zeros(N_points,Nx,Ny);
    fit_pts_tot = zeros(N_points,Nx,Ny);
    baseline_pts_tot = zeros(N_points,Nx,Ny);
    SNR = zeros(Nx,Ny);
    FWHM = zeros(Nx,Ny);
    for i = 1:Nx
        for j = 1:Ny
            if (obj.Final_mask(slice,i,j) == 1)
                disp([num2str(i), ' ' num2str(j)])
                filename = [folder_name '@' num2str(i) '_' num2str(j) '.table'];
                fid = fopen(fullfile(directory_table,filename));
                if isfile(fullfile(directory_table,filename))
                    % Getting the concentrations
                    [C, ~] = textscan(fid, '%f %s %f %s', 'headerlines', 5);
                    fclose(fid);
                    if length(C{1,3}) >= N_met
                        CRB_temp = cell(N_met);
                        for Met = 1:N_met
                            abs_conc(Met,i,j) = C{1,1}(Met);
                            rel_conc(Met,i,j) = C{1,3}(Met);
                            CRB_temp(Met) = C{1,2}(Met);
                        end
                        CRLB(:,i,j) = CRB_c(CRB_temp);
                    end
                end
                coord_filename = [folder_name '@' num2str(i) '_' num2str(j) '.coord'];
                fid = fopen(fullfile(directory_data,coord_filename));
                if isfile(fullfile(directory_data,coord_filename))
                    [S, ~] = textscan(fid, '%s','Delimiter', '\n');
                    fclose(fid);
                    S = S{1,1};
                    if(~isempty(find(contains(S, 'points on ppm-axis = NY'))))
                        S2 = cellfun(@str2double,cellfun(@strsplit, ...
                            S(idx_data_points+1:idx_baseline_points+diff_points-1),'UniformOutput',false),'UniformOutput',false);
                        % Getting the corrected data points from coord file
                        data_pts_tot(:,i,j) = cat(2,S2{1:diff_points-1}).';
                        % Getting the fitted data points from coord file
                        fit_pts_tot(:,i,j) = cat(2,S2{diff_points+1:2*diff_points-1}).';
                        % Getting the baseline points from coord file
                        baseline_pts_tot(:,i,j) = cat(2,S2{2*diff_points+1:end}).';
                        % Getting the SNR & FWHM
                        S2 = cellfun(@str2double,cellfun(@strsplit, ...
                            S(idx_SNR_FWHM),'UniformOutput',false),'UniformOutput',false);
                        SF = rmmissing(S2{:});
                        SNR(i,j) = SF(2);
                        FWHM(i,j) = SF(1);
                    else
                        data_pts_tot(:,i,j) = zeros(N_points,1);
                        % Getting the fitted data points from coord file
                        fit_pts_tot(:,i,j) = zeros(N_points,1);
                        % Getting the baseline points from coord file
                        baseline_pts_tot(:,i,j) = zeros(N_points,1);
                        % Getting the SNR & FWHM
                        SNR(i,j) = 0;
                        FWHM(i,j) = 0;
                    end
                end
            end
        end
    end
    met_map = struct;
    met_map.abs_conc = abs_conc;
    met_map.rel_conc = rel_conc;
    met_map.CRB = CRLB;
    met_map.met_names = string(name_conc_tot);
    obj.Met_names = string(name_conc_tot);
    met_map.ppm_scale = ppm_scale;
    met_map.data_pts_tot = data_pts_tot;
    met_map.fit_pts_tot = fit_pts_tot;
    met_map.baseline_pts_tot = baseline_pts_tot;
    met_map.SNR = SNR;
    met_map.FWHM = FWHM;
    if pa == 1
        obj.Final_met_map = met_map;
    else
        obj.Final_met_map(slice) = met_map;
    end
    save(fullfile(directory_results,'metab_map.mat'),"met_map")
end

%% Additional functions
function CRB_f = CRB_c(CRB_temp)
CRB_f = zeros(length(CRB_temp),1);
for i = 1:length(CRB_temp)
    a = CRB_temp{i};
    CRB_f(i) = str2double(a(1:end-1))/100;
end