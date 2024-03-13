%%  read_table.m 

% Copyright All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, MRS4Brain research group @ CIBM MRI EPFL AIT, 2024
% See the LICENSE.TXT file for more details.

% Guillaume Briand, CIBM - MRS4Brain group, 2023
% 
% USAGE : DWS Class public method 
% out = obj.read_table()
% 
% DESCRIPTION :
% Read table files to extract metabolite concentration values as well as
% CRLB, SNR and FWHM
%
% INPUTS :
% obj       = DWS Class object with properties and methods
%
% OUTPUT :
% msg       = Error message
function msg = read_table(obj)
msg = {''};
N_met = 24;

idx_met = find([obj.DWS_struct.met_bool]);
folder_name = fullfile(obj.result_dir,obj.foldername,'quantified');

% Find the number of metabolite calculated by lcmodel
study_met = obj.DWS_struct(idx_met(1)).sum_processed_study;
filenameRAW = study_met.filename(1:end-4);
lcmodel_folder = fullfile(folder_name,filenameRAW);
table_filename = fullfile(lcmodel_folder,[filenameRAW,'.table']);
fid = fopen(table_filename);
if isfile(table_filename)
    % Number of metabolite calculated
    [C, ~] = textscan(fid, '%f %s %f %s', 'headerlines', 5);
    N_met = length(C{1,1});
    name_conc_tot = C{1,4};
    fclose(fid);
end

% Merged results
abs_conc_tot = zeros(N_met,length(idx_met));
rel_conc_tot = zeros(N_met,length(idx_met));
crlb_tot = zeros(N_met,length(idx_met));
for m = 1:length(idx_met)
    study_met = obj.DWS_struct(idx_met(m)).sum_processed_study;
    filenameRAW = study_met.filename(1:end-4);
    lcmodel_folder = fullfile(folder_name,filenameRAW);
    table_filename = fullfile(lcmodel_folder,[filenameRAW,'.table']);

    % Results initialisation
    abs_conc = zeros(N_met,1);
    rel_conc = zeros(N_met,1);
    CRLB = zeros(N_met,1);
    SNR = 0;
    FWHM = 0;

    fid = fopen(table_filename);
    if isfile(table_filename)
        % Getting the concentrations
        [C, ~] = textscan(fid, '%f %s %f %s', 'headerlines', 5);
        fclose(fid);
        if length(C{1,3}) >= N_met
            CRLB_temp = cell(N_met);
            for Met = 1:N_met
                abs_conc(Met) = C{1,1}(Met);
                rel_conc(Met) = C{1,3}(Met);
                CRLB_temp(Met) = C{1,2}(Met);
            end
            CRLB = CRLB_c(CRLB_temp);
            % Getting the SNR & FWHM
            fid = fopen(table_filename);
            [Cprime, ~] = textscan(fid, '%s', 'headerlines',5);
            table_int = cell2table(Cprime{1,1});
            index = find(ismember(table_int.Var1,'S/N') == 1) + 2;
            index_FWHM = find(ismember(table_int.Var1,'FWHM') == 1) + 2;
            SNR = str2double(Cprime{1}{index});
            FWHM = str2double(Cprime{1}{index_FWHM});
            fclose(fid);
        end
    end
    quantif_res = struct;
    quantif_res.abs_conc = abs_conc;
    abs_conc_tot(:,m) = abs_conc;
    quantif_res.rel_conc = rel_conc;
    rel_conc_tot(:,m) = rel_conc;
    quantif_res.CRLB = CRLB;
    crlb_tot(:,m) = CRLB;
    quantif_res.met_names = string(name_conc_tot);
    quantif_res.SNR = SNR;
    quantif_res.FWHM = FWHM;
    obj.DWS_struct(idx_met(m)).quantif_res = quantif_res;
end
    obj.Met_names = string(name_conc_tot);
    obj.abs_conc_tot = abs_conc_tot;
    obj.rel_conc_tot = rel_conc_tot;
    obj.crlb_tot = crlb_tot;
end

%% Additional functions
function CRB_f = CRLB_c(CRB_temp)
CRB_f = zeros(length(CRB_temp),1);
for i = 1:length(CRB_temp)
    a = CRB_temp{i};
    CRB_f(i) = str2double(a(1:end-1))/100;
end
end