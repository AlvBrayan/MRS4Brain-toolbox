%% DWS_quantification.m

% Copyright All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, MRS4Brain research group @ CIBM MRI EPFL AIT, 2024
% See the LICENSE.TXT file for more details.

% Guillaume Briand, CIBM - MRS4Brain group, 2023
% 
% USAGE : DWS Class public method 
% msg = obj.DWS_quantification()
% 
% DESCRIPTION :
% Fitting and quantification of data using LCModel software
%
% INPUTS :
% obj       = DWS Class object with properties and methods
% prog_dbox = MRS4Brain Toolbox progress dialog box
%
% OUTPUT :
% msg       = Error message
function msg = DWS_quantification(obj,prog_dbox)

folder_name = fullfile(obj.result_dir,obj.foldername,'quantified');
if ~exist(folder_name,"dir")
    mkdir(folder_name);
end

prog_dbox.Message = 'Create RAW files to quantify spectra';
prog_dbox.Value = 1/5;
msg = obj.data2RAW(); % Create RAW files to quantify spectra

idx_met = find([obj.DWS_struct.met_bool]);
idx_ref = find([obj.DWS_struct.ref_bool]);
if isempty(idx_ref)
    idx_ref = 0;
end
current_dir = pwd;
prog_dbox.Message = 'Perform parallel LCModel quantification';
prog_dbox.Value = 2/5;
dws_param = obj.DWS_param;
dws_struct = obj.DWS_struct;
parfor t = 1:length(idx_met)
    BatchLCModel(dws_param,dws_struct,folder_name,idx_met(t),idx_ref);
end
cd(current_dir)

prog_dbox.Message = 'Read table and coord files for metabolite concentration and LCModel fitting';
prog_dbox.Value = 3/5;
pause(0.5)
obj.read_table(); % Read the table file for absolute and relative concentrations, CRLB, SNR and FWHM
obj.read_coord(); % Read the coord file for ppm_range, corrected and fitted data and the baseline of the fit
prog_dbox.Message = 'Apply biological modelling for diffusion coefficients';
prog_dbox.Value = 4/5;
obj.DWS_fitting(); % Apply biological model fitting for diffusion coefficients

study = obj.DWS_struct(1).raw_study;
filename = [study.day '_' study.timeacq '_DWS_experiment.mat'];
filename = fullfile(obj.result_dir,obj.foldername,filename);
DWS_experiment = obj;
save(filename,'DWS_experiment')
end

%% BatchLCModel.m 
% Guillaume Briand, CIBM - MRS4Brain group, 2023
% 
% USAGE :
% [status, cmdout] = BatchLCModel(dws_param,dws_struct,folder_name,idx_met,idx_ref)
% 
% DESCRIPTION :
% Perform LCModel quantification on the selected metabolite fid spectrum
% with or without reference fid spectrum. A LCModel control file is created
% to give information to LCModel software.
%
% INPUTS :
% dws_param     = DWS parameters
% dws_struct    = DWS structure data
% folder_name   = Name of quantification folder
% idx_met       = Index of the metabolite fid in the obj.DWS_struct
% idx_ref       = Index of the reference  fid in the obj.DWS_struct
%
% OUTPUT :
% status    = Status of the LCModel software
% cmdout    = Output of the LCModel software
function [status, cmdout] = BatchLCModel(dws_param,dws_struct,folder_name,idx_met,idx_ref)

% Controlfile
Cfile_name = create_control_file(dws_param,dws_struct,folder_name,idx_met,idx_ref);

% Run LCModel
if strcmp(computer,'MACI64')
    cd(dws_param.LCModel_path)
    [status, cmdout] = system(['./lcmodel <  "' Cfile_name '"']);
elseif strcmp(computer,'PCWIN64')
    cd(dws_param.LCModel_path)
    [status, cmdout] = system(['LCModel.exe < "' Cfile_name '"']);
elseif strcmp(computer,'GLNXA64')
    cd(svs_param.LCModel_path)
    [status, ~] = system(['/home/lcmodel/.lcmodel/bin/lcmodel < "' Cfile_name '"']);
end

end

%% create_control_file.m 
% Guillaume Briand, CIBM - MRS4Brain group, 2023
% 
% USAGE :
% cfile = create_control_file(dws_param,dws_struct,folder_name,idx_met,idx_ref)
% 
% DESCRIPTION :
% LCModel control file is created to give information to LCModel software
% using the DWS parameters given by the user
%
% INPUTS :
% dws_param     = DWS parameters
% dws_struct    = DWS structure data
% folder_name   = Name of quantification folder
% idx_met       = Index of the metabolite fid in the obj.DWS_struct
% idx_ref       = Index of the reference  fid in the obj.DWS_struct
%
% OUTPUT :
% cfile     = Name of the created control file (dir + filename)
function cfile = create_control_file(dws_param,dws_struct,folder_name,idx_met,idx_ref)
% Input initialization
study_met = dws_struct(idx_met).sum_processed_study;
filenameRAW = study_met.filename(1:end-4);
if idx_ref ~= 0
    study_wat = dws_struct(idx_ref).sum_processed_study;
    filenameH20 = study_wat.filename(1:end-4);
end
LCModel_results_folder = fullfile(folder_name,filenameRAW);
if ~exist(LCModel_results_folder,"dir")
    mkdir(LCModel_results_folder)
end
basis_set = dws_param.Basis_set;
cfile = fullfile(LCModel_results_folder,[filenameRAW,'.CONTROL']);
fileid = fopen(cfile,'w');
fprintf(fileid,' $LCMODL\n');
fprintf(fileid,' ATTH2O = 1.0\n'); % attenuation of the NMR-visible water signal
fprintf(fileid,'NRATIO = %.2f\n',dws_param.NRATIO); % NRATIO parameter : number of soft constraints on concentration ratios (default = 12)
% if ~dws_param.NRATIO
%     fprintf(fileid,'NRATIO = 0\n'); % NRATIO parameter : number of soft constraints on concentration ratios (default = 12)
% end
if ~dws_param.NSIMUL
    fprintf(fileid,'NSIMUL = 0\n'); % NSIMUL parameter : number of Basis Spectra that you will simulate (default = 13)
end
COMB = dws_param.CHCOMB;  % Combination of metabolites
fprintf(fileid,' NCOMBI = %i\n',length(COMB));
for j = 1:length(COMB)
    fprintf(fileid," CHCOMB(%i) = '%s'\n",j,COMB{j});
end
OMIT = dws_param.CHOMIT;
fprintf(fileid,' NOMIT = %i\n',length(OMIT));
for j = 1:length(OMIT)
    fprintf(fileid," CHOMIT(%i) = '%s'\n",j,OMIT{j});
end
USE1 = dws_param.CHUSE1; % Basic spectra in primary analysis
fprintf(fileid,' NUSE1 = %i\n',length(USE1));
for j = 1:length(USE1)
    fprintf(fileid," CHUSE1(%i) = '%s'\n",j,USE1{j});
end
fprintf(fileid,' CONREL = %.2f\n',dws_param.Rel_conc); % Relative metabolite concentration
fprintf(fileid,' DELTAT = %.2i\n',1./study_met.spectralwidth); % 1/samplerate
fprintf(fileid,' DKNTMN = %1.2f\n',dws_param.DKNTMN); % DKNTMN
fprintf(fileid,' DOECC = F\n'); % DO Eddy current correction
if idx_ref ~= 0
    fprintf(fileid,' DOWS = T\n'); % DO water scaling
else
    fprintf(fileid,' DOWS = F\n'); % DO water scaling
end
fprintf(fileid,' DOREFS = F\n'); % DO Cr referencing
fprintf(fileid,' FWHMBA = 0.0050\n');
fprintf(fileid,' HZPPPM = %.3f\n',study_met.resfreq); % NMRfreq
fprintf(fileid,' LCOORD = 9\n'); % Save coord file, Y : LCOORD = 9, N : LCOORD = 0
fprintf(fileid,' LTABLE = 7\n'); % Save table file, Y : LTABLE = 7, N : LTABLE = 0
fprintf(fileid," NAMREL = '%s'\n",dws_param.Rel_met); % Relative metabolite name
fprintf(fileid,' NCALIB = 0\n');
fprintf(fileid,' NUNFIL = %i\n',2048); % FidPoints
fprintf(fileid," PGNORM = 'US'\n");
fprintf(fileid,' PPMEND = %1.2f\n',dws_param.PPMend); % right edge of the window
fprintf(fileid,' PPMST = %1.2f\n',dws_param.PPMstart); % left edge of the window
fprintf(fileid,' RFWHM = 1.8\n');
if dws_param.VITRO
    fprintf(fileid,' VITRO = T\n');
else
    fprintf(fileid,' VITRO = F\n');
end
fprintf(fileid,' WCONC = %1.2f\n',dws_param.WCONC); % the NMR-visible water concentration (mM) in the voxel
fprintf(fileid,' NEACH = 999\n');
fprintf(fileid,' SHIFMN = -0.2,-0.1\n');
fprintf(fileid,' SHIFMX = 0.3,0.3\n');
fprintf(fileid,' KEY = %i\n',dws_param.KEY); % Licence KEY
fprintf(fileid," OWNER = '%s'\n",dws_param.OWNER); % Licence OWNER
fprintf(fileid," FILBAS = '%s'\n",basis_set);
fprintf(fileid,' DEGPPM = 0\n');
fprintf(fileid,' DEGZER = 0\n');
fprintf(fileid,' SDDEGP = 10\n');
fprintf(fileid,' SDDEGZ = 10\n');
add_title = char(datetime);
fprintf(fileid," TITLE = '%s %s'\n",filenameRAW,add_title);
fprintf(fileid," FILPS = '%s'\n",fullfile(LCModel_results_folder,[filenameRAW,'.ps']));
fprintf(fileid," FILCOO = '%s'\n",fullfile(LCModel_results_folder,[filenameRAW,'.coord']));
fprintf(fileid," FILTAB = '%s'\n",fullfile(LCModel_results_folder,[filenameRAW,'.table']));
fprintf(fileid," FILRAW = '%s'\n",fullfile(folder_name,'raw_lcmodel',[filenameRAW,'.RAW']));
if idx_ref ~= 0
    fprintf(fileid," FILH2O = '%s'\n",fullfile(folder_name,'raw_lcmodel',[filenameH20,'.RAW']));
end
fprintf(fileid,' $END');

fclose(fileid);
end