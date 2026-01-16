%%  quantify_data.m 

% Copyright All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, MRS4Brain research group @ CIBM MRI EPFL AIT, 2024
% See the LICENSE.TXT file for more details.

% Brayan Alves, CIBM - MRS4Brain group, 2022
% Guillaume Briand, CIBM - MRS4Brain group, 2023
% 
% USAGE : MRSI Class public method 
% msg = obj.quantify_data(mrsi_params)
% 
% DESCRIPTION :
% Fit and Quantify the saved RAW files using LCModel quantification tool
% (a control file is created to use LCModel)
%
% INPUTS :
% obj           = MRSI Class object with properties and methods
% mrsi_params   = Quantification parameters
%
% OUTPUT :
% msg           = Error message
function msg = quantify_data(obj,mrsi_params)
% Function to quantify the saved raw file using LCModel quantification tool
% a control file is created to use LCModel
msg = {''};
acq_params = obj.acq_params;
obj.mrsi_params = mrsi_params;

control_path = fullfile(obj.results_folder,'ccfiles');
if ~isfolder(control_path)
    mkdir(control_path)
end

slice_directories = dir(fullfile(obj.results_folder,'Slice_N*'));
slice_directories = {slice_directories.name};

quantify_paths = cell(length(slice_directories),1);
qp = 1;
for k = 1:length(slice_directories)
    if isfolder(fullfile(obj.results_folder,slice_directories{k},'Data'))
        quantify_paths{qp} = fullfile(obj.results_folder,slice_directories{k});
        qp = qp + 1;
    end
end
qp = qp - 1;
if qp == 0
    msg = {'No existing file for quantify data','Please check the study path and name'};
    return
end
quantify_paths = quantify_paths(~cellfun('isempty',quantify_paths));

current_dir = pwd;
for pa = 1:qp
    % First the ordinary data
    MRSI_datadir = fullfile(quantify_paths{pa},'Data');
    if ~isfolder(fullfile(MRSI_datadir,'tables'))
        mkdir(fullfile(MRSI_datadir,'tables'));
    end

    src_dir = dir(fullfile(MRSI_datadir,'*.RAW'));
    name_list = {src_dir.name}; % Metabolite and reference files
    size_name = length(name_list)/2;

    parfor t = 1:size_name
        size(name_list);
        filename = name_list{2*t-1};
        filename = filename(1:end-4);
        BatchLCModel(filename,MRSI_datadir,mrsi_params,acq_params,control_path);
    end
    if any(~cellfun(@isempty,msg))
        return
    end
    
    % Quantify missed voxels if any during parfor
    for t = 1:size_name
        filename = name_list{2*t-1};
        filename = filename(1:end-4);
        if ~isfile(fullfile(MRSI_datadir,[filename '.ps'])) % Check if any voxels aren't quantified
            BatchLCModel(filename,MRSI_datadir,mrsi_params,acq_params,control_path);
        end
    end
end
cd(current_dir)
end

%% Additional functions

%%  BatchLCModel.m 
% Brayan Alves, CIBM - MRS4Brain group, 2022
% Guillaume Briand, CIBM - MRS4Brain group, 2023
% 
% USAGE : MRSI Class public method 
% msg = BatchLCModel(filename,MRSIdataDir,mrsi_params,acq_params,control_path)
% 
% DESCRIPTION :
% Quantification of MRSI data by running LCModel with a CONTROL file
%
% INPUTS :
% filename           = Name of the file to be fitted by LCModel
% MRSIdataDir        = Directory of the MRSI data (metabolite and reference)
% mrsi_params        = Quantification parameters
% acq_params         = Acquisition parameters
% control_path       = Directory of the control files
%
% OUTPUT :
% msg               = Error message
function msg = BatchLCModel(filename,MRSIdataDir,mrsi_params,acq_params,control_path)

% Controlfile
Cfile_name = create_control_file(filename,MRSIdataDir,mrsi_params,acq_params,control_path); 

% Run LCModel
if strcmp(computer,'MACI64')
    cd(mrsi_params.LCModel_path)
    [status, ~] = system(['./lcmodel <  "' Cfile_name '"']);
elseif strcmp(computer,'PCWIN64')
    cd(mrsi_params.LCModel_path)
    [status, ~] = system(['LCModel.exe < "' Cfile_name '"']);
elseif strcmp(computer,'GLNXA64')
    cd(mrsi_params.LCModel_path)
    [status, ~] = system(['/home/lcmodel/.lcmodel/bin/lcmodel < "' Cfile_name '"']);
end

if status ~= 0
    msg = {'Error with LCModel quantification, Please check your LCModel path'};
end
end

%%  create_control_file.m 
% Brayan Alves, CIBM - MRS4Brain group, 2022
% Guillaume Briand, CIBM - MRS4Brain group, 2023
% 
% USAGE : MRSI Class public method 
% cfile = create_control_file(filename,MRSIdataDir,mrsi_params,acq_params,cfile_dir)
% 
% DESCRIPTION :
% Create control file from quantification and acquistion parameters
%
% INPUTS :
% filename           = Name of the file to be fitted by LCModel
% MRSIdataDir        = Directory of the MRSI data (metabolite and reference)
% mrsi_params        = Quantification parameters
% acq_params         = Acquisition parameters
% cfile_dir          = Directory of the control files
%
% OUTPUT :
% cfile              = Control filename
function cfile = create_control_file(filename,MRSIdataDir,mrsi_params,acq_params,cfile_dir)

% Input initialization
basis_set = mrsi_params.Basis_set;
if nargin < 5
    cfile_dir = MRSIdataDir;
end

cfile = fullfile(cfile_dir,[filename,'.CONTROL']);
fileid = fopen(cfile,'w');
fprintf(fileid,' $LCMODL\n');
fprintf(fileid,' ATTH2O = 1.0\n'); % attenuation of the NMR-visible water signal
fprintf(fileid,'NRATIO = %i\n',mrsi_params.NRATIO); % NRATIO parameter : number of soft constraints on concentration ratios
% if ~mrsi_params.NRATIO
%     fprintf(fileid,'NRATIO = 0\n'); % NRATIO parameter : number of soft constraints on concentration ratios
% end
if ~mrsi_params.NSIMUL
    fprintf(fileid,'NSIMUL = 0\n'); % NSIMUL parameter : number of Basis Spectra that you will simulate (default = 13)
end
COMB = mrsi_params.CHCOMB;  % Combination of metabolites
fprintf(fileid,' NCOMBI = %i\n',length(COMB));
for j = 1:length(COMB)
    fprintf(fileid," CHCOMB(%i) = '%s'\n",j,COMB{j});
end
OMIT = mrsi_params.CHOMIT;
fprintf(fileid,' NOMIT = %i\n',length(OMIT));
for j = 1:length(OMIT)
    fprintf(fileid," CHOMIT(%i) = '%s'\n",j,OMIT{j});
end
USE1 = mrsi_params.CHUSE1; % Basic spectra in primary analysis
fprintf(fileid,' NUSE1 = %i\n',length(USE1));
for j = 1:length(USE1)
    fprintf(fileid," CHUSE1(%i) = '%s'\n",j,USE1{j});
end
fprintf(fileid,' CONREL = %.2f\n',mrsi_params.Rel_conc); % Relative metabolite concentration
fprintf(fileid,' DELTAT = %.2i\n',1./acq_params.spectralwidth); % 1/samplerate
fprintf(fileid,' DKNTMN = %1.2f\n',mrsi_params.DKNTMN); % DKNTMN
%  DOECC = F

fprintf(fileid,' DOWS = T\n'); % DO water scaling

%  DOREFS = T
fprintf(fileid,' FWHMBA = 0.0050\n');
fprintf(fileid,' HZPPPM = %.3f\n',acq_params.resfreq); % NMRfreq
fprintf(fileid,' LCOORD = 9\n'); % Save coord file, Y : LCOORD = 9, N : LCOORD = 0
fprintf(fileid,' LTABLE = 7\n'); % Save table file, Y : LTABLE = 7, N : LTABLE = 0
fprintf(fileid," NAMREL = '%s'\n",mrsi_params.Rel_met); % Relative metabolite name
fprintf(fileid,' NCALIB = 0\n');

fprintf(fileid,' NUNFIL = %i\n',acq_params.np_met); % FidPoints

fprintf(fileid," PGNORM = 'US'\n");
fprintf(fileid,' PPMEND = %1.2f\n',mrsi_params.PPMend); % right edge of the window
fprintf(fileid,' PPMST = %1.2f\n',mrsi_params.PPMstart); % left edge of the window
fprintf(fileid,' RFWHM = 1.8\n');
if mrsi_params.VITRO
    fprintf(fileid,' VITRO = T\n');
else
    fprintf(fileid,' VITRO = F\n');
end
fprintf(fileid,' WCONC = %1.2f\n',mrsi_params.WCONC); % the NMR-visible water concentration (mM) in the voxel
fprintf(fileid,' NEACH = 999\n');
fprintf(fileid,' SHIFMN = -0.2,-0.1\n');
fprintf(fileid,' SHIFMX = 0.3,0.3\n');
fprintf(fileid,' KEY = %i\n',mrsi_params.KEY); % Licence KEY
fprintf(fileid," OWNER = '%s'\n",mrsi_params.OWNER); % Licence OWNER
fprintf(fileid," FILBAS = '%s'\n",basis_set);
fprintf(fileid,' DEGPPM = %1.2f\n', mrsi_params.DEGPPM);
fprintf(fileid,' DEGZER = %1.2f\n', mrsi_params.DEGZER);
fprintf(fileid,' SDDEGP = %1.2f\n', mrsi_params.SDDEGP);
fprintf(fileid,' SDDEGZ = %1.2f\n', mrsi_params.SDDEGZ);
add_title = char(datetime);
fprintf(fileid," TITLE = '%s %s'\n",filename,add_title);
fprintf(fileid," FILPS = '%s'\n",fullfile(MRSIdataDir,[filename,'.ps']));
fprintf(fileid," FILCOO = '%s'\n",fullfile(MRSIdataDir,[filename,'.coord']));
fprintf(fileid," FILTAB = '%s'\n",fullfile(MRSIdataDir,'tables',[filename,'.table']));
fprintf(fileid," FILRAW = '%s'\n",fullfile(MRSIdataDir,[filename,'.RAW']));
fprintf(fileid," FILH2O = '%s'\n",fullfile(MRSIdataDir,[filename,'w.RAW']));
fprintf(fileid,' $END');

fclose(fileid);
end