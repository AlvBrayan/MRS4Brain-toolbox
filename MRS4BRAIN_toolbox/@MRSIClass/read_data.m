%%  read_data.m 

% Copyright All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, MRS4Brain research group @ CIBM MRI EPFL AIT, 2024
% See the LICENSE.TXT file for more details.

% Brayan Alves, CIBM - MRS4Brain group, 2022
% Guillaume Briand, CIBM - MRS4Brain group, 2023
% 
% USAGE : MRSI Class public method 
% msg = obj.read_data()
% 
% DESCRIPTION :
% Read Bruker experiment folders and extract MRSI data and acquisition
% parameters
%
% INPUTS :
% obj       = MRSI Class object with properties and methods
%
% OUTPUT :
% msg       = Error message
function msg = read_data(obj)
%% Function to read Brucker data from the scanner on the MRSI matrix
msg = {''}; % Output initialization

fields = {'acqtype','path','fulldate','day','timeacq','sequence','acquisition_time','nucleus', ...
    'Bo','resfreq','matrix_sz','np_met','np_ref','spectralwidth','acq_delay',...
    'rep_time','nav','nrep','acq_freqshift','ppm_ref','ppm_workoffset','scale_ppm','grpdly',...
    'gain','acq_time_spec','voxs'};
c = cell(length(fields),1);
acq_params = cell2struct(c,fields);

try
    acq_params.path = fullfile(obj.data_folder,num2str(obj.metab_expnb));
    methodfilename = fullfile(obj.data_folder,num2str(obj.metab_expnb),'method');
    if ~isfile(methodfilename)
        msg = {'No existing method file in the metabolite folder'};
        return
    end
    methodfile = fileread(methodfilename);
    
    %date and time
    year = 1999;
    startind_date = [];
    while isempty(startind_date) % Find the experiment year
        startind_date = strfind(methodfile,['$$ ',num2str(year)]);
        year = year + 1;
    end
    timecode = methodfile(startind_date+3:startind_date+24);
    timecode = datetime(timecode,'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
    timeacq = hour(timecode);
    acq_params.fulldate = timecode;
    acq_params.day = datestr(timecode,'yyyymmdd'); %#ok<DATST>
    if timeacq > 11
        acq_params.timeacq = 'PM';
    else
        acq_params.timeacq = 'AM';
    end
    
    %literal names = acquisition time spectrum (sec), spectral width (Hz), sequence,
    % number of points metabolites (complex_nb*2), nucleus, MRSI matrix size,
    % resonance frequency B0, acquisition freq shift, ppm receiver offset,
    % ppm freq work offset, repetition time (sec), RG, voxel dimensions (mm),
    % nrep, nav, acq delay time (sec)
    if(obj.MRSI3D)
        BRUKERparamnames = {"PVM_SpecAcquisitionTime=","##$PVM_ScanTime=","##$PVM_SpecSWH=( 1 )", ...
            "##$Method=","##$PVM_SpecMatrix=( 1 )","##$PVM_Nucleus1Enum=", ...
            "##$PVM_Matrix=( 3 )" + newline, ...
            "##$PVM_FrqRef=( 8 )" + newline, ...
            "##$PVM_FrqWorkOffset=( 8 )" + newline, ...
            "##$PVM_FrqRefPpm=( 8 )" + newline, ...
            "##$PVM_FrqWorkOffsetPpm=( 8 )" + newline, ...
            "##$PVM_RepetitionTime=","##$PVM_RgValue=", ...
            "##$PVM_VoxArrSize=( 1, 3 )" + newline, ...
            "PVM_NRepetitions=","PVM_NAverages=","PVM_EchoTime="};
    else
        BRUKERparamnames = {"PVM_SpecAcquisitionTime=","##$PVM_ScanTime=","##$PVM_SpecSWH=( 1 )", ...
            "##$Method=","##$PVM_SpecMatrix=( 1 )","##$PVM_Nucleus1Enum=", ...
            "##$PVM_Matrix=( 2 )" + newline, ...
            "##$PVM_FrqRef=( 8 )" + newline, ...
            "##$PVM_FrqWorkOffset=( 8 )" + newline, ...
            "##$PVM_FrqRefPpm=( 8 )" + newline, ...
            "##$PVM_FrqWorkOffsetPpm=( 8 )" + newline, ...
            "##$PVM_RepetitionTime=","##$PVM_RgValue=", ...
            "##$PVM_VoxArrSize=( 1, 3 )" + newline, ...
            "PVM_NRepetitions=","PVM_NAverages=","PVM_EchoTime="};
    end
    
    GUIparamnames = {"acq_time_spec","acquisition_time","spectralwidth","sequence","np_met","nucleus",...
        "matrix_sz","resfreq","acq_freqshift","ppm_ref","ppm_workoffset",...
        "rep_time","gain","voxs","nrep","nav","acq_delay"};
    
    acq_params = fill_params(acq_params,methodfile,BRUKERparamnames,GUIparamnames);

    % Group Delay - acqus file
    if isfile(fullfile(acq_params.path,'acqus'))
        acqusfile = fileread(fullfile(acq_params.path,'acqus'));
        acq_params = fill_params(acq_params,acqusfile,{"##$GRPDLY= "},{"grpdly"});
    elseif isfile(fullfile(acq_params.path,'acqp'))
        acqpfile = fileread(fullfile(acq_params.path,'acqp'));
        acq_params = fill_params(acq_params,acqpfile,...
            {"##$ACQ_RxFilterInfo=( 2 )" + newline + "("},{"grpdly"});
    else
        msg = {'No existing acquistion parameters file in the metabolite folder'};
        return
    end

    % Number of points reference - method file ref
    methodfilename = fullfile(obj.data_folder,num2str(obj.ref_expnb),'method');
    if ~isfile(methodfilename)
        msg = {'No existing method file in the reference folder'};
        return
    end
    methodfile = fileread(methodfilename);
    acq_params = fill_params(acq_params,methodfile,{"##$PVM_SpecMatrix=( 1 )"},{"np_ref"});

    acq_params = params_adjustements(acq_params);
    assignin("base",'acq_params',acq_params)

    obj.acq_params = acq_params;

    % Paravision 3.3
    if ~isfile(fullfile(obj.data_folder,num2str(obj.metab_expnb),'fid'))
        metab_filename = fullfile(obj.data_folder,num2str(obj.metab_expnb),'pdata',num2str(obj.reco_expnb),'fid_proc.64');
        ref_filename = fullfile(obj.data_folder,num2str(obj.ref_expnb),'pdata',num2str(obj.reco_expnb),'fid_proc.64');
    % Paravision 1.1
    else
        if (obj.reco_expnb > 1)
            metab_filename = fullfile(obj.data_folder,num2str(obj.metab_expnb),'fid');
            ref_filename = fullfile(obj.data_folder,num2str(obj.ref_expnb),'fid');
        else
            metab_filename = fullfile(obj.data_folder,num2str(obj.metab_expnb),'fid_nofilter');
            ref_filename = fullfile(obj.data_folder,num2str(obj.ref_expnb),'fid_nofilter');
        end
    end
    
    fileid = fopen(metab_filename,'r','ieee-le'); % read binary format
    refid = fopen(ref_filename,'r','ieee-le');
    if fileid == -1
        msg = {'Cannot open the FID file, please check the correct experience number'};
        return
    elseif refid == -1
        msg = {'Cannot open the reference signal file, please check the correct experience number'};
        return
    end

    Time_t = (0 :(acq_params.np_met-1))'/acq_params.spectralwidth;
    Freqshift_1t = exp(2 * pi * 1j * Time_t * acq_params.acq_freqshift);
    
    buffer = fread(fileid,'double'); %note: CSI Bruker format is double 
    buffer_c = buffer(1:2:end) + 1j*buffer(2:2:end);
    
    buffer_ref = fread(refid,'double'); %note: CSI Bruker format is double 
    buffer_ref_c = buffer_ref(1:2:end) + 1j*buffer_ref(2:2:end);
    
    MatSize = acq_params.matrix_sz;
    fid_mat_c = reshape(buffer_c,[acq_params.np_met,MatSize(1)*MatSize(2),obj.Nslices]);
    ref_mat_c = reshape(buffer_ref_c,[acq_params.np_ref,MatSize(1)*MatSize(2),obj.Nslices]);
    
    grpdly = acq_params.grpdly;
    fid_mat_c_shift = [fid_mat_c(grpdly:end,:,:);0*fid_mat_c(1:grpdly-1,:,:)];
    ref_mat_c_shift = [ref_mat_c(grpdly:end,:,:);0*ref_mat_c(1:grpdly-1,:,:)];
    
    fid_mat_tkkn = conj(fft(fft(reshape(fid_mat_c_shift,[acq_params.np_met,MatSize(1),MatSize(2),obj.Nslices]),[],2),[],3)); % Restructuring the data + going to the k-space
    ref_mat_tkkn = conj(fft(fft(reshape(ref_mat_c_shift,[acq_params.np_ref,MatSize(1),MatSize(2),obj.Nslices]),[],2),[],3));
    
    obj.fid_mat_tkkn = fid_mat_tkkn .* Freqshift_1t; % We apply the frequency shift to center the water peak on the receiveroffset (ppm scale)
    obj.ref_mat_tkkn = ref_mat_tkkn .* Freqshift_1t; % Same for the reference signal
catch ME
    msg = {'Error while reading the spectroscopy data : ',ME.message};
end
end

%% Additional functions

%%  fill_study.m 
% Guillaume Briand, CIBM - MRS4Brain group, 2023
% 
% USAGE :
% acq_params = fill_params(acq_params,file,paramnames,GUIparamnames,substruct)
% 
% DESCRIPTION :
% Fill the structure with acquistion parameters from given file with specific
% acquistion parameter name and with or without substruct
%
% INPUTS :
% acq_params    = Structure with acquistion parameters
% file          = String of the given file with specific acquistion parameters 
% paramnames    = Literal names of parameters in file
% GUIparamnames = Names of parameters in study struct
% substruct     = Name of substrcture in study struct
%
% OUTPUT :
% acq_params    = Structure with acquistion parameters
function acq_params = fill_params(acq_params,file,paramnames,GUIparamnames,substruct)
if nargin < 5
    substruct = '';
end
for par = 1:length(paramnames)
    startind = strfind(file,paramnames{par}); % find the starting index in method file
    lparam = strlength(paramnames{par}); % length of the parameter in method file
    if ~isempty(startind)
        startind=startind(1);
        kparam = 1;
        while ~strcmp(file(startind + lparam + kparam),{'#','$',newline})
            kparam = kparam + 1;
        end
        param = string(file(startind + lparam : startind + lparam + kparam - 1));
        splitparam = strsplit(param);
        if ~isnan(str2double(param))
            param = str2double(param);
        elseif ~isnan(str2double(splitparam))
            param = str2double(splitparam);
        elseif ~isnan(str2double(splitparam{1}))
            param = str2double(splitparam{1});
        else
            param = char(param);
        end
        varname = matlab.lang.makeValidName(GUIparamnames{par}); % Ensure the fieldname is compatible
        if isempty(substruct)
            acq_params.(varname{1}) = param;
        else
            acq_params.(substruct).(varname{1}) = param;
        end
    end
end
end

%%  study_adjustements.m 
% Guillaume Briand, CIBM - MRS4Brain group, 2023
% 
% USAGE :
% [acq_params,msg] = params_adjustements(acq_params)
% 
% DESCRIPTION :
% Apply adjustements on specific acquistion parameters
%
% INPUTS :
% acq_params    = Structure with acquistion parameters
%
% OUTPUT :
% acq_params    = Structure with acquistion parameters updated
function acq_params = params_adjustements(acq_params)
% Adjustements
acq_params.acq_time_spec = acq_params.acq_time_spec*10^-3;          % in seconds
acq_params.acquisition_time = acq_params.acquisition_time*10^-3;    % in seconds
acq_params.acq_delay = acq_params.acq_delay*10^-3;                  % in seconds
acq_params.rep_time = acq_params.rep_time/1000;                     % in seconds 
acq_params.acqtype = 'MRSI/CSI';
acq_params.resfreq = acq_params.resfreq(1);
acq_params.acq_freqshift = abs(acq_params.acq_freqshift(1));
acq_params.ppm_ref = acq_params.ppm_ref(1);
acq_params.ppm_workoffset = acq_params.ppm_workoffset(1);
acq_params.grpdly = round(acq_params.grpdly) + 1;
acq_params.nucleus(regexp(acq_params.nucleus,'[<>]')) = [];
fmax = acq_params.spectralwidth/2;
acq_params.scale_f = fmax:-2*fmax/(acq_params.np_met-1):-fmax;
acq_params.scale_ppm = acq_params.scale_f./acq_params.resfreq + acq_params.ppm_ref;
if strcmp(acq_params.nucleus,'1H')
    acq_params.Bo = acq_params.resfreq / 42.576;
elseif strcmp(acq_params.nucleus,'2H')
    acq_params.Bo = acq_params.resfreq / 6.536;
elseif strcmp(acq_params.nucleus,'13C')
    acq_params.Bo = acq_params.resfreq / 10.7084;
elseif strcmp(acq_params.nucleus,'19F')
    acq_params.Bo = acq_params.resfreq / 40.078;
elseif strcmp(acq_params.nucleus,'31P')
    acq_params.Bo = acq_params.resfreq / 17.235;
end
end