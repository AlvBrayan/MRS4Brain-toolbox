%%  create_Bruker_rawstudy.m 

% Copyright All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, MRS4Brain research group @ CIBM MRI EPFL AIT, 2024
% See the LICENSE.TXT file for more details.

% Jessie Mosso, CIBM - MRS4Brain group, LIFMET, 2021
% Guillaume Briand, CIBM - MRS4Brain group, 2023
% 
% USAGE : DWS Class public method 
% msg = obj.create_Bruker_rawstudy(prog_dbox)
% 
% DESCRIPTION :
% Create from a Bruker experiment folder and experiment number a study
% structure
%
% INPUTS :
% obj       = DWS Class object with properties and methods
% prog_dbox = MRS4Brain Toolbox progress dialog box
% OUTPUT :
% msg       = Error message
function msg = create_Bruker_rawstudy(obj,prog_dbox)
%% From a Bruker exp. folder, generates MATLAB structures with necessary 
msg = {''};
refscanbool = false;
if(~exist(fullfile(obj.result_dir,obj.foldername),"dir"))
    mkdir(fullfile(obj.result_dir,obj.foldername));
end
for i = 1:length(obj.DWS_struct)
    prog_dbox.Message = ['Create raw study for : ',obj.DWS_struct(i).exp_name];
    prog_dbox.Value = i/(length(obj.DWS_struct)+1);
    study = struct;
    folderexp = obj.DWS_struct(i).folder;
    expnb = obj.DWS_struct(i).expnb;
    %% Parameters - method file
    study.path = fullfile(folderexp,num2str(expnb));
    methodfilename = fullfile(folderexp,num2str(expnb),'method');
    if ~isfile(methodfilename)
        msg = {'No existing method file in the directory'};
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
    study.fulldate = timecode;
    study.day = datestr(timecode,'yyyymmdd'); %#ok<DATST>
    if timeacq > 11
        study.timeacq = 'PM';
    else
        study.timeacq = 'AM';
    end
    
    
    %literal names = acquisition time (sec), spectral width (Hz), sequence, np
    % (complex nb *2), nucleus, resonance frequency B0, ppm receiver offset,
    % ppm freq work offset, tr(sec), RG, voxel dimensions (mm), acq type (MRS),
    % nrep, nav, echo time
    BRUKERparamnames = {"PVM_SpecAcquisitionTime=","##$PVM_SpecSWH=( 1 )", ...
        "##$Method=","##$PVM_SpecMatrix=( 1 )","##$PVM_Nucleus1Enum=", ...
        "##$PVM_FrqRef=( 8 )" + newline, ...
        "##$PVM_FrqRefPpm=( 8 )" + newline, ...
        "##$PVM_FrqWorkOffsetPpm=( 8 )" + newline, ...
        "##$PVM_RepetitionTime=","##$PVM_RgValue=", ...
        "##$PVM_VoxArrSize=( 1, 3 )" + newline,"##$PVM_EncSpectroscopy=", ...
        "PVM_NRepetitions=","PVM_NAverages=","PVM_EchoTime="};
    
    %literal names = delta, b-value, diffusion time (sec), mixing time (sec),
    % diffusion gradient in (x,y,z), diffusion direction in (x,y,z)
    BRUKER_DWSparamnames = {"##$Delta=","##$Bvalue=","##$DiffusionTime=","##$MixingTime=", ...
        "##$GdiffX=","##$GdiffY=","##$GdiffZ=","##$DWdirX=","##$DWdirY=","##$DWdirZ="};
    
    GUIparamnames = {"acq_time","spectralwidth","sequence","np","nucleus",...
        "resfreq","ppm_ref","ppm_workoffset","tr","gain","voxs","acqtype","nrep","nav","te"};
    
    GUI_DWSparamnames = {"delta","bvalue","diff_time","mixing_time","gdiff_x","gdiff_y", ...
        "gdiff_z","dwdir_x","dwdir_y","dwdir_z"};
    
    % Fill the study struct
    study = fill_study(study,methodfile,BRUKERparamnames,GUIparamnames);
    study = fill_study(study,methodfile,BRUKER_DWSparamnames,GUI_DWSparamnames,'dwsparams');
    
    % Make adjustements to the study struct and check that it is MRS
    [study,msg] = study_adjustements(study);
    if ~any(cellfun(@isempty,msg))
        msg = {msg{1},['Experiment name : ' obj.DWS_struct(i).exp_name]};
        return
    end
    
    % Group Delay - acqus file
    if isfile(fullfile(folderexp,num2str(expnb),'acqus'))
        acqusfile = fileread(fullfile(folderexp,num2str(expnb),'acqus'));
        study = fill_study(study,acqusfile,{"##$GRPDLY= "},{"grpdly"},'params');
    elseif isfile(fullfile(folderexp,num2str(expnb),'acqp'))
        acqpfile = fileread(fullfile(folderexp,num2str(expnb),'acqp'));
        study = fill_study(study,acqpfile,{"##$ACQ_RxFilterInfo=( 2 )" + newline + "("},{"grpdly"},'params');
    end
    
    % ADCoverflow? - acqp file
    acqpfile = fileread(fullfile(folderexp,num2str(expnb),'acqp'));
    study = fill_study(study,acqpfile,{"##$ACQ_adc_overflow=( 2 )" + newline},{"adcoverflow"},'params');
    adcoverflow = strsplit(study.params.adcoverflow);
    if or(strcmp(adcoverflow{1},'Yes'),strcmp(adcoverflow{2},'Yes'))
        study.params.adcoverflow = 'Yes';
        % f=msgbox(['ADC overflow during acquisition E' num2str(expnb)]);
    else
        study.params.adcoverflow = 'No';
    end 
    
    %% Read data
    
    %job0
    if isfile(fullfile(folderexp,num2str(expnb),'ser'))
        fileid = fopen(fullfile(folderexp,num2str(expnb),'ser'),'r','ieee-le'); %read binary format
        if fileid == -1
            msg = {'Cannot open the FID file',['Experiment name : ' obj.DWS_struct(i).exp_name]};
            return
        end
    elseif isfile(fullfile(folderexp,num2str(expnb),'fid'))
        if study.nrep == 1
            dimtype="1D";
        else
            dimtype="2D";
        end 
        disp (dimtype + " fid, " + "nav=" + num2str(study.nav) + " nrep=" + num2str(study.nrep)+ " saved")
        fileid = fopen(fullfile(folderexp,num2str(expnb),'fid'),'r','ieee-le'); %read binary format
        if fileid == -1
            msg = {'Cannot open the FID file',['Experiment name : ' obj.DWS_struct(i).exp_name]};
            return
        end
    elseif isfile(fullfile(folderexp,num2str(expnb),'fid','pdata','1','fid_proc.64'))
        fileid = fopen(fullfile(folderexp,num2str(expnb),'fid','pdata','1','fid_proc.64'),'r','ieee-le'); %read binary format
        if fileid == -1
            msg = {'Cannot open the FID file',['Experiment name : ' obj.DWS_struct(i).exp_name]};
            return
        end
    end
    
    buffer = fread(fileid,'int32'); 
    nbptsfid = length(buffer)/(2*study.nrep);
    
    buffer_ser = zeros(study.nrep,nbptsfid*2);
    for rep = 1:study.nrep
        buffer_ser(rep,:) = buffer((rep-1)*(nbptsfid*2)+1:rep*(nbptsfid*2))';
    end
    ser_c = buffer_ser(:,1:2:end) + 1j*buffer_ser(:,2:2:end);
    fclose(fileid);
    
    %offset 
    if study.ppm_workoffset ~= 0
        offset_hz = study.ppm_workoffset*study.resfreq;
        dw = 1/study.spectralwidth;
        timecode = 0:dw:(study.np/2-1)*dw;
        tmat = repmat(timecode,study.nrep,1);
        ser_c_shift = ser_c.*exp(1i.*(2*pi*offset_hz).*tmat);
        study.data.real(:,1,:) = real(ser_c_shift);
        study.data.imag(:,1,:) = -imag(ser_c_shift); %flips the spectrum
    else
        study.data.real(:,1,:) = real(ser_c);
        study.data.imag(:,1,:) = -imag(ser_c); %flips the spectrum
    end 
    
    %filename and liststring
    if study.nrep>1
        filename=['Bruker_' study.day '_' study.timeacq '_' num2str(expnb) '_ser.mat'];
        study.params.nt = study.nrep*study.nav;
    
    else
        filename=['Bruker_' study.day '_' study.timeacq '_' num2str(expnb) '_fid.mat'];
        study.params.nt = study.nav;
    end 
    
    %normalize: xNA/RG/voxelsize
    voxvol = study.params.vox1 .* study.params.vox2 .* study.params.vox3;
    study.data.real = study.data.real .* study.nav ./ study.params.gain ./ voxvol;
    study.data.imag = study.data.imag .* study.nav ./ study.params.gain ./ voxvol;
    
    study.filename = filename;
    study.liststring = char(fullfile(obj.result_dir,obj.foldername,'raw',study.filename));
    
    %multiplicity
    study.multiplicity = study.nrep;
    
    %process
    study.process.lsfid = round(study.params.grpdly) + 1;
    study.process.apodizefct = 'exponential'; %default
    study.process.apodparam1 = zeros(1,study.multiplicity); %default
    study.process.apodparam2 = zeros(1,study.multiplicity); %default
    study.process.transfsize = 0; %default
    study.process.appltoarray1 = 0; %default
    study.process.appltoarray2 = 0; %default
    study.process.phasecorr0 = zeros(1,study.multiplicity); %default
    study.process.phasecorr1 = zeros(1,study.multiplicity); %default
    study.process.DCoffset = 0; %default
    study = rmfield(study,'nav');
    study = rmfield(study,'nrep');
    
    raw_study = study;
    
    if(~exist(fullfile(obj.result_dir,obj.foldername,'raw'),"dir"))
        mkdir(fullfile(obj.result_dir,obj.foldername,'raw'));
    end
    obj.DWS_struct(i).raw_study = raw_study;
    
    save(study.liststring,'raw_study')
end
%% Refscan

if and(isfile(fullfile(folderexp,num2str(expnb),'fid.refscan')),refscanbool)

    %literal names = navrefscan,
    BRUKERparamnames_refscan={"PVM_RefScanNA=","##$PVM_RefScanRG="};
    localparamnames_refscan={"NArefscan","RGrefscan"};
    refscan = struct;
    refscan = fill_study(refscan,methodfile,BRUKERparamnames_refscan,localparamnames_refscan);

    %store data
    fileid=fopen(fullfile(folderexp,num2str(expnb),'fid.refscan'),'r','ieee-le'); %read binary format
    if fileid == -1
        msg = {'Cannot open the FID file',['Experiment name : ' obj.DWS_struct(i).exp_name]};
        return
    end
    buffer = fread(fileid,'int32'); 
    % nbptsfid = length(buffer)/2;
    buffer_c_ref = buffer(1:2:end) + 1j*buffer(2:2:end); 

    fclose(fileid);

    study.data.real = zeros(1,1,study.np/2);
    study.data.imag = zeros(1,1,study.np/2);
    study.data.real(1,1,:) = real(buffer_c_ref);
    study.data.imag(1,1,:) = -imag(buffer_c_ref); %flips the spectrum

    %filename and liststring
    filename=['Bruker_' study.day '_' study.timeacq '_' num2str(expnb) '_fidrefscan.mat'];
    study.filename = filename;
    study.liststring = char(fullfile(folderexp,study.filename));


    %normalize: xNA/RG/voxelsize 
    study.data.real = study.data.real .* refscan.NArefscan ./ refscan.RGrefscan ./ voxvol;
    study.data.imag = study.data.imag .* refscan.NArefscan ./ refscan.RGrefscan ./ voxvol;

    %multiplicity
    study.multiplicity = 1; %1D

    %process
    study.process.lsfid=round(study.params.grpdly) + 1;
    study.process.apodizefct = 'exponential'; %default
    study.process.apodparam1 = 0;  %default
    study.process.apodparam2 = 0; %default
    study.process.transfsize = 0; %default
    study.process.appltoarray1 = 0; %default
    study.process.appltoarray2 = 0; %default
    study.process.phasecorr0 = 0; %default
    study.process.phasecorr1 = 0; %default
    study.process.DCoffset = 0; %default
    study.process.B0 = zeros(1,study.np/2);

    %nt
    study.params.nt = refscan.NArefscan;
    
    raw_study = study;
    
    save(raw_study.liststring,'raw_study');

end
end
    
%% Additional functions

%%  fill_study.m 
% Guillaume Briand, CIBM - MRS4Brain group, 2023
% 
% USAGE :
% study = fill_study(study,file,paramnames,GUIparamnames,substruct)
% 
% DESCRIPTION :
% Fill the study structure with parameters from given file with specific
% acquistion parameter name and with or without substruct
%
% INPUTS :
% study         = Study structure with acquistion parameters and data
% file          = String of the given file with specific acquistion parameters 
% paramnames    = Literal names of parameters in file
% GUIparamnames = Names of parameters in study struct
% substruct     = Name of substrcture in study struct
%
% OUTPUT :
% study         = Study structure with acquistion parameters and data
function study = fill_study(study,file,paramnames,GUIparamnames,substruct)
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
            study.(varname{1}) = param;
        else
            study.(substruct).(varname{1}) = param;
        end
    end
end
end

%%  study_adjustements.m 
% Guillaume Briand, CIBM - MRS4Brain group, 2023
% 
% USAGE :
% [study,msg] = study_adjustements(study)
% 
% DESCRIPTION :
% Apply adjustements on specific study structure parameters
%
% INPUTS :
% study         = Study structure with acquistion parameters and data
%
% OUTPUT :
% study         = Study structure with acquistion parameters and data
% msg           = Error message for correct acquistion type
function [study,msg] = study_adjustements(study)
msg = {''};
% Adjustements
study.acq_time = study.acq_time*10^-3;  % in seconds
study.np = study.np*2;                  % real/imag
study.tr = study.tr/1000;               % in seconds 
if strcmp(study.acqtype,'Yes')
    study.acqtype = 'MRS';
else 
    msg = {"Can't open this type of data"};
end
study.resfreq = study.resfreq(1);
study.ppm_ref = study.ppm_ref(1);
study.ppm_workoffset = study.ppm_workoffset(1);
study.format = 'Matlab';
study.nucleus(regexp(study.nucleus,'[<>]')) = [];
study.params.sfrq = study.resfreq;
if strcmp(study.nucleus,'1H')
    study.params.Bo = study.resfreq / 42.576;
elseif strcmp(study.nucleus,'2H')
    study.params.Bo = study.resfreq / 6.536;
elseif strcmp(study.nucleus,'13C')
    study.params.Bo = study.resfreq / 10.7084;
elseif strcmp(study.nucleus,'19F')
    study.params.Bo = study.resfreq / 40.078;
elseif strcmp(study.nucleus,'31P')
    study.params.Bo = study.resfreq / 17.235;
end
study.params.arraydim = 1; %MRS
study.params.np = study.np;
study.params.sw = study.spectralwidth;
study.params.tr = study.tr;
study=rmfield(study,'tr');
study.params.te=study.te; 
study=rmfield(study,'te');
study.params.gain=study.gain;
study=rmfield(study,'gain');
study.params.vox1 = study.voxs(1);
study.params.vox2 = study.voxs(2);
study.params.vox3 = study.voxs(3);
study=rmfield(study,'voxs');
end
