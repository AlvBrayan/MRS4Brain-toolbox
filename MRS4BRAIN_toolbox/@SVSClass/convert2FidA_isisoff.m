%%  convert2FidA_isisoff.m 

% Copyright All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, MRS4Brain research group @ CIBM MRI EPFL AIT, 2024
% Copyright 2020 Jamie Near
% See the LICENSE.TXT file for more details.

% Jessie Mosso, CIBM - MRS4Brain group, LIFMET, 2021
% Guillaume Briand, CIBM - MRS4Brain group, 2023
% 
% USAGE : SVS Class public method 
% out = obj.convert2FidA_isisoff(nb)
% 
% DESCRIPTION :
% Convert a study structure to FidA format without ISIS
%
% INPUTS :
% obj       = SVS Class object with properties and methods
% idx_fid   = Index of the studied fid in obj.SVS_struct
%
% OUTPUT :
% out       = FidA study format
function out = convert2FidA_isisoff(obj,nb)
study = obj.SVS_struct(nb).raw_study;
% Load the study struct and convert it to FidA format
out.fids = squeeze(study.data.real) + 1i*squeeze(study.data.imag);
grpdly = round(study.params.grpdly) + 1; 
np = study.np/2;
out.fids = [out.fids(:,grpdly:end),zeros(size(out.fids,1),grpdly-1)];
out.specs = fftshift(fft(out.fids,[],2),2);
out.fids = out.fids.'; 
out.specs = out.specs.'; 
out.spectralwidth = round(study.spectralwidth,0);
out.dwelltime = 1/out.spectralwidth;
out.t = 0:out.dwelltime:(np-1)*out.dwelltime;
fmax = out.spectralwidth/2;
f = fmax:-2*fmax/(np-1):-fmax;
ppmscale = f/study.resfreq + study.ppm_ref - study.ppm_workoffset;
out.sz = [size(out.fids,1),size(out.fids,2),1,1];
out.date = study.day;
out.averages = size(out.fids,2);
out.rawAverages = size(out.fids,2);
out.subspecs = 1; 
out.rawSubspecs = 1; 
out.txfrq = study.resfreq;
% Magnetic field Bo
if strcmp(study.nucleus,'1H')
    out.Bo = study.resfreq / 42.576;
elseif strcmp(study.nucleus,'2H')
    out.Bo = study.resfreq / 6.536;
elseif strcmp(study.nucleus,'13C')
    out.Bo = study.resfreq / 10.7084;
elseif strcmp(study.nucleus,'19F')
    out.Bo = study.resfreq / 40.078;
elseif strcmp(study.nucleus,'31P')
    out.Bo = study.resfreq / 17.235;
end
out.n = np;
out.dims.t = 1;
out.dims.averages = 2;
out.dims.coils = 0; 
out.dims.subSpecs = 0;
out.dims.extras = 0;
out.ppm = ppmscale;
out.flags.writtentostruct = 1;
out.flags.gotparams = 1;
out.flags.filtered = 0;
out.flags.zeropadded = 0;
out.flags.freqcorrected = 0;
out.flags.phasecorrected = 0;
out.flags.averaged = 0;
out.flags.addedrcvrs = 1; 
out.flags.Subtracted = 1;
out.flags.Writtentotext = 0;
out.flags.Downsampled = 0;
out.flags.avgNormalized = 0;
out.flags.isISIS = 0;
end
 