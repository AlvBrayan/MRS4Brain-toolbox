%%  plot_signal.m 

% Copyright All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, MRS4Brain research group @ CIBM MRI EPFL AIT, 2024
% See the LICENSE.TXT file for more details.

% Guillaume Briand, CIBM - MRS4Brain group, 2023
% 
% USAGE : SVS Class public method 
% [scale_ppm,data_spec] = obj.plot_signal(exp,multi,type)
% 
% DESCRIPTION :
% Return ppm scale and data spectrum for current experiment and average
%
% INPUTS :
% obj       = SVS Class object with properties and methods
% exp       = Experiment number
% multi     = For multi data fids, the one wanted
% type      = type of study (raw,processed,sum+processed)
%
% OUTPUT :
% scale_ppm = ppm scale to plot data
% data_spec = Fourier transformed data fids
function [scale_ppm,data_spec] = plot_signal(obj,exp,multi,type)

if nargin < 4
    type = 'sum_processed';
end

if strcmp(type,'raw')
    study = obj.SVS_struct(exp).raw_study;
elseif strcmp(type,'processed')
    study = obj.SVS_struct(exp).processed_study;
elseif strcmp(type,'sum_processed')
    study = obj.SVS_struct(exp).sum_processed_study;
else
    study = obj.SVS_struct(exp).sum_processed_study;
end
fmax = study.params.sw/2;
% Data ppm scale
f_data = fmax:-2*fmax/(floor(study.np/2) - 1):-fmax;
scale_ppm = f_data/study.params.sfrq + study.ppm_ref - study.ppm_workoffset;

data_conj = squeeze(study.data.real(multi,:,:)) + 1j*squeeze(study.data.imag(multi,:,:));
data_spec = fftshift(fft(data_conj,[],1),1);
end