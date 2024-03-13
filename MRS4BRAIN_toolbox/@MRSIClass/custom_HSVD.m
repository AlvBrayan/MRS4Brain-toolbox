%%  custom_HSVD.m 

% Copyright All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, MRS4Brain research group @ CIBM MRI EPFL AIT, 2024
% See the LICENSE.TXT file for more details.

% Brayan Alves, CIBM - MRS4Brain group, 2022
% Guillaume Briand, CIBM - MRS4Brain group, 2023
% 
% USAGE : MRSI Class public method 
% msg = obj.custom_HSVD(prog_dbox)
% 
% DESCRIPTION :
% apply Hankel Lanczos Singular value decomposition (HLSVD) on MRSI data
%
% INPUTS :
% obj       = MRSI Class object with properties and methods
% prog_dbox = MRS4Brain Toolbox progress dialog box
%
% OUTPUT :
% msg       = Error message
function msg = custom_HSVD(obj,prog_dbox)
msg = {''}; % Output initialization

obj.HSVD_fid_tkkn = obj.fid_mat_tkkn; % metabolite HSVD fid initialization
samplerate = obj.acq_params.spectralwidth;
FiltParam.Comp = 16;
FiltParam.Water_minFreq = -obj.acq_params.resfreq/2;
FiltParam.Water_maxFreq = obj.acq_params.resfreq/2;

for ii = 1:obj.Nslices
    MRSIData_tkk = squeeze(obj.fid_mat_tkkn(:,:,:,ii)); % Fourier space of 1 slice
    Brain_mask = squeeze(obj.Brain_mask(ii,:,:));
    MRSIData_trr = ifft(ifft(MRSIData_tkk,[],2),[],3); % Real space of 1 slice
    try
        for b = 1:size(MRSIData_tkk,3)
            parfor a = 1:size(MRSIData_tkk,2)
                % disp(['Applying HSVD on voxel (', num2str(a), ',', num2str(b), ')'])
                if Brain_mask(a,b) == 1 % Check if in the brain
                    [msg_err, MRSIData_trr(:,a,b)] = Fast_HSVD_Filter(...
                        MRSIData_trr(:,a,b),samplerate,FiltParam);
                    if ~isempty(msg_err)
                        disp(msg_err)
                        msg = {msg_err}; %#ok<PFTUSW>
                    end
                end
            end
            if any(~cellfun(@isempty,msg))
                return
            end
            if prog_dbox.CancelRequested
                msg = {'Unfinished business due to cancellation'};
                return
            end
        end
        filtMRSI = fft(fft(MRSIData_trr,[],2),[],3);
    catch ME
        msg = {'HSVD computation error',ME.message};
        return
    end
%     assignin("base",'filtMRSI',filtMRSI)
    if exist('filtMRSI','var')
        obj.HSVD_fid_tkkn(:,:,:,ii) = filtMRSI;
    else
        obj.HSVD_fid_tkkn(:,:,:,ii) = MRSIData_tkk;
    end
end
end

%% HSVD function

function [msg, filtMRSI] = Fast_HSVD_Filter(mrsiData,tempSamplerate,FiltParam)
% Output variables
msg = '';
filtMRSI = zeros(size(mrsiData));

minFreq = FiltParam.Water_minFreq; %in Hz
maxFreq = FiltParam.Water_maxFreq; %in Hz
try
    [freqs,~,basis,amps] = HSVD(mrsiData, tempSamplerate, FiltParam.Comp);
    indx = find((freqs >= minFreq) & (freqs <= maxFreq));
    filtMRSI = mrsiData - sum(basis(:, indx) * diag(amps(indx), 0), 2);
catch ME
    msg = ME.mesage;
    return
end
end

function [frequencies, dampings, basis, ahat] = HSVD(y, fs, K)
% Decompose the signal y using the method of Barkhuijsen, et al. (1987)
%
% Obligatory arguments:
% 'y' is the FID (a linear combination of complex exponentials).
% 'fs' is the sampling frequency (bandwidth).
% 'K' is the desired model order (expected number of complex
% exponentials comprising the signal)
%
% Outputs:
% 'frequencies' - frequencies of components in signal
% 'dampings' - damping of components in signal
% 'basis' - basis vectors, one for each component
% 'ahat' - amplitudes of each basis in signal
%
% Author: Greg Reynolds (remove.this.gmr001@bham.ac.uk)
% Date: August 2006 

N = length(y);
L = floor(0.5*N);
% M = N+1-L;

% H is the LxM Hankel LP matrix
H       = hankel(y(1:L), y(L:N));

% compute the singular value decomposition
[U,~,~] = svd(H);

% construct H of rank K
Uk = U(:,1:K);

% find Ukt and Ukb
Ukt = Uk(2:end,:);
Ukb = Uk(1:end-1,:);
Zp  = pinv(Ukb)*Ukt;

% find the poles
[Q,~]       = eig(Zp);
Z           = Q\Zp*Q; % inv(Q)*Zp*Q
q           = log(diag(Z));
dt          = 1/fs;
dampings    = real(q)/dt;
frequencies = imag(q)/(2*pi)/dt;

%dampings(dampings>10)
dampings(dampings>10) = 10;

% construct the basis
t     = (0:dt:(length(y)-1)*dt);
basis = exp(t.'*(dampings.' + 2*pi*1i*frequencies.'));

% compute the amplitude estimates
ahat = pinv(basis(1:length(y),:))*y;
end