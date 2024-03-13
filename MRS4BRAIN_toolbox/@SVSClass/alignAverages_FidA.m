%% alignAverages_FidA.m

% Copyright All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, MRS4Brain research group @ CIBM MRI EPFL AIT, 2024
% Copyright 2020 Jamie Near
% See the LICENSE.TXT file for more details.

% Jamie Near, McGill University 2014.
% Jessie Mosso, CIBM - MRS4Brain group, LIFMET, 2022 - modif
% 
% USAGE: DWS class method
% [out,fs,phs] = obj.op_alignAverages_fd(in,minppm,maxppm,tmax,med,ref);
% 
% DESCRIPTION:
% Perform time-domain spectral registration using a limited range of
% frequencies to correct frequency and phase drifts.  As described in Near
% et al.  Frequency and phase drift correction of magnetic resonance 
% spectroscopy data by spectral registration in the time domain. Magn Reson 
% Med 2015; 73(1):44-50.
% 
% INPUTS:
% obj       = DWS class object
% in        = Input data structure.
% minppm	= Minimum of frequency range (ppm).
% maxppm	= Maximum of frequnecy range (ppm).
% tmax      = Maximum time (s) in time domain to use for alignment.
% med       = Align averages to the median of the averages? ('y','n', 'a' or 
%             'r').  If you select 'n', all averages will be aligned to a 
%             single average.  The average chosen as the reference 
%             average will be the one with the lowest 'unlikeness' metric 
%             (see 'op_rmbadaverages.m').  If you select 'y', all
%             averages will be aligned to the median of the averages.  If
%             you select 'a', all averages will be aligned to the average
%             of the averages.  If you select 'r', all averages will be 
%             aligned to an externally provided reference spectrum.
% ref       = An externally provided reference spectrum that you would like
%             to align everything to (Required only if med = 'r').  
%
% OUTPUTS:
% out       = Output following alignment of averages.  
% fs        = Vector of frequency shifts (in Hz) used for alignment.
% phs       = Vector of phase shifts (in degrees) used for alignment.

function [out,fs,phs] = alignAverages_FidA(~,in,minppm,maxppm,tmax,med,ref)

if ~in.flags.addedrcvrs
    error('ERROR:  I think it only makes sense to do this after you have combined the channels using op_addrcvrs.  ABORTING!!');
end

if (strcmp(med,'r') || strcmp(med,'R'))
    if nargin < 7
        error('ERROR:  If using the ''r'' option for input variable ''med'', then a 6th input argument must be provided');
    end
else
    if nargin < 7
        ref=struct();
    end
end

parsFit = [0,0];

if in.dims.subSpecs == 0
    B = 1;
else
    B = in.sz(in.dims.subSpecs);
end

fs = zeros(in.sz(in.dims.averages),B);
phs = zeros(in.sz(in.dims.averages),B);
fids = zeros(in.sz(in.dims.t),1,B);
for m = 1:B
    if med == 'y' || med == 'Y'
        disp('Aligning all averages to the Average of the averages.');
        base = op_averaging(in);
        base = op_freqrange(base,minppm,maxppm);
        base = [real(base.fids(base.t>=0 & base.t<tmax,m));imag(base.fids(base.t>=0 & base.t<tmax,m))];
        ind_min = 0;
    elseif med == 'n' || med == 'N'
        %First find the average that is most similar to the total average:
        inavg=op_median(in);
        metric = zeros(in.sz(in.dims.averages),B);
        for k = 1:in.sz(in.dims.averages)
            for l = 1:B
                metric(k,l) = sum((real(in.fids(in.t >= 0 & in.t <= tmax,k,l)) -  ...
                    (real(inavg.fids(inavg.t >= 0 & inavg.t <= tmax,l)))).^2);
            end
        end
        [~,ind_min] = min(metric(:,m));
        
        %Now set the base function using the index of the most similar
        %average:
        disp(['Aligning all averages to average number ' num2str(ind_min) '.']);
        base = op_freqrange(in,minppm,maxppm);
        base = [real(base.fids(base.t>=0 & base.t<tmax,ind_min,m));imag(base.fids(base.t>=0 & base.t<tmax,ind_min,m))];
        fids(:,ind_min,m) = in.fids(:,ind_min,m);
    elseif med == 'r' || med == 'R'
        disp('Aligning all averages to an externally provided reference spectrum.');
        base = ref;
        base = op_freqrange(base,minppm,maxppm);
        base = [real(base.fids(base.t>=0 & base.t<tmax,m));imag(base.fids(base.t>=0 & base.t<tmax,m))];
        ind_min = 0;
    end
    for n = 1:in.sz(in.dims.averages)
        if n ~= ind_min
            parsGuess = parsFit;
            %parsGuess(1) = parsGuess(1);
            %disp(['fitting subspec number ' num2str(m) ' and average number ' num2str(n)]);
            datarange = op_freqrange(in,minppm,maxppm);
            start = datarange.fids(datarange.t >= 0 & datarange.t < tmax,n,m);
            parsFit = nlinfit(start,base,@op_freqPhaseShiftComplexRangeNest,parsGuess);
            fids(:,n,m) = op_freqPhaseShiftNest(parsFit,in.fids(:,n,m));
            fs(n,m) = parsFit(1);
            phs(n,m) = parsFit(2);
            %plot(in.ppm,fftshift(ifft(fids(:,1,m))),in.ppm,fftshift(ifft(fids(:,n,m))));
        end
    end
end

fids = conj(fids);
%re-calculate Specs using fft
specs = fftshift(ifft(fids,[],in.dims.t),in.dims.t);

%FILLING IN DATA STRUCTURE
out = in;
out.fids = fids;
out.specs = specs;

%FILLING IN THE FLAGS
out.flags = in.flags;
out.flags.writtentostruct = 1;
out.flags.freqcorrected = 1;

    function y = op_freqPhaseShiftComplexRangeNest(pars,input)
        f = pars(1);     %Frequency Shift [Hz]
        p = pars(2);     %Phase Shift [deg]

        dwelltime = datarange.dwelltime;
        t = 0:dwelltime:(length(input)-1)*dwelltime;
        fid = input(:);
        
        shifted = addphase(fid.*exp(1i*t'*f*2*pi),p);
        
        y = [real(shifted);imag(shifted)];
    end

    function y = op_freqPhaseShiftNest(pars,input)
        f = pars(1);     %Frequency Shift [Hz]
        p = pars(2);     %Phase Shift [deg]
  
        dwelltime = in.dwelltime;
        t = 0:dwelltime:(length(input)-1)*dwelltime;
        fid = input(:);
        
        y = addphase(fid.*exp(1i*t'*f*2*pi),p);
    end
end

%% Additional functions

%% op_averaging.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% out=op_averaging(in);
% 
% DESCRIPTION:
% Combine the averages in a scan by adding the averages together and then 
% dividing by the number of averages.
% 
% INPUTS:
% in	= input data in matlab structure format.
%
% OUTPUTS:
% out   = Output following averaging.  

function out = op_averaging(in)
if in.flags.averaged || in.averages<2
    %DO NOTHING
    disp('WARNING: No averages found. Returning input without modification!');
    return;
end
if in.dims.averages == 0
    %DO NOTHING
    disp('WARNING: No averages found. Returning input without modification!');
    out = in;
    return;
else
    %add the spectrum along the averages dimension;
    fids = sum(in.fids,in.dims.averages);
    fids = squeeze(fids);
    fids = fids/in.sz(in.dims.averages); %divide by number of averages;
    
    %re-calculate Specs using fft
    specs = fftshift(ifft(fids,[],in.dims.t),in.dims.t);
    
    %change the dims variables.
    if in.dims.t > in.dims.averages
        dims.t = in.dims.t-1;
    else
        dims.t = in.dims.t;
    end
    if in.dims.coils > in.dims.averages
        dims.coils = in.dims.coils-1;
    else
        dims.coils = in.dims.coils;
    end
    dims.averages = 0;
    if in.dims.subSpecs > in.dims.averages
        dims.subSpecs = in.dims.subSpecs-1;
    else
        dims.subSpecs = in.dims.subSpecs;
    end
    if in.dims.extras > in.dims.averages
        dims.extras = in.dims.extras-1;
    else
        dims.extras = in.dims.extras;
    end
       
    %re-calculate the sz variable
    sz = size(fids);
    
    %FILLING IN DATA STRUCTURE
    out = in;
    out.fids = fids;
    out.specs = specs;
    out.sz = sz;
    out.dims = dims;
    out.averages = 1;
    
    %FILLING IN THE FLAGS
    out.flags = in.flags;
    out.flags.writtentostruct = 1;
    out.flags.averaged = 1;
end
end

%% op_freqrange.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% out=op_freqrange(in,ppmmin,ppmmax);
% 
% DESCRIPTION:
% Output only a specified frequency range of the input spectrum.
% 
% INPUTS:
% in         = input data in matlab structure format.
% ppmmin     = minimum extent of frequency range in ppm.
% ppmmax     = maximum extent of frequency range in ppm.
%
% OUTPUTS:
% out        = Output following frequency range selection.

function out = op_freqrange(in,ppmmin,ppmmax)
%Calculate Specs using fft
fullspecs = fftshift(ifft(in.fids,[],in.dims.t),in.dims.t);
%now take only the specified range of the spectrum
specs = fullspecs(in.ppm > ppmmin & in.ppm < ppmmax,:,:);
%convert back to time domain
%if the length of Fids is odd, then you have to do a circshift of one to
%make sure that you don't introduce a small frequency shift into the fids
%vector.
if mod(size(specs,in.dims.t),2) == 0
    %disp('Length of vector is even.  Doing normal conversion');
    fids = fft(fftshift(specs,in.dims.t),[],in.dims.t);
else
    %disp('Length of vector is odd.  Doing circshift by 1');
    fids = fft(circshift(fftshift(specs,in.dims.t),1),[],in.dims.t);
end

%calculate the size;
sz = size(fids);

%calculate the ppm scale
ppm = in.ppm(in.ppm > ppmmin & in.ppm < ppmmax);

%calculate the new spectral width and dwelltime:
dppm = abs(ppm(2)-ppm(1));
ppmrange = abs((ppm(end)-ppm(1))) + dppm;
spectralwidth = ppmrange*in.txfrq;
dwelltime = 1/spectralwidth;

%calculate the time scale
t = 0:dwelltime:(sz(1)-1)*dwelltime;

%FILLING IN DATA STRUCTURE
out = in;
out.fids = fids;
out.specs = specs;
out.sz = sz;
out.ppm = ppm;  
out.t = t; 
out.spectralwidth = spectralwidth;
out.dwelltime = dwelltime;

%FILLING IN THE FLAGS
out.flags = in.flags;
out.flags.writtentostruct = 1;
out.flags.freqranged = 1;
end

%% op_median.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% out=op_median(in);
% 
% DESCRIPTION:
% Combine the averages in a scan by calculating the median of all averages.
% 
% INPUTS:
% in	= input data in matlab structure format.
%
% OUTPUTS:
% out   = Output dataset following median calculation.

function out = op_median(in)
if in.flags.averaged || in.dims.averages == 0 || in.averages < 2
    error('ERROR:  Averaging has already been performed!  Aborting!');
end
%add the spectrum along the averages dimension;
fids = median(real(in.fids),in.dims.averages)+(1i*median(imag(in.fids),in.dims.averages));
fids = squeeze(fids);
%re-calculate Specs using fft
specs = fftshift(ifft(fids,[],in.dims.t),in.dims.t);
%change the dims variables.  
if in.dims.t > in.dims.averages
    dims.t = in.dims.t-1;
else
    dims.t = in.dims.t;
end
if in.dims.coils > in.dims.averages
    dims.coils = in.dims.coils-1;
else
    dims.coils = in.dims.coils;
end
dims.averages = 0;
if in.dims.subSpecs > in.dims.averages
    dims.subSpecs = in.dims.subSpecs-1;
else
    dims.subSpecs = in.dims.subSpecs;
end
if in.dims.extras > in.dims.averages
    dims.extras = in.dims.extras-1;
else
    dims.extras = in.dims.extras;
end
%re-calculate the sz variable
sz = size(fids);

%FILLING IN DATA STRUCTURE
out = in;
out.fids = fids;
out.specs = specs;
out.sz = sz;
out.dims = dims;
out.averages = 1;

%FILLING IN THE FLAGS
out.flags = in.flags;
out.flags.writtentostruct = 1;
out.flags.averaged = 1;
end

%% addphase.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% PhasedSpecs=addphase(specs,AddedPhase);
% 
% DESCRIPTION:  
% Add equal amount of complex phase to each point of a vector.  This
% function operates on a vector (fid or spectrum), not on a FID-A data
% structure.  For a phase shifting function that operates on a FID-A data
% structure, see 'op_addphase.m'.
% 
% INPUTS:
% specs          = Input vector.
% AddedPhase     = Amount of phase (degrees) to add.
%
% OUTPUTS: 
% PhasedSpecs    = Output vector (0th order phased version of the input). 

function PhasedSpecs = addphase(specs,AddedPhase)
PhasedSpecs = specs.*(ones(size(specs))*exp(1i*AddedPhase*pi/180));
end
