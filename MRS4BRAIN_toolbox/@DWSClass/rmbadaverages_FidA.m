%% rmbadaverages_FidA.m

% Copyright All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, MRS4Brain research group @ CIBM MRI EPFL AIT, 2024
% Copyright 2020 Jamie Near
% See the LICENSE.TXT file for more details.

% Jamie Near, McGill University 2014.
% Jessie Mosso, CIBM - MRS4Brain group, LIFMET, 10/12/2022 - modif
%
% USAGE: DWS class method
% [out,metric,badAverages] = obj.rmbadaverages_FidA(in,nsd,domain,plot_figure);
%
% DESCRIPTION:
% Removes motion corrupted averages from a dataset containing multiple
% averages.  Bad averages are identified by calculating a 'likeness' metric
% for each average.  This is done by subtracting each average from the
% median of the averages, and then calculating the root mean squared of
% this difference spectrum.  Averages whose likeness metrics are greater
% than 'nsd' above the mean are discarded.
%
% INPUTS:
% obj         = DWS class object
% in          = input data in matlab structure format
% nsd         = number of standard deviations to use a rejection threshold
% domain      = domain in which to perform calculations ('t' or 'f')
% plot_figure = Display figures in the function
%
% OUTPUTS:
% out         = Output dataset following removal of motion corrupted averages.
% metric      = Vector of unlikeness metrics corresponding to all input
%               averages. 
% badAverages = Indices of the averages that were removed. 

function [out,metric,badAverages] = rmbadaverages_FidA(~,in,nsd,domain,plot_figure)
if in.flags.averaged
    error('ERROR:  Averaging has already been performed!  Aborting!');
end

if ~in.flags.addedrcvrs
    error('ERROR:  Receivers should be combined first!  Aborting!');
end

if nargin < 5
    plot_figure = 0;
    if nargin < 4
        domain = 't';
        if nargin < 3
            nsd = 3;
        end
    end
end

%first, make a metric by subtracting all averages from the first average, 
%and then taking the sum of all all the spectral points.  
if in.dims.subSpecs > 0
    SS = in.sz(in.dims.subSpecs);
else
    SS = 1;
end
if domain == 't' || domain == 'T'
    infilt = in;
    tmax = 0.4;
    trange = (infilt.t >= 0 & infilt.t <= tmax);
elseif domain == 'f' || domain == 'F'
    filt = 10;
    infilt = op_filter(in,filt);
end
inavg = op_median(infilt);

metric = zeros(in.sz(in.dims.averages),SS);
for n = 1:in.sz(in.dims.averages)
    for m = 1:SS
        if domain == 't' || domain == 'T'
            metric(n,m) = sum((real(infilt.fids(trange,n,m))-(real(inavg.fids(trange,m)))).^2);
        elseif domain == 'f' || domain == 'F'
            metric(n,m) = sum((real(infilt.specs(:,n,m))-(real(inavg.specs(:,m)))).^2);
        end
    end
end

%find the average and standard deviation of the metric
avg = mean(metric);
stdev = std(metric);

%Now z-transform the metric so that it is centered about zero, and they
%have a standard deviation of 1.0.  
zmetric = (metric-avg) ./ stdev;

P = zeros(SS,3);
for m = 1:SS
    P(m,:) = polyfit((1:in.sz(in.dims.averages))',zmetric(:,m),2);
    if plot_figure
        figure('position',[0 (m-1)*500 560 420]);
        plot(1:in.sz(in.dims.averages),zmetric(:,m),'.',...
            1:in.sz(in.dims.averages),polyval(P(m,:),1:in.sz(in.dims.averages)),...
            1:in.sz(in.dims.averages),(polyval(P(m,:),1:in.sz(in.dims.averages))' + nsd),':');
        xlabel('Scan Number');
        ylabel('Unlikeness Metric (z-score)');
        title('Metric for rejection of motion corrupted scans');
    end
end

%Now make a mask that represents the locations of the averages 
%whose metric values are more than nsd standard deviations away from the 
%mean metric value.
mask = zeros(in.sz(in.dims.averages),SS);
for n = 1:SS
    %mask(:,n) = metric(:,n)>(avg(n)+(nsd*stdev(n))) | metric(:,n)<(avg(n)-(nsd*stdev(n)));
    %mask(:,n) = metric(:,n)>(polyval(P(n,:),[1:in.sz(in.dims.averages)])'+(nsd*stdev(n))) | metric(:,n)<(polyval(P(n,:),[1:in.sz(in.dims.averages)])'-(nsd*stdev(n)));
    mask(:,n) = zmetric(:,n)>(polyval(P(n,:),1:in.sz(in.dims.averages))'+nsd);
end

%Unfortunately, if one average is corrupted, then all of the subspecs
%corresponding to that average have to be thrown away.  Therefore, take the
%minimum intensity projection along the subspecs dimension to find out
%which averages contain at least one corrupted subspec:
if size(mask,2) > 1
    mask = sum(mask,2)' > 0;
end

%now the corrupted and uncorrupted average numbers are given by:
badAverages = find(mask);
goodAverages = find(~mask);

%make a new fids array containing only good averages
fids = in.fids(:,goodAverages,:,:);
fids = conj(fids);
%%re-calculate Specs using fft
specs = fftshift(ifft(fids,[],in.dims.t),in.dims.t);

%re-calculate the sz variable
sz = size(fids);

%FILLING IN DATA STRUCTURE
out = in;
out.fids = fids;
out.specs = specs;
out.sz = sz;
out.averages = length(goodAverages) * in.rawSubspecs;

%FILLING IN THE FLAGS
out.flags = in.flags;
out.flags.writtentostruct = 1;
end

%% Additional functions

%% op_filter.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% [out,lor]=op_filter(in,lb);
% 
% DESCRIPTION:
% Perform line broadening by multiplying the time domain signal by an
% exponential decay function.  
% 
% INPUTS:
% in     = input data in matlab structure format.
% lb     = line broadening factor in Hz.
%
% OUTPUTS:
% out    = Output following alignment of averages.  
% lor    = Exponential (time domain) filter envelope that was applied.

function [out,lor] = op_filter(in,lb)

if lb == 0
    out=in;
else    
    if in.flags.filtered
        %cont=input('WARNING:  Line Broadening has already been performed!  Continue anyway?  (y or n)','s');
        cont='y';
        if cont=='y'
            %continue;
        else
            error('STOPPING');
        end
    end
    fids = in.fids;
    
    t2 = 1/(pi*lb);
    
    %Create an exponential decay (lorentzian filter):
    lor = exp(-in.t/t2);
    %plot(in.t,lor);
    
    %first make a bunch of vectors of ones that are the same lengths as each of
    %the dimensions of the data.  Store them in a cell array for ease of use.
    p = cell(length(in.sz),1);
    for n = 1:length(in.sz)
        p{n} = ones(in.sz(n),1);
    end
    
    %Now, now take the lorentzian filter vector that we made earlier (lor) and use it
    %to populate a filter array that has the same dimensions as the data.  To
    %do this, we have to use the ndgrid function, which is essentially the same
    %as the meshgrid function, except in multiple dimensions.  b, c, and d are
    %dummy variables and are not used.
    if length(in.sz) == 1
        fil = lor;
    end
    if length(in.sz) == 2
        [fil,~] = ndgrid(lor,p{2});
    end
    if length(in.sz) == 3
        [fil,~,~] = ndgrid(lor,p{2},p{3});
    end
    if length(in.sz) == 4
        [fil,~,~,~] = ndgrid(lor,p{2},p{3},p{4});
    end
    if length(in.sz) == 5
        [fil,~,~,~,~] = ndgrid(lor,p{2},p{3},p{4},p{5});
    end
    
    %Now multiply the data by the filter array.
    fids = fids.*fil;
  
    %re-calculate Specs using fft
    specs = fftshift(ifft(fids,[],in.dims.t),in.dims.t);
    
    %FILLING IN DATA STRUCTURE
    out = in;
    out.fids = fids;
    out.specs = specs;
    
    %FILLING IN THE FLAGS
    out.flags = in.flags;
    out.flags.writtentostruct = 1;
    out.flags.filtered = 1;
end
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