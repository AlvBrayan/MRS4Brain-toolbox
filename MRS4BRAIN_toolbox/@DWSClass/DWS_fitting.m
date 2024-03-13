%%  DWS_fitting.m 

% Copyright All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, MRS4Brain research group @ CIBM MRI EPFL AIT, 2024
% See the LICENSE.TXT file for more details.

% Jessie Mosso, CIBM - MRS4Brain group, LIFMET, 2021
% Guillaume Briand, CIBM - MRS4Brain group, 2023
% 
% USAGE : DWS Class public method 
% out = obj.DWS_fitting()
% 
% DESCRIPTION :
% Convert processed data from study structure to RAW files in order to be
% quantified by LCModel
%
% INPUTS :
% obj       = DWS Class object with properties and methods
%
% OUTPUT :
% msg       = Error message
function DWS_fitting(obj)

% gdiff
idx_met = find([obj.DWS_struct.met_bool]);
b_value = zeros(length(idx_met),1);
for t = 1:length(idx_met)
    params = obj.DWS_struct(t).raw_study.dwsparams;
    b_value(t) = params.bvalue * 10^(-3);
end
obj.b_value = b_value;

% Oriented sticks model
pa_sticks = @(Da,b) sqrt(pi./(4*b.*Da)).*erf(sqrt(b.*Da));
x0_pasticks = 0.45; % Da initial guess

% Kurtosis model
kurtosis = @(b,D,K) exp(-b.*D + 1/6.*b.^2*D.^2.*K);
x0_kurtosis = [0.07 1]; % [D K] initial guess
opts = optimset('Display','off'); % lscurvefit options

rc = 19; %21; % radius of convergence of kurtosis model
excluded_bval = 30; % excluded points above this b-value

Da_fit = NaN*ones(length(obj.Met_names),1);
adc = Da_fit; akc = Da_fit;
for i = 1:length(obj.Met_names) % going through the metabolite list
    b = obj.b_value; % setting the b-value
    [b,idx_b] = sort(b);
    abs_conc = obj.abs_conc_tot(i,idx_b);
    S1 = squeeze(abs_conc)/ squeeze(abs_conc(1));
    index = b < excluded_bval;
    b = b(~isnan(S1.') & index); S1 = S1(~isnan(S1.') & index);
    if ~isempty(b)
        try
            disp(b)
            disp(S1)
            x_sol = lsqcurvefit(@(x,b) pa_sticks(x,b),x0_pasticks,b,S1.',0,1,opts);
            % solving the Dcoeff with the model and the data given
            Da_fit(i) = x_sol; % put the Dcoeff computed in the Dcoeff matrix
        catch
        end
    end
    b = obj.b_value; % setting the b-value
    [b,idx_b] = sort(b);
    index = b < rc;
    abs_conc = obj.abs_conc_tot(i,idx_b);
    S1 = squeeze(abs_conc)/ squeeze(abs_conc(1));
    b = b(~isnan(S1.') & index); S1 = S1(~isnan(S1.') & index);
    if ~isempty(b)
        try
        x_sol = lsqcurvefit(@(x,b) kurtosis(b,x(1),x(2)),x0_kurtosis,b,S1.',[0 0],[1 3],opts);
        adc(i) = x_sol(1);
        akc(i) = x_sol(2);
        catch
        end
    end
end
obj.Da_fit = Da_fit;
obj.adc = adc;
obj.akc = akc;
end