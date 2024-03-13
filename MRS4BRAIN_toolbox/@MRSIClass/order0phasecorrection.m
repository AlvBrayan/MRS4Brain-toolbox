% Copyright All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, MRS4Brain research group @ CIBM MRI EPFL AIT, 2024
% See the LICENSE.TXT file for more details.

function [corr_spect,phase0,integ_val] = order0phasecorrection(~,spect,ppm,ppm_max,ppm_min,err,eta)

ppm_range = ppm(logical((ppm < ppm_max).*(ppm > ppm_min)));
coeff_angle = 2*pi/360;

% First, let's try a gradient descent approach
% the loss function is defined as
% -integral(real(spect*phase_corr),ppm_min,ppm_max)
% Thus the derivative becomes
% -coeff_angle*integral(real(-sqrt(-1)*spect*phase_corr),ppm_min,ppm_max)

spect_range = real(-sqrt(-1)*spect(logical((ppm < ppm_max).*(ppm > ppm_min))));
init_grad = -coeff_angle*trapezoid_intergal(spect_range,ppm_range);

theta_init = 0;
theta_plus1 = 0 - eta*init_grad;

new_spect = spect;
N_inter=0;
while (abs(theta_plus1-theta_init) > err)
    new_eta = eta;
    % disp(['We are at the ' num2str(N_inter+1) ' iteration now !'])
    phase_corr = exp(-sqrt(-1)*coeff_angle*theta_plus1);
    new_spect = spect.*phase_corr;
    spect_range = real(-sqrt(-1)*new_spect(logical((ppm < ppm_max).*(ppm > ppm_min))));
    grad = -coeff_angle*trapezoid_intergal(spect_range,ppm_range);
    % while(abs(eta*grad)>180)
    %     new_eta = 90/grad;
    % end

    theta_init = theta_plus1;
    theta_plus1 = theta_plus1 - new_eta*grad;
    % if(theta_plus1 < -180)
    %     theta_plus1 = mod(theta_plus1,360);
    % end
    % if(theta_plus1 > 180)
    %     theta_plus1 = mod(theta_plus1,-360);
    % end
    N_inter = N_inter +1;
end

% disp([num2str(N_inter) ' iterations were needed !'])
phase0 = theta_init;
corr_spect = new_spect;
integ_val = trapezoid_intergal(real(new_spect(logical((ppm < ppm_max).*(ppm > ppm_min)))),ppm_range);

%% Additional functions 

function val_integral = trapezoid_intergal(array_f,scale)

dx = abs(scale(1:end-1)-scale(2:end));

val_integral = 0;

for kk=1:length(dx)
    val_integral = val_integral + dx(kk)*(array_f(kk)+array_f(kk+1))*0.5;
end
