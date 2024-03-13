%%  plot_signal.m 

% Copyright All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, MRS4Brain research group @ CIBM MRI EPFL AIT, 2024
% See the LICENSE.TXT file for more details.

% Guillaume Briand, CIBM - MRS4Brain group, 2023
% 
% USAGE : MRSI Class public method 
% [scale_ppm,data_spec] = obj.plot_signal(pos_x,pos_y)
% 
% DESCRIPTION :
% Return ppm scale and data spectrum for current voxel
%
% INPUTS :
% obj       = MRSI Class object with properties and methods
% pos_x     = x position of the voxel
% pos_y     = y position of the voxel
%
% OUTPUT :
% scale_ppm = ppm scale to plot data
% data_spec = Fourier transformed data fids
function plot_signal(obj,pos_x,pos_y)

if nargin < 4
    if nargin < 3
        if nargin < 2
            pos_x = 15;
        end
        pos_y = 15;
    end
end

fmax = obj.acq_params.spectralwidth/2;
% Data ppm scale
f_data = fmax:-2*fmax/(obj.acq_params.np_met - 1):-fmax;
scale_ppm_data = f_data/obj.acq_params.resfreq + obj.acq_params.ppm_ref;
% Reference ppm scale
f_ref = fmax:-2*fmax/(obj.acq_params.np_ref - 1):-fmax;
scale_ppm_ref = f_ref/obj.acq_params.resfreq + obj.acq_params.ppm_ref;

for ii = 1:obj.Nslices
    
    if obj.Lipsup
        slice_tkk = squeeze(obj.HSVD_lipsup_fid_tkkn(:,:,:,ii));
    else
        if isempty(obj.HSVD_fid_tkkn)
            slice_tkk = squeeze(obj.fid_mat_tkkn(:,:,:,ii));
        else
            slice_tkk = squeeze(obj.HSVD_fid_tkkn(:,:,:,ii));
        end
    end
    
    ref_slice_tkk = squeeze(obj.ref_mat_tkkn(:,:,:,ii));

    slice_trr = ifft(ifft(slice_tkk,[],2),[],3);
    ref_slice_trr = ifft(ifft(ref_slice_tkk,[],2),[],3);
    
    data = squeeze(slice_trr(:,pos_y,pos_x));
    ref_data = squeeze(ref_slice_trr(:,pos_y,pos_x));
    
    data_spec = fftshift(fft(data,[],1),1);
    ref_spec = fftshift(fft(ref_data,[],1),1);

    im1 = figure('Visible','on');
    plot(scale_ppm_data,real(data_spec))
    %plot(scale_ppm,real(data_spec))
    %plot(abs(data_spec))
    %ylim([-1 1]*1e5)
    xlim([0 5])
    %plot(f,real(data_spec(end:-1:1)))
    set(gca,'Xdir','reverse')
    title([' Data signal (Slice ' num2str(ii) ') : Pix_x = ',num2str(pos_x),' / Pix_y = ',num2str(pos_y)])
    
    im2 = figure('Visible','on');
    plot(scale_ppm_ref,abs(ref_spec))
    %plot(scale_ppm,abs(ref_spec))
    %plot(abs(ref_spec))
    set(gca,'Xdir','reverse')
    title([' Reference signal (Slice ' num2str(ii) ') : Pix_x = ',num2str(pos_x),' / Pix_y = ',num2str(pos_y)])

    im1.Position(1:2) = [50 ,350];
    im2.Position(1:2) = [1350 ,350];

end
end