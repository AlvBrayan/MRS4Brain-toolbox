% Copyright All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, MRS4Brain research group @ CIBM MRI EPFL AIT, 2024
% See the LICENSE.TXT file for more details.

function msg = Linewidth_map_comp(obj,save_figure)
% Input and Output initialization
msg = {''};
if nargin < 2
    save_figure = true;
end

try
    MatSize = obj.acq_params.matrix_sz; % MRSI Matrix Size
    fmax = obj.acq_params.spectralwidth/2;
    f = fmax:-2 * fmax/(obj.acq_params.np_met - 1):-fmax;
    
    % Linewidth map initialization
    obj.Linewidth_map = zeros(obj.Nslices,MatSize(1),MatSize(2));
    
    for ii = 1:obj.Nslices
        FWHM_map = zeros(MatSize(1),MatSize(2)); % Linewidth map init
        mean_FWHM = 0;
        ft_ref_temp = fftshift(fft(ifft(ifft(squeeze(obj.ref_mat_tkkn(:,:,:,ii)),[],2),[],3),[],1),1);

        % 0 order-phase correction applied
        [ft_ref,~] = obj.MRSI_0orderphasecorrection_water(ft_ref_temp, ...
            obj.acq_params.scale_ppm,4.78,4.62,squeeze(obj.Brain_mask(ii,:,:)),0.005,1E-6);
        for x = 1:MatSize(1)
            for y = 1:MatSize(2)
                %First, we find the maximum value of the absolute spectrum
                spect = squeeze(ft_ref(:,x,y));
                [val,~] = max(abs(real(spect)));
    
                %Then, we find the indexes corresponding to a value close to the half the maximum
                int_spect = abs(real(spect)) - ones(size(spect,1),1)*val*0.5;
                list = find(int_spect>=0);
                left_ind = list(1);
                right_ind = list(end);

                if(left_ind==1)
                    left_ind = 2;
                end
                if(right_ind==length(spect))
                    right_ind = length(spect)-1;
                end
    
                %Then, we translate those indexes into frequencies (Hz)
                left_Hz = f(left_ind);
                left_Hz_under = f(left_ind-1);
                right_Hz = f(right_ind);
                right_Hz_under = f(right_ind+1);
    
                %Finally, we use a linear interpolation to find the expected
                %frequency where we are at half maximum
    
                left_alpha = int_spect(left_ind-1)/(int_spect(left_ind)-int_spect(left_ind-1));
                left_ans = (1+left_alpha)*left_Hz_under - left_alpha*left_Hz;
    
                right_alpha = int_spect(right_ind)/(int_spect(right_ind+1)-int_spect(right_ind));
                right_ans = (1+right_alpha)*right_Hz - right_alpha*right_Hz_under;
    
                FWHM_map(x,y) = (left_ans-right_ans)*obj.Brain_mask(ii,x,y); % Linewidth map for slice ii
                mean_FWHM = mean_FWHM+(left_ans-right_ans)*obj.Brain_mask(ii,x,y);
            end
        end
    
        mean_FWHM = mean_FWHM/(sum(sum(obj.Brain_mask(ii,x,y))));
    
        ppm_map = FWHM_map/obj.acq_params.resfreq;
        mean_ppm = mean_FWHM/obj.acq_params.resfreq;
    
        obj.Linewidth_map(ii,:,:) = FWHM_map;
    
        
        % Save figure on the data folder
        if save_figure
            % Plot
            f = figure('Visible','off');
            imagesc(FWHM_map)
            title(['Linewidth map (mean = ' num2str(mean_FWHM) ' Hz) : NSlice = ' num2str(ii)])
            clim([0 60])
            colorbar

            g = figure('Visible','off');
            imagesc(ppm_map)
            title(['Linewidth map (mean = ' num2str(mean_ppm) ' ppm) : NSlice = ' num2str(ii)])
            clim([0 0.1])
            colorbar
            % Save
            if(~exist(fullfile(obj.data_folder,'Linewidths'),'dir'))
                mkdir(obj.data_folder,'Linewidths');
            end
            set(f,'CreateFcn','set(gcbo,''Visible'',''on'')'); 
            filename = fullfile(obj.data_folder,'Linewidths', ...
                ['E' num2str(obj.ref_expnb) '_LinewidthHz_Slice' num2str(ii)]);
            saveas(f,filename,'fig');
            saveas(f,filename,'png');
            set(g,'CreateFcn','set(gcbo,''Visible'',''on'')'); 
            filename = fullfile(obj.data_folder,'Linewidths', ...
                ['E' num2str(obj.ref_expnb) '_Linewidthppm_Slice' num2str(ii)]);
            saveas(g,filename,'fig');
            saveas(g,filename,'png');

            close([f,g]);
        end
        
    end
catch ME
    msg = {'Linewidth map calculus error :',ME.message};
end
end