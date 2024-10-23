%%  SNR_map_comp.m 

% Copyright All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, MRS4Brain research group @ CIBM MRI EPFL AIT, 2024
% See the LICENSE.TXT file for more details.

% Brayan Alves, CIBM - MRS4Brain group, 2022
% Guillaume Briand, CIBM - MRS4Brain group, 2023
% 
% USAGE : MRSI Class public method 
% msg = obj.SNR_map_comp(save_figure)
% 
% DESCRIPTION :
% Compute the SNR map for the MRSI data acquired
%
% INPUTS :
% obj           = MRSI Class object with properties and methods
% save_figure   = T/F to save the figures in the data_folder
%
% OUTPUT :
% msg           = Error message
function msg = SNR_map_comp(obj,save_figure)
% Input and output initialization
msg = {''};
if nargin < 2
    save_figure = true;
end
try
    MatSize = obj.acq_params.matrix_sz;
    
    maxSNR_ind_column = zeros(1,obj.Nslices);
    maxSNR_ind_line = zeros(1,obj.Nslices);
    
    obj.SNR_map = zeros(obj.Nslices,MatSize(1),MatSize(2));
    obj.avg_SNR = zeros(obj.Nslices,1);
    
    for ii = 1:obj.Nslices
        slice_number = obj.Slices_number(ii); % On which slice we are
        if obj.Lipsup
            slice_tkk = squeeze(obj.HSVD_lipsup_fid_tkkn(:,:,:,ii)); % CHANGE XNUCLEI
        else
            slice_tkk = squeeze(obj.HSVD_fid_tkkn(:,:,:,ii)); % CHANGE XNUCLEI
        end
        ref_tkk = squeeze(obj.ref_mat_tkkn(:,:,:,ii)); 
        slice_frr = fftshift(fft(ifft(ifft(slice_tkk,[],2),[],3),[],1),1);
        ref_frr = fftshift(fft(ifft(ifft(ref_tkk,[],2),[],3),[],1),1);
        
        % CHANGE HERE FOR X NUCLEI %!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%
        ppm_range_NAA = logical((obj.acq_params.scale_ppm < 2.1).*(obj.acq_params.scale_ppm > 1.9));
        ppm_range_wat = logical((obj.acq_params.scale_ppm < 4.8).*(obj.acq_params.scale_ppm > 4.6));

        map_data = squeeze(max(abs(slice_frr(ppm_range_wat,:,:)),[],1));
        map_ref = squeeze(max(abs(ref_frr(ppm_range_wat,:,:)),[],1));
        max_signal = squeeze(max(abs(slice_frr(ppm_range_NAA,:,:)),[],1));
    
        SNR_data = squeeze(max_signal./squeeze(mean(abs(slice_frr(end-100:end,:,:)),1)));
        
        obj.SNR_map(ii,:,:) = SNR_data;
        ratio_dat_ref = (map_data./map_ref)*100;
        ratio_res_ref = (max_signal./map_data)*100;
        
        total = 0;
        tot_SNR = 0;
        tot_res_ref = 0;
        count = 0;
        for x = 1:MatSize(1)
            for y = 1:MatSize(2)
                if(ratio_dat_ref(x,y) > 100)
                    ratio_dat_ref(x,y) = 100;
                end
                total = total + ratio_dat_ref(x,y) * obj.Brain_mask(ii,x,y);
                tot_SNR = tot_SNR + SNR_data(x,y) * obj.Brain_mask(ii,x,y);
                tot_res_ref = tot_res_ref + ratio_res_ref(x,y) * obj.Brain_mask(ii,x,y);
                count = count + obj.Brain_mask(ii,x,y);
            end
        end
        
        Average = total/count;
        Average_SNR = tot_SNR/count;
        Average_Res_ref = tot_res_ref/count;

        SNR_pert_data = SNR_data./(sqrt(obj.acq_params.acquisition_time));
        Average_SNR_pert = Average_SNR./(sqrt(obj.acq_params.acquisition_time));

        obj.SNR_map(ii,:,:) = SNR_data;
        obj.avg_SNR(ii) = Average_SNR;
    
        [~, maxSNR_ind_column(ii)] = max(max(SNR_data .* squeeze(obj.Brain_mask(ii,:,:))));
        [~, ind_line_prime] = max(SNR_data .* squeeze(obj.Brain_mask(ii,:,:)));
        maxSNR_ind_line(ii) = ind_line_prime(maxSNR_ind_column(ii));
                

        if save_figure
            % Plot
            fid_1 = figure('Visible','off');
            imagesc(ratio_res_ref(end:-1:1,:) .* squeeze(obj.Brain_mask(ii,end:-1:1,:)))
            title(['Signal-to-reference Ratio (in %) / Slice ' ...
                num2str(slice_number) ' / Mean = ' num2str(Average_Res_ref)])
            clim([0 100])
            colorbar;

            fid_2 = figure('Visible','off');
            imagesc(ratio_dat_ref(end:-1:1,:) .* squeeze(obj.Brain_mask(ii,end:-1:1,:)))
            title(['Reference Amplitude AFTER HSVD (Ratio in %) / Slice ' ...
                num2str(slice_number) ' / Mean = ' num2str(Average)])
            clim([0 2]);
            colorbar;

            fid_3 = figure('Visible','off');
            imagesc(SNR_data .* squeeze(obj.Brain_mask(ii,end:-1:1,:)))
            title(['SNR Map / Slice ' num2str(slice_number) ' / Mean = ' num2str(Average_SNR)])
            clim([6 15])
            colorbar;
            set(findall(gcf,'-property','FontSize'),'FontSize',18)
            set(findall(gcf,'-property','FontWeight'),'FontWeight','bold')

            fid_4 = figure('Visible','off');
            imagesc(SNR_pert_data(end:-1:1,:) .* squeeze(obj.Brain_mask(ii,end:-1:1,:)))
            title(['SNR per time / Slice ' num2str(ii) ' / Mean = ' num2str(Average_SNR_pert)])
            clim([0 0.55])
            colorbar;
            set(findall(gcf,'-property','FontSize'),'FontSize',18)
            set(findall(gcf,'-property','FontWeight'),'FontWeight','bold')
            
            % Save
            if(~exist(fullfile(obj.data_folder,'SNR_maps'),"dir"))
                mkdir(obj.data_folder,'SNR_maps');
            end
            set(fid_3,'CreateFcn','set(gcbo,''Visible'',''on'')'); 
            filename = fullfile(obj.data_folder, 'SNR_maps',...
                ['E' num2str(obj.metab_expnb) '_SNRMap_Slice' num2str(slice_number)]);
            saveas(fid_3,filename,'fig');
            saveas(fid_3,filename,'png');
            set(fid_4,'CreateFcn','set(gcbo,''Visible'',''on'')'); 
            filename = fullfile(obj.data_folder, 'SNR_maps',...
                ['E' num2str(obj.metab_expnb) '_SNRpertMap_Slice' num2str(slice_number)]);
            saveas(fid_4,filename,'fig');
            saveas(fid_4,filename,'png');

            close([fid_1,fid_2,fid_3,fid_4]);
        end
    
    end
catch ME
    msg = {'SNR map calculation error',ME.message};
end
end