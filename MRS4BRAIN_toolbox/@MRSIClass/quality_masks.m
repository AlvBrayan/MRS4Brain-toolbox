% Copyright All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, MRS4Brain research group @ CIBM MRI EPFL AIT, 2024
% See the LICENSE.TXT file for more details.

function msg = quality_masks(obj,save_figure)
% Input and Output initialization
msg = {''};
if nargin < 2
    save_figure = true;
end

try
    for ii = 1:obj.Nslices
        slice_number = obj.Slices_number(ii); % current MRSI slab
        FWHM_mask = (obj.Linewidth_map(ii,:,:) < 60) .* obj.Brain_mask(ii,:,:); % Linewidth mask
        SNR_mask = (obj.SNR_map(ii,:,:) >= 6) .* obj.Brain_mask(ii,:,:); % SNR mask
        int_mask = FWHM_mask;
        vol_mask = sum(sum(int_mask));
        vol_brain = sum(sum(obj.Brain_mask(ii,:,:)));
        
        ratio_voxels = (vol_mask/vol_brain)*100;
        
        tot_mask = obj.Brain_mask(ii,:,:) + FWHM_mask; % Total mask 
        
        if save_figure
            mask_1 = figure('Visible','off');
            imagesc(squeeze(FWHM_mask))
            title('Linewidth mask (<60 Hz)')
            
            mask_2 = figure('Visible','off');
            imagesc(squeeze(SNR_mask))
            title('SNR mask (>6)')
            
            mask_3 = figure('Visible','off');
            imagesc(squeeze(tot_mask))
            title(['Comparison Brain map (', num2str(vol_brain), ...
                ') VS quality map (', num2str(vol_mask), ') / Ratio : ',  ...
                num2str(ratio_voxels), '%'])
            

            if(~exist(fullfile(obj.data_folder,'Masks'),"dir"))
                mkdir(obj.data_folder,'Masks');
            end
            set(mask_1,'CreateFcn','set(gcbo,''Visible'',''on'')'); 
            filename = fullfile(obj.data_folder, 'Masks',...
                ['E' num2str(obj.metab_expnb) '_LinewidthMask_Slice' num2str(slice_number)]);
            saveas(mask_1,filename,'fig');
            saveas(mask_1,filename,'png');
            set(mask_2,'CreateFcn','set(gcbo,''Visible'',''on'')'); 
            filename = fullfile(obj.data_folder, 'Masks',...
                ['E' num2str(obj.metab_expnb) '_SNRMask_Slice' num2str(slice_number)]);
            saveas(mask_2,filename,'fig');
            saveas(mask_2,filename,'png');
            set(mask_3,'CreateFcn','set(gcbo,''Visible'',''on'')'); 
            filename = fullfile(obj.data_folder, 'Masks',...
                ['E' num2str(obj.metab_expnb) '_TotMask_Slice' num2str(slice_number)]);
            saveas(mask_3,filename,'fig');
            saveas(mask_3,filename,'png');

            close([mask_1,mask_2,mask_3]);
        end
        
    end
catch ME
    msg = {'Quality control mask calculation',ME.message};
end
end