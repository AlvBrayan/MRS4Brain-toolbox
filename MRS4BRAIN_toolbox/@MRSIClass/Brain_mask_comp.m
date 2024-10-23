%%  Brain_mask_comp.m 

% Copyright All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, MRS4Brain research group @ CIBM MRI EPFL AIT, 2024
% See the LICENSE.TXT file for more details.

% Brayan Alves, CIBM - MRS4Brain group, 2022
% Guillaume Briand, CIBM - MRS4Brain group, 2023
% 
% USAGE : MRSI Class public method 
% msg = obj.Brain_mask_comp(brain_mask,save_figure)
% 
% DESCRIPTION :
% Compute the Brain mask in the MRSI resolution from a brain mask with
% different resolution
%
% INPUTS :
% obj           = MRSI Class object with properties and methods
% brain_mask    = MRI Brain mask computed by registration with atlas
% save_figure   = T/F to save the figures in the data_folder
%
% OUTPUT :
% msg           = Error message
function msg = Brain_mask_comp(obj,brain_mask,save_figure)
% Input and Output initialization
msg = {''};
brain_mask_exist = 1;
if nargin < 3
    save_figure = true;
    if nargin < 2
        brain_mask_exist = 0;
    end
end

try
    MatSize = obj.acq_params.matrix_sz; % MRSI matrix size
    % Brain mask and Power mask Initialization
    obj.Brain_mask = ones(obj.Nslices,MatSize(1),MatSize(2));
    obj.Power_map = ones(obj.Nslices,MatSize(1),MatSize(2));
    
    for ii = 1:obj.Nslices
    %% Power map computation

        slice_number = obj.Slices_number(ii); % On which slice we are
        if ~isempty(obj.ref_mat_tkkn)
            power_on = true;
            ft_ref = fftshift(fft(ifft(ifft( ...
                squeeze(obj.ref_mat_tkkn(:,:,:,ii)),[],2),[],3),[],1),1); % Fourier Transform of the matrix signal
            [~,Nx,Ny] = size(squeeze(obj.ref_mat_tkkn(:,:,:,ii)));
            power = squeeze(sum(ft_ref.*conj(ft_ref),1));
            power = reshape(power,Nx,Ny);
            obj.Power_map(ii,:,:) = power;
        else
            power_on = false;
        end

    %% Brain mask computation
        
        if brain_mask_exist
            [Nx,Ny,~] = size(brain_mask);
            min_slice = slice_number - obj.Slice_range;
            max_slice = slice_number + obj.Slice_range;
            slice = rot90(squeeze(sum(double(brain_mask(:,:,min_slice:max_slice)),3)),-1);
            slice = double(logical(slice));
            super_image = kron(slice,ones(MatSize(1),MatSize(2)));
            R = reshape(super_image, Nx, MatSize(1), Ny, MatSize(2));
            S = sum(sum(R,1),3)/(Nx * Ny);
            final_slice = reshape(S, MatSize(1), MatSize(2));
            obj.Brain_mask(ii,:,:) = (final_slice > 0.8);
        else
            obj.Brain_mask(ii,:,:) = obj.Power_map(ii,:,:);
        end
    
    %% Figure plot and save
       
        if save_figure        
            % Plot
            f = figure('Visible','off');
            imagesc(squeeze(obj.Brain_mask(ii,:,:)))
            title(['NSlice = ' num2str(ii)])
            colorbar
            
            if power_on
                g = figure('Visible','off');
                imagesc(power)
                title(['NSlice = ' num2str(ii)])
                colorbar
            end
            
            % Save
            if(~exist(fullfile(obj.data_folder,'PowerMaps'),"dir"))
                mkdir(obj.data_folder,'PowerMaps');
            end
            set(f,'CreateFcn','set(gcbo,''Visible'',''on'')'); 
            filename = fullfile(obj.data_folder, 'PowerMaps',...
                ['E' num2str(obj.ref_expnb) '_BrainMask_Slice' num2str(ii)]);
            saveas(f,filename,'fig');
            saveas(f,filename,'png');
            if power_on
                set(g,'CreateFcn','set(gcbo,''Visible'',''on'')');
                filename = fullfile(obj.data_folder, 'PowerMaps', ...
                    ['E' num2str(obj.ref_expnb) '_PowerMap_Slice' num2str(ii)]);
                saveas(g,filename,'fig');
                saveas(g,filename,'png');
            end

            close(f);
            if power_on
                close(g);
            end
        end
    end
catch ME
    msg = {'Brain mask calculation error',ME.message};
end

obj.Final_mask = obj.Brain_mask; % Final mask as the Brain mask
% Save the brain mask in the result folder
brain_mask_filename = ['brain_mask_' num2str(MatSize(1)) '_' num2str(MatSize(2)) '.mat'];
if ~isfolder(fullfile(obj.results_folder,'Registration'))
    mkdir(fullfile(obj.results_folder,'Registration'))
end
if ~isfile(fullfile(obj.results_folder,'Registration',brain_mask_filename))
    Brain_mask = obj.Brain_mask;
    save(fullfile(obj.results_folder,'Registration',brain_mask_filename),'Brain_mask')
end
end