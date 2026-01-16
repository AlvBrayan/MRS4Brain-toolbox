%%  save_data.m 

% Copyright All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, MRS4Brain research group @ CIBM MRI EPFL AIT, 2024
% See the LICENSE.TXT file for more details.

% Brayan Alves, CIBM - MRS4Brain group, 2022
% Guillaume Briand, CIBM - MRS4Brain group, 2023
% 
% USAGE : MRSI Class public method 
% msg = obj.save_data()
% 
% DESCRIPTION :
% Save current MRSI voxels data as RAW files to be then quantified
% with LCModel
%
% INPUTS :
% obj       = MRSI Class object with properties and methods
%
% OUTPUT :
% msg       = Error message
function msg = save_data(obj)
% save_data.m to save FIDs as RAW files to be quantified with LCModel

    keyParams.data_folder       = obj.data_folder;
    keyParams.results_folder    = obj.results_folder;
    keyParams.metab_expnb       = obj.metab_expnb;
    keyParams.ref_expnb         = obj.ref_expnb;
    keyParams.Slices_number     = obj.Slices_number;
    keyParams.reco_expnb        = obj.reco_expnb;
    keyParams.mrsi_params       = obj.mrsi_params;
    keyParams.LipsupOnOff       = obj.Lipsup;
    keyParams.FillgapOnOff      = obj.Fillgaps;

    if ~isempty(dir(fullfile(obj.results_folder,'Nbasis.mat'))) 
        load([obj.results_folder,'\Nbasis.mat']); delete([obj.results_folder,'\Nbasis.mat']); 
        keyParams.NLipid        = Nbasis;
    else 
        keyParams.NLipid        = [];
    end
    save([keyParams.results_folder,'\Processing_Pars.mat'],'keyParams');

msg = {''};

try
    % Save data in RAW files
    if ~isfolder(fullfile(obj.results_folder))
        mkdir(fullfile(obj.results_folder));
    end
    for ii = 1:obj.Nslices % goes through all the slices
        if and(obj.Fillgaps,obj.Lipsup)
            metab_data_tkk = obj.HSVD_lipsup_filled_fid_tkkn(1:obj.acq_params.np_met,:,:,ii); % WITH FILLGAPS & LIPSUP
        elseif and(obj.Fillgaps,~obj.Lipsup)
            metab_data_tkk = obj.HSVD_lipsup_filled_fid_tkkn(1:obj.acq_params.np_met,:,:,ii); % WITH FILLGAPS & LIPSUP
        elseif and(obj.Lipsup,~obj.Fillgaps)
            metab_data_tkk = obj.HSVD_lipsup_fid_tkkn(:,:,:,ii); % WITH LIPSUP
        else
            metab_data_tkk = obj.HSVD_fid_tkkn(:,:,:,ii); % WITHOUT LIPSUP
        end       
        metab_data_trr = ifft(ifft(metab_data_tkk,[],2),[],3); % Real space
        reference_signal = ifft(ifft(obj.ref_mat_tkkn(:,:,:,ii),[],2),[],3); % Real space
        output_name = ['Slice_N' num2str(ii)]; % MRSI slice number
        folder_name = fullfile(obj.results_folder,output_name);
        if isfolder(folder_name)
            rmdir(folder_name,'s');
        end
        mkdir(folder_name);
        data_folder = fullfile(folder_name,'Data');
        if ~isfolder(data_folder)
            mkdir(data_folder);
        end
        msg = obj.data2RAW(metab_data_trr,reference_signal,output_name,ii); % call data2RAW function
        if any(~cellfun(@isempty,msg))
            return;
        end
    end

    % Save the Final mask
    Final_mask = obj.Final_mask;
    mkdir(obj.results_folder,'Mask')
    save(fullfile(obj.results_folder,'Mask',['Final_mask_',output_name,'.mat']),"Final_mask")

catch ME
    msg = {'Error saving the data',ME.message};
end
end