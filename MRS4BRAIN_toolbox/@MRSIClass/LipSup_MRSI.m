% Copyright All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, MRS4Brain research group @ CIBM MRI EPFL AIT, 2024
% See the LICENSE.TXT file for more details.

function msg = LipSup_MRSI(obj,PercentThres)
% LipSup_MRSI allows to apply the Lipid suppression method
if nargin < 2
    PercentThres = 0.8;
end

% Output initialization
msg = {''};
obj.HSVD_lipsup_fid_tkkn = obj.HSVD_fid_tkkn;
if isempty(obj.HSVD_fid_tkkn)
    msg = {'Please apply HSVD before Lipid suppression'};
    return
end

try
    obj.Lipsup = true;
    for ii = 1:obj.Nslices
        slice_tkk = obj.HSVD_fid_tkkn(:,:,:,ii); % Current slice metabolite fid
        BrainMap = squeeze(obj.Brain_mask(ii,:,:)); % Brain mask of current slice
        mrsiReconParams = obj.Create_mrsiReconParams_BA(slice_tkk,PercentThres,BrainMap);

        %ADDED THE 30/10/2024
        if ~exist([mrsiReconParams.Log_Dir '/LipidSuppression/']);mkdir([mrsiReconParams.Log_Dir '/LipidSuppression/']);end

        cd([mrsiReconParams.Log_Dir '/LipidSuppression/']);
        
        mrsiReconParams.Log_Dir = [mrsiReconParams.Log_Dir '/LipidSuppression/'];
        %END
        
        % Nbasis
        if(mrsiReconParams.L2SVDparams.PercentThres < 2)
            Nbasis = obj.FindNBasisAllCoil_BA(slice_tkk,mrsiReconParams);
        else
            Nbasis = round(mrsiReconParams.L2SVDparams.PercentThres);
        end

        [newMRSI_tkk,~,~,~,~] = obj.ProjSVDLipidSuppression_BA(slice_tkk,slice_tkk,mrsiReconParams,Nbasis,'LipidComponents');

        obj.HSVD_lipsup_fid_tkkn(:,:,:,ii) = newMRSI_tkk; % update the HSVD lipsup matrix
    end
catch ME
    obj.Lipsup = false;
    msg = {ME.message};
end

end