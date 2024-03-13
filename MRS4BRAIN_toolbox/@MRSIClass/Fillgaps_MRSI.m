% Copyright All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, MRS4Brain research group @ CIBM MRI EPFL AIT, 2024
% See the LICENSE.TXT file for more details.

function msg = Fillgaps_MRSI(obj)
%FUNCTION ADAPTED FROM ANTOINE'S CODE (MRSI_CSSENSELR_Recon.m) to fill the
%gaps due to the acquisition delay
msg = {''};
try 
    obj.Fillgaps = true;
    obj.HSVD_lipsup_filled_fid_tkkn = obj.HSVD_lipsup_fid_tkkn;
    NbMissingPoints = round(obj.MRSI_parameters.samplerate * obj.MRSI_parameters.acq_delay);

    for ii = 1:obj.Nslices
        % current slice
        if obj.Lipsup
            slice_tkk = squeeze(obj.HSVD_lipsup_fid_tkkn(:,:,:,ii));
        else
            slice_tkk = squeeze(obj.HSVD_fid_tkkn(:,:,:,ii));
        end
%         VSize = size(slice_tkk,1) + NbMissingPoints;
        
        Size_data = size(slice_tkk);
        NaNData = [NaN([NbMissingPoints,Size_data(2)*Size_data(3)]); ...
            reshape(slice_tkk,Size_data(1),[])];
        slice_fg_tkk = fft(reshape(fillgaps(double(NaNData)),...
            [size(NaNData,1) Size_data(2:3)]),[],1); 
        % Using matlab fillgaps function (Signal processing toolbox)
        obj.HSVD_lipsup_filled_fid_tkkn(:,:,:,ii) = ifft(slice_fg_tkk,[],1);
    end
catch ME
    msg = {'Fillgaps error :',ME.message};
    obj.Fillgaps = false;
end
end