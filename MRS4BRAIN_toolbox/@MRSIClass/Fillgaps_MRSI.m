% Copyright All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, MRS4Brain research group @ CIBM MRI EPFL AIT, 2024
% See the LICENSE.TXT file for more details.

function msg = Fillgaps_MRSI(obj)
%FUNCTION ADAPTED FROM ANTOINE'S CODE (MRSI_CSSENSELR_Recon.m) to fill the
%gaps due to the acquisition delay
msg = {''};
try 
    obj.Fillgaps = true;
    NbMissingPoints = round(obj.acq_params.acq_delay/(1./obj.acq_params.spectralwidth));
    obj.HSVD_lipsup_filled_fid_tkkn = zeros([obj.acq_params.np_met+NbMissingPoints,obj.acq_params.matrix_sz,obj.Nslices]);

    for ii = 1:obj.Nslices
        % current slice
        if obj.Lipsup
            slice_tkk = squeeze(obj.HSVD_lipsup_fid_tkkn(:,:,:,ii));
        else
            slice_tkk = squeeze(obj.HSVD_fid_tkkn(:,:,:,ii));
        end

        % VSize = size(slice_tkk,1) + NbMissingPoints;
        % OLD METHOD
        % Size_data = size(slice_tkk);
        % NaNData = [NaN([NbMissingPoints,Size_data(2)*Size_data(3)]); ...
        %     reshape(slice_tkk,Size_data(1),[])];
        % slice_fg_tkk = fft(reshape(fillgaps(double(NaNData)),...
        %     [size(NaNData,1) Size_data(2:3)]),[],1); 
        % VSize=size(OriginalFullData_trr,1)+NbMissingPoints; % Vsize is the new number of data points: existing ones + missing ones
        % ppm =(-4.7+((1:VSize)/(VSize*obj.acq_params.Bo*obj.acq_params.spectralwidth))); %rescale freq. range considering the new time span
        
        OriginalFullData_trr = ifft(ifft(slice_tkk,[],2),[],3);

        Size_data = size(OriginalFullData_trr); 
        NaNData = [NaN([NbMissingPoints,Size_data(2)*Size_data(3)]); reshape(OriginalFullData_trr,Size_data(1),[])];
        OriginalData_frr = fft(reshape(fillgaps(double(NaNData),300,10),[size(NaNData,1),Size_data(2:3)]),[],1);

        % Using matlab fillgaps function (Signal processing toolbox)
        obj.HSVD_lipsup_filled_fid_tkkn(:,:,:,ii) = ifft(fft(fft(OriginalData_frr,[],2),[],3),[],1);
    end
catch ME
    msg = {'Fillgaps error :',ME.message};
    obj.Fillgaps = false;
end
end