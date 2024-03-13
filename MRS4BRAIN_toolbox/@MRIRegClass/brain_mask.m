%%  brain_mask.m 

% Copyright All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, MRS4Brain research group @ CIBM MRI EPFL AIT, 2024
% See the LICENSE.TXT file for more details.

% Guillaume Briand, CIBM - MRS4Brain group, 2023
% 
% USAGE : MRIReg Class public method 
% msg = obj.brain_mask()
% 
% DESCRIPTION :
% Obtain the Brain mask from an exisiting Nifti brain mask file
%
% INPUTS :
% obj       = MRIReg Class object with properties and methods
%
% OUTPUT :
% msg       = Error message
function msg = brain_mask(obj)
msg = {''};
brain_mask_fn = fullfile(obj.registration_folder,'brain_mask.nii.gz');
if isfile(brain_mask_fn)

    Brain_mask = logical(niftiread(brain_mask_fn));
    obj.Brain_mask = Brain_mask;

    if ~isfolder(fullfile(obj.result_folder,'Registration'))
        mkdir(fullfile(obj.result_folder,'Registration'))
    end

    if ~isfile(fullfile(obj.result_folder,'Registration','brain_mask_full.mat'))
        save(fullfile(obj.result_folder,'Registration','brain_mask_full.mat'),'Brain_mask')
    end
else
    msg = {'No existing file'};
end
end