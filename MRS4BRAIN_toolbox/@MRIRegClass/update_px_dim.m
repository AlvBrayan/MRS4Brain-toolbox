%% update_px_dim.m 

% Copyright All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, MRS4Brain research group @ CIBM MRI EPFL AIT, 2024
% See the LICENSE.TXT file for more details.

% Guillaume Briand, CIBM - MRS4Brain group, 2023
% 
% USAGE : MRIReg Class public method 
% msg = obj.update_px_dim()
% 
% DESCRIPTION :
% Update the pixel dimension by applying a gain coefficient. Used to
% correspond to human brain dimensions
%
% INPUTS :
% obj       = MRIReg Class object with properties and methods
% G         = Gain coefficient of the pixel dimension 
%
% OUTPUT :
% msg       = Error message
function update_px_dim(obj,G,head_prone)
% Change the pixel dimension to correspond to human brain dimensions
if nargin < 3
    head_prone = 0;
end
if nargin < 2
    G = 10;
end

info = niftiinfo(obj.Nifti_image_filename);
% Open the image slices
V = niftiread(info);

if head_prone
    V = V(:,end:-1:1,:);
end

info.PixelDimensions = G*info.PixelDimensions;

info.raw.pixdim(1,2:4) = info.raw.pixdim(1,2:4)*G;
info.raw.srow_x(1) = info.raw.srow_x(1)*G;
info.raw.srow_y(2) = info.raw.srow_x(2)*G;
info.raw.srow_z(3) = info.raw.srow_x(3)*G;

info.Transform.T(1,1) = info.Transform.T(1,1)*G;
info.Transform.T(2,2) = info.Transform.T(2,2)*G;
info.Transform.T(3,3) = info.Transform.T(3,3)*G;

if(~exist(obj.registration_folder,"dir"))
    mkdir(obj.registration_folder);
end

% Save the Nifti with header and image in compressed file with gunzip
if isfile(fullfile(obj.registration_folder,'image_study_xG.nii.gz'))
    info2 = niftiinfo(fullfile(obj.registration_folder,'image_study_xG.nii.gz'));
    % Open the image slices
    V2 = niftiread(info2);
    if size(V) ~= size(V2)
        niftiwrite(V,fullfile(obj.registration_folder,'image_study_xG'),info,'Compressed',true);
    end
else
    niftiwrite(V,fullfile(obj.registration_folder,'image_study_xG'),info,'Compressed',true);
end
obj.Nifti_image_xG_filename = fullfile(obj.registration_folder,'image_study_xG.nii.gz');
end