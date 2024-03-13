%%  image_segmentation.m 

% Copyright All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, MRS4Brain research group @ CIBM MRI EPFL AIT, 2024
% See the LICENSE.TXT file for more details.

% Guillaume Briand, CIBM - MRS4Brain group, 2023
% 
% USAGE : MRIReg Class public method 
% msg = obj.image_segmentation()
% 
% DESCRIPTION :
% Apply image segmentation using transformation found by image registration
%
% INPUTS :
% obj       = MRIReg Class object with properties and methods
%
% OUTPUT :
% msg       = Error message
function msg = image_segmentation(obj)
msg = {''};
try
    delfiles = fullfile(obj.registration_folder,'brain_labels.nii.gz');
    if isfile(delfiles)
        delete(delfiles)
    end

    needed_files_seg = {obj.registered_image_fn,...
        fullfile(obj.tmp_registration_path,'tmp_registration','space_transformation.nii.gz'),...
        fullfile(obj.tmp_registration_path,'tmp_registration','affine_transformation.mat')};
    
    if any(~isfile(needed_files_seg))
        msg = {'Please register the image with the template before processing the segmentation !'};
        return
    end

    % Filenames
    [~, temp_label_fn, ~] = fileparts(obj.template_labels_filename);
    [~, temp_label_fn, ~] = fileparts(temp_label_fn);
    [~, image_fn, ~] = fileparts(obj.Nifti_image_xG_filename);
    [~, image_fn, ~] = fileparts(image_fn);

    % Move files to registration folder
    if ~isfile(fullfile(obj.tmp_registration_path,'tmp_registration',[temp_label_fn '.nii.gz']))
        copyfile(obj.template_labels_filename,fullfile(obj.tmp_registration_path,'tmp_registration'),'f');
    end

    toolbox_folder = pwd;
    cd(obj.tmp_registration_path)
    start_docker = '';
    if strcmp(computer,'MACI64')
        start_docker = '/usr/local/bin/docker run -v "./tmp_registration:/data" --rm antsx/ants';
    elseif strcmp(computer,'PCWIN64')
        start_docker = 'docker run -v ".\tmp_registration:/data" antsx/ants';
    end
    disp([start_docker ' /opt/ants/bin/antsApplyTransforms -d 3' ...
            ' -i /data/' temp_label_fn '.nii.gz' ...
            ' -r /data/' image_fn '_r.nii.gz' ...
            ' -o brain_labels.nii.gz' ...
            ' -t /data/space_transformation.nii.gz' ...
            ' -t /data/affine_transformation.mat' ...
            ' -n NearestNeighbor'])
    [~,~] = system([start_docker ' /opt/ants/bin/antsApplyTransforms -d 3' ...
            ' -i /data/' temp_label_fn '.nii.gz' ...
            ' -r /data/' image_fn '_r.nii.gz' ...
            ' -o brain_labels.nii.gz' ...
            ' -t /data/space_transformation.nii.gz' ...
            ' -t /data/affine_transformation.mat' ...
            ' -n NearestNeighbor'],"-echo");
    cd(toolbox_folder)

    time_elapsed = 0;
    while and(~isfile(fullfile(obj.tmp_registration_path,'tmp_registration','brain_labels.nii.gz')),time_elapsed < 10)
        pause(0.5)
        time_elapsed = time_elapsed + 0.5;
    end
    if time_elapsed >= 10
        msg = {'Error while applying the space and affine transformations with ANTs'};
    end

    copyfile(fullfile(obj.tmp_registration_path,'tmp_registration'),obj.registration_folder,'f');

    rmdir(fullfile(obj.tmp_registration_path,'tmp_registration'),'s')

catch ME
    msg = {ME.message};
        rmdir(fullfile(obj.tmp_registration_path,'tmp_registration'),'s')
end

end