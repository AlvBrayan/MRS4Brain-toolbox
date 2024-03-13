%%  image_registration.m 

% Copyright All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, MRS4Brain research group @ CIBM MRI EPFL AIT, 2024
% See the LICENSE.TXT file for more details.

% Guillaume Briand, CIBM - MRS4Brain group, 2023
% 
% USAGE : MRIReg Class public method 
% msg = obj.image_registration(d)
% 
% DESCRIPTION :
% Apply image registration between MRI rat brain (fixed) and anatomical
% template (moving)
%
% INPUTS :
% obj       = MRIReg Class object with properties and methods
% prog_dbox = MRS4Brain Toolbox progress dialog box
%
% OUTPUT :
% msg       = Error message
function msg = image_registration(obj,d)
msg = {''};
try
    obj.tmp_registration_path = char(java.lang.System.getProperty('user.home'));

    delfiles = {fullfile(obj.registration_folder, 'reg_InverseWarped.nii.gz'), ...
    fullfile(obj.registration_folder, 'reg_0GenericAffine.mat'), ...
    fullfile(obj.registration_folder, 'reg_Warped.nii.gz'), ...
    fullfile(obj.registration_folder, 'reg_1Warp.nii.gz'), ...
    fullfile(obj.registration_folder, 'reg_1InverseWarp.nii.gz')};

    for i = 1:length(delfiles)
        if isfile(delfiles{i})
            delete(delfiles{i})
        end
    end

    if(~exist(fullfile(obj.tmp_registration_path,'tmp_registration'),"dir"))
        mkdir(fullfile(obj.tmp_registration_path,'tmp_registration'));
    end

    % Filenames
    [~, temp_anat_fn, ~] = fileparts(obj.template_anatomical_filename);
    [~, temp_anat_fn, ~] = fileparts(temp_anat_fn);
    [~, temp_bm_fn, ~] = fileparts(obj.template_brain_mask_filename);
    [~, temp_bm_fn, ~] = fileparts(temp_bm_fn);
    [~, temp_label_fn, ~] = fileparts(obj.template_labels_filename);
    [~, temp_label_fn, ~] = fileparts(temp_label_fn);
    [~, image_fn, ~] = fileparts(obj.Nifti_image_xG_filename);
    [~, image_fn, ~] = fileparts(image_fn);

    % Move files to registration folder
    if ~isfile(fullfile(obj.registration_folder,[temp_anat_fn '.nii.gz']))
        copyfile(obj.template_anatomical_filename, obj.registration_folder,'f');
    end
    if ~isfile(fullfile(obj.registration_folder,[temp_bm_fn '.nii.gz']))
        copyfile(obj.template_brain_mask_filename, obj.registration_folder,'f');
    end
    if ~isfile(fullfile(obj.registration_folder,[image_fn '.nii.gz']))
        copyfile(obj.Nifti_image_xG_filename,obj.registration_folder,'f');
    end
    if ~isfile(fullfile(obj.registration_folder,[temp_label_fn '.nii.gz']))
        copyfile(obj.template_labels_filename, obj.registration_folder,'f');
    end

    copyfile(obj.registration_folder,fullfile(obj.tmp_registration_path,'tmp_registration'),'f');
  
    d.Message = 'Processing the registration ...';
    % Do ANTs registration with current image as fixed and template as
    % moving using Docker with antsx/ants installed on it

    toolbox_folder = pwd;
    cd(obj.tmp_registration_path)
    start_docker = '';
    if strcmp(computer,'MACI64')
        start_docker = '/usr/local/bin/docker run -v "./tmp_registration:/data" --rm antsx/ants';
    elseif strcmp(computer,'PCWIN64')
        start_docker = 'docker run -v ".\tmp_registration:/data" antsx/ants';
    end
    [~,cmdout] = system([start_docker ' /opt/ants/bin/antsRegistrationSyNQuick.sh -d 3' ...
        ' -f /data/' image_fn '.nii.gz' ...
        ' -m /data/' temp_anat_fn '.nii.gz' ...
        ' -o ' 'reg_' ...
        ' -t s'],"-echo");

    time_elapsed = 0;
    check_files_reg1 = {fullfile(obj.tmp_registration_path,'tmp_registration','reg_Warped.nii.gz'),...
        fullfile(obj.tmp_registration_path,'tmp_registration','reg_InverseWarped.nii.gz')};
    
    while any(~isfile(check_files_reg1))
        pause(0.5)
        if d.CancelRequested
            system(['kill' cmdout]);
            msg = {'Unfinished business due to cancellation'};
            rmdir(fullfile(obj.tmp_registration_path,'tmp_registration'),'s')
            return
        end
        time_elapsed = time_elapsed + 0.5;
    end
    cd(toolbox_folder)
    
    % Warped image
    movefile(fullfile(obj.tmp_registration_path,'tmp_registration','reg_Warped.nii.gz'), ...
        fullfile(obj.tmp_registration_path,'tmp_registration',[image_fn '_r.nii.gz']),'f');
    obj.registered_image_fn = fullfile(obj.tmp_registration_path,'tmp_registration',[image_fn '_r.nii.gz']);
    
    % Affine transformation
    movefile(fullfile(obj.tmp_registration_path,'tmp_registration','reg_0GenericAffine.mat'), ...
        fullfile(obj.tmp_registration_path,'tmp_registration','affine_transformation.mat'),'f');

    % Space transformation
    if isfile(fullfile(obj.tmp_registration_path,'tmp_registration','reg_1Warp.nii.gz'))
        movefile(fullfile(obj.tmp_registration_path,'tmp_registration','reg_1Warp.nii.gz'), ...
            fullfile(obj.tmp_registration_path,'tmp_registration','space_transformation.nii.gz'),'f');
    end
    
    % Apply the transformation to remove the skull and fat from brain and
    % create the brain mask
    toolbox_folder = pwd;
    cd(obj.tmp_registration_path)
    if isfile(fullfile(obj.Bruker_exp_folder,'space_transformation.nii.gz'))
        [~,~] = system([start_docker ' /opt/ants/bin/antsApplyTransforms -d 3' ...
            ' -i /data/' temp_bm_fn '.nii.gz' ...
            ' -r /data/' image_fn '_r.nii.gz' ...
            ' -o brain_mask.nii.gz' ...
            ' -t /data/space_transformation.nii.gz' ...
            ' -t /data/affine_transformation.mat' ...
            ' -n NearestNeighbor'],"-echo");
    else
       [~,~] = system([start_docker ' /opt/ants/bin/antsApplyTransforms -d 3' ...
            ' -i /data/' temp_bm_fn '.nii.gz' ...
            ' -r /data/' image_fn '_r.nii.gz' ...
            ' -o brain_mask.nii.gz' ...
            ' -t /data/affine_transformation.mat' ...
            ' -n NearestNeighbor'],"-echo");
    end
    cd(toolbox_folder)
    
    % Brain mask
    Brain_mask = logical(niftiread(fullfile(obj.tmp_registration_path,'tmp_registration','brain_mask.nii.gz')));
    obj.Brain_mask = Brain_mask;
    
    % Skull stripping
    obj.skstr_image = double(obj.original_image) .* obj.Brain_mask;

    if ~isfolder(fullfile(obj.result_folder,'Registration'))
        mkdir(fullfile(obj.result_folder,'Registration'))
    end
    save(fullfile(obj.result_folder,'Registration','brain_mask_full.mat'),'Brain_mask')

catch ME
    msg = {ME.message};
    rmdir(fullfile(obj.tmp_registration_path,'tmp_registration'),'s')
end

% Delete other files
delfiles = {fullfile(obj.tmp_registration_path,'tmp_registration', 'reg_InverseWarped.nii.gz'), ...
    fullfile(obj.tmp_registration_path,'tmp_registration', 'reg_Warped.nii.gz'), ...
    fullfile(obj.tmp_registration_path,'tmp_registration', 'reg_1InverseWarp.nii.gz'), ...
    fullfile(obj.tmp_registration_path,'tmp_registration', 'reg_1Warp.nii.gz'), ...
    fullfile(obj.tmp_registration_path,'tmp_registration', 'reg_0GenericAffine.mat')};

for i = 1:length(delfiles)
    if isfile(delfiles{i})
        delete(delfiles{i})
    end
end

end