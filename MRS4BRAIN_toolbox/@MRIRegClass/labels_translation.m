%%  labels_translation.m 

% Copyright All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, MRS4Brain research group @ CIBM MRI EPFL AIT, 2024
% See the LICENSE.TXT file for more details.

% Guillaume Briand, CIBM - MRS4Brain group, 2023
% 
% USAGE : MRIReg Class public method 
% msg = obj.labels_translation()
% 
% DESCRIPTION :
% Translate the labelled MRI image with the regions names and calculate the
% volume of each brain region
%
% INPUTS :
% obj       = MRIReg Class object with properties and methods
%
% OUTPUT :
% msg       = Error message
function msg = labels_translation(obj)
% Function to translate the labels and compute the volumetry
msg = {''};
try
    label_file = fullfile(obj.registration_folder, 'brain_labels.nii.gz');
    if isfile(label_file)
        % Volume of a single voxel
        info_label = niftiinfo(label_file);
        v = info_label.PixelDimensions / 10;
        obj.voxel_volume = v(1) * v(2) * v(3);
        V_label = niftiread(info_label);

        % Load Labels names and numbers
        mlabel_temp = load(obj.template_labels_table);
        fns_mlabel_temp = fieldnames(mlabel_temp);
        mlabel = mlabel_temp.(fns_mlabel_temp{1});
        segmentation_oldfieldnames = fieldnames(mlabel);
        segmentation_newfieldnames = {'Name','Left','Right'};
        for k = 1:numel(segmentation_newfieldnames)
            if ~strcmp(segmentation_newfieldnames{k},segmentation_oldfieldnames{k})
                [mlabel.(segmentation_newfieldnames{k})] = mlabel.(segmentation_oldfieldnames{k});
                mlabel = orderfields(mlabel,[1:(k-1),4,k:3]);
                mlabel = rmfield(mlabel,segmentation_oldfieldnames{k});
            end
        end
        for j = 1:length(mlabel)
            mlabel(j).Name = convertCharsToStrings(mlabel(j).Name);
        end

        volumetry = mlabel;
        
        % Update the volumetry and mlabel struct
        n_regions = length(mlabel);
        for r = 1:n_regions
            V_left = zeros(size(V_label));
            V_right = zeros(size(V_label));
            V_left(V_label == mlabel(r).Left) = 1;
            V_right(V_label == mlabel(r).Right) = 1;
            mlabel(r).Left = rot90(V_left,-1);
            mlabel(r).Right = rot90(V_right,-1);
            volumetry(r).Left = squeeze(sum(sum(V_left,1),2)) * obj.voxel_volume;
            volumetry(r).Right = squeeze(sum(sum(V_right,1),2)) * obj.voxel_volume;
            volumetry(r).LR = volumetry(r).Left + volumetry(r).Right;
        end

        obj.Segmentation = mlabel;
        obj.Volumetry = volumetry;

    else
        msg = {['No exisiting segmentation, ' ...
            'please procede the image segmentation before translating the labels !']};
        return
    end

catch ME
    msg = {ME.message};
end
end