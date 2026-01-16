% Copyright All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, MRS4Brain research group @ CIBM MRI EPFL AIT, 2024
% See the LICENSE.TXT file for more details.

classdef MRIRegClass < handle
    % RegMRIClass is a class object to register rat brain MRI images with an
    % anatomical template and apply segmentation of brain regions
    % This Class requires an external software ANTs that can be installed
    % with the document that can been found in the MRS4Brain Toolbox folder
    properties(Access = public)
        % MRI registration parameters
        MRIReg_parameters               % Registration MRI parameters

        template_anatomical_filename    % Anatomical rat brain template
        template_brain_mask_filename    % Anatomical rat brain mask
        template_labels_filename        % Rat brain labelled template
        template_labels_table           % Rat brain template labels table
        
        % Initial image filenames
        Nifti_image_filename            % Rat brain MRI image 
        Nifti_image_xG_filename         % Rat brain MRI image with updated pixel dimension
        
        % Final Data
        Bruker_exp_folder               % Bruker experiment folder directory
        result_folder                   % MRSI/MRI registration Result folder directory
        registration_folder             % Folder with registration materials
        tmp_registration_path           % Folder with temporary registration materials

        % Images
        original_image                  % MRI rat brain original image
        registered_image_fn             % MRI registered image to template 
        skstr_image                     % Skull stripped original image

        % Labels
        voxel_volume                    % Volume of an MRI voxel
        Brain_mask                      % Mask of the rat brain
        Segmentation                    % Segmentation of the rat brain

        % Volumetry
        Volumetry                       % Volumetry of each brain regions (Left and Right)
    end

    methods(Access = public)
        function obj = MRIRegClass(Nifti_image,mrireg_param,Bruker_folder,result_folder,head_prone)
            
            % MRI REGISTRATION PARAMETERS 
            obj.MRIReg_parameters = mrireg_param;
            
            % INPUT IMAGE
            obj.Nifti_image_filename = Nifti_image;
            
            % TEMPLATE
            obj.template_anatomical_filename = mrireg_param.template_anatomical;
            obj.template_brain_mask_filename = mrireg_param.template_brain_mask;
            obj.template_labels_filename = mrireg_param.template_seg;
            obj.template_labels_table = mrireg_param.labels_seg;
            
            % OUTPUT FOLDERS
            obj.Bruker_exp_folder = Bruker_folder;
            obj.registration_folder = fullfile(Bruker_folder,'registration');
            obj.result_folder = result_folder;
            
            % CHANGE PX DIMENSION AND READ INPUT IMAGE
            obj.update_px_dim(10,head_prone);
            obj.original_image = double(niftiread(obj.Nifti_image_xG_filename));
        end
        
        varargout = update_px_dim(varargin);            % Change the dimension of pixels to correspond to human brain
        varargout = image_registration(varargin);       % Register the MRI rat brain to the anatomical template
        varargout = brain_mask(varargin);               % Create a brain mask from a segmentation Nifti file
        varargout = image_segmentation(varargin);       % Apply brain segmentation using the transformations found
        varargout = labels_translation(varargin);       % Translate the labels number to corresponding regions
    end
end