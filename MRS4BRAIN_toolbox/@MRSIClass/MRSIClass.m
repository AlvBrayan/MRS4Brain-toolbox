% Copyright All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, MRS4Brain research group @ CIBM MRI EPFL AIT, 2024
% See the LICENSE.TXT file for more details.

classdef MRSIClass < handle
    % MRSClass is a class with the properties and methods to process and
    % quantify MRSI data from a Brucker experiment 
    % The 
    properties(Access = public)
        % MRSI Properties
        mrsi_params                     % MRSI parameters (LCModel)
        acq_params                      % Acquistion parameters
        Nslices                         % Number of MRSI slices
        Slice_range                     % Slice range from the central slice number
        Slices_number                   % MRI central slices
        Lipsup = false                  % Lipid suppression enabled / disabled
        Fillgaps = false                % Fillgaps enabled / disabled

        % Initial Data
        data_folder                     % Directory of the Brucker data
        metab_expnb                     % Folder with the metabolite fid
        ref_expnb                       % Folder with the reference fid
        
        % Intermediate Data
        fid_mat_tkkn                    % Fid of the metabolites in Fourier domain
        ref_mat_tkkn                    % Fid of the reference in Fourier domain
        Brain_mask                      % Brain mask computed from a MRI segmentation
        Power_map                       % Powermap based on Linewidth
        Linewidth_map                   % Water Linewidth map
        SNR_map                         % SNR map
        HSVD_fid_tkkn                   % Applied custom_HSVD on metabolites fid_mat_tkkn
        HSVD_lipsup_fid_tkkn            % Applied Lipsup_MRSI on metabolites HSVD_fid_tkkn
        HSVD_lipsup_filled_fid_tkkn     % Applied Fillgaps_MRSI on metabolites HSVD__lipsup_fid_tkkn
        
        % Final data
        results_folder                  % Directory of the result folder
        Final_mask                      % Final mask used to quantify the data
        Final_met_map                   % Metabolite map based on the quantified spectra
        Met_names                       % Name of the metabolites in the LCModel basis set
    end


    methods(Access = public)
        
        % MRSIClass class constructor
        function obj = MRSIClass(data_folder, metab_expnb, ref_expnb, ...
                mrsi_params, Slice_range, Slices_number)
            obj.data_folder = data_folder;
            obj.metab_expnb = metab_expnb;
            obj.ref_expnb = ref_expnb;
            obj.Slices_number = Slices_number;
            obj.mrsi_params = mrsi_params;
            obj.Slice_range = Slice_range;
            obj.Nslices = length(Slices_number);
        end

        varargout = check_prexisting_results(varargin);     % Check if the data have already been quantified or already saved as Raw files
        varargout = read_data(varargin);                    % Read the Brucker fid data files and convert them 
        varargout = Brain_mask_comp(varargin);              % Compute the Brain mask from a MRI Brain mask with higher resolution
        varargout = Linewidth_map_comp(varargin);           % Compute the Linewidth map for the reference spectra
        varargout = custom_HSVD(varargin);                  % Apply the HSVD method on the metabolite fid data
        varargout = SNR_map_comp(varargin);                 % Compute the SNR map
        varargout = quality_masks(varargin);                % Compute the quality masks and apply it to the Final mask
        varargout = LipSup_MRSI(varargin);                  % Apply the Lipid suppression method on the metabolite data
        varargout = Fillgaps_MRSI(varargin);                % Apply the Fillgaps method on the metabolite data
        varargout = save_data(varargin);                    % Save the data as RAW files for the quantification
        varargout = quantify_data(varargin);                % Quantification with LCModel 
        varargout = read_coord_tables(varargin);            % Read LCModel table/coord files to create metabolite maps
        varargout = plot_signal(varargin);                  % Plot a spectra at specific coordinates

    end

    methods(Access = private)
        
        % Lipid suppression functions
        varargout = Create_mrsiReconParams_BA(varargin);
        varargout = FindNBasisAllCoil_BA(varargin);
        varargout = ProjSVDLipidSuppression_BA(varargin);
        % Save as RAW data
        varargout = data2RAW(varargin);
        % 0 Order Phase correction functions
        varargout = MRSI_0orderphasecorrection(varargin);
        varargout = MRSI_0orderphasecorrection_water(varargin);
        varargout = order0phasecorrection(varargin);
        % Anatomical Map Combiner
        varargout = Anatomical_map_combiner(varargin);

    end


end