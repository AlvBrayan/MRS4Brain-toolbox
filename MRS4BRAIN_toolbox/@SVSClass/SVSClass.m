% Copyright All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, MRS4Brain research group @ CIBM MRI EPFL AIT, 2024
% See the LICENSE.TXT file for more details.

classdef SVSClass < handle
    % SVSClass is a class object to process Single Voxel Spectroscopy
    % experiments from Bruker format, fit and quantify the spectra by using
    % LCModel
    properties(Access = public)
        % SVS Properties
        SVS_param                       % SVS parameters for preprocessing and quantification
        
        % Result folder
        result_dir                      % Directory of result folder
        foldername                      % Name of the result folder

        % Data 
        SVS_struct                      % Structure with each experiment, their type and with study structures

        % Metabolite concentrations
        Met_names                       % Names of quantified metabolites
        abs_conc_tot                    % Absolute metabolite concentration for all metabolite fid
        rel_conc_tot                    % Relative metabolite concentration for all metabolite fid
        crlb_tot                        % CRLB of each metabolite for all metabolite fid (Cramer-Rao Lower Bound)

        % LCModel data fitting
        ppm_scale                       % ppm scale from LCModel quantification
        data_pts_tot                    % Data points 
        fit_pts_tot                     % LCModel fit on data points
        baseline_pts_tot                % Baseline calculated by LCModel on data points with BASIS set used
    end


    methods(Access = public)
        
        % Class constructor
        function obj = SVSClass(result_dir,foldername,SVS_struct,SVS_parameters)
            obj.result_dir = result_dir;
            obj.foldername = foldername;
            obj.SVS_struct = SVS_struct;
            obj.SVS_param = SVS_parameters;
        end

        varargout = create_Bruker_rawstudy(varargin);           % Read the Brucker fid data files and convert them 
        varargout = preprocessing_SVS_FidA(varargin);           % Apply FidA preprocessing on fid data
        varargout = SVS_quantification(varargin);               % Fitting and quantification of spectra with LCModel
        varargout = plot_signal(varargin);                      % Return ppm scale and spectrum

    end

    methods(Access = private)
        
        varargout = convert2FidA_isison(varargin);              % Convert study structure format to FidA format with ISIS
        varargout = convert2FidA_isisoff(varargin);             % Convert study structure format to FidA format without ISIS
        varargout = preprocessing_SVS_FidA_isison(varargin);    % Apply FidA preprocessing on fid data with ISIS
        varargout = preprocessing_SVS_FidA_isisoff(varargin);   % Apply FidA preprocessing on fid data without ISIS
        varargout = alignAverages_FidA(varargin);               % FidA function to align Averages
        varargout = rmbadaverages_FidA(varargin);               % FidA function to remove badAverages
        varargout = read_table(varargin);                       % Read table file with metabolite concentrations
        varargout = read_coord(varargin);                       % Read coord file with LCModel fitting, baseline
        varargout = data2RAW(varargin);                         % Save study struct data as RAW file for LCModel quantif

    end

end