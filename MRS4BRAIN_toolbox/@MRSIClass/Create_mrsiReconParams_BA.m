% Copyright All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, MRS4Brain research group @ CIBM MRI EPFL AIT, 2024
% See the LICENSE.TXT file for more details.

function mrsiReconParams = Create_mrsiReconParams_BA(obj,mrsiData_tkk,PercentThres,BrainMap)
mrsiReconParams.results_folder = obj.results_folder;

MatSize = obj.acq_params.matrix_sz;

mrsiReconParams.Log_Dir = obj.data_folder;
mrsiReconParams.InitialB0MapFromWater = 1; % yes / no, Compute a initial B0 map from water signal (generally 1 but might be incorrect in phantoms)

mrsiReconParams.AcqDelay = obj.acq_params.acq_delay; % 0 for a Spin-echo. For FID-MRSI sequence, put the acquisition delay in second. First missing point of the FID will be predicted to correct the 1st order phase

mrsiReconParams.MinPPM = -5.3; % begining of the frequency window for low-rank Recon. Put to -1000 for full window
mrsiReconParams.MaxPPM = 1000.0; % End of the frequency window for low-rank Recon.Put to 1000 for full window

mrsiReconParams.DoHomogeneityCorrection = 0; % yes / no % use the Water signal to uniform the metabolite signal during reconstruction. This factor is also applied to the resulting Water signal so the final quantification is not affected by this option  .DoHomogeneityCorrection leads to instabilities with high B1 inhomogeneity and enhance lowest slice signal by strong factor. Not recommended
mrsiReconParams.HammingFidFilter = 0;% Filter applied on the recon data fidelity terms. 0 - no filter , 1 max filtering. It helps the convergence in case of very low SNR Data

mrsiReconParams.Threshold_LipMask = 0.5; % Threshold to reduce the size of the skull mask for each coil indivudually (Used for lipid suppression only)

mrsiReconParams.UndersamplingF = 1; % Retrospective undesampling (for acceleration simulation purpose)

mrsiReconParams.ZeroPaddingF = 0;% Pad the end of the data: 1 add 100%, 0.5 add 50%, 0 no padding .(LCModel works somehow better with padded data)

mrsiReconParams.L2SVDparams.PercentThres = PercentThres; % strength of the Lipid suppresion
mrsiReconParams.L2SVDparams.Nfit = 20; %
mrsiReconParams.L2SVDparams.NBasisMax = 128; % Max allowed components for Lipid suppression
mrsiReconParams.L2SVDparams.NBasisMin = 5; % Min allowed components for Lipid suppression

mrsiReconParams.LipidMinPPM = -3; % Range of lipid in ppm (used to compute lipid contamination in lipid suppression)
mrsiReconParams.LipidMaxPPM = 0.0; % Range of lipid in ppm (used to compute lipid contamination in lipid suppression)
mrsiReconParams.NbPtForWaterPhAmp = 5; % First pt of the time series taken for Amplitude and Phase calculation
mrsiReconParams.ESPIRIT_kernel = [6,6]; % Size of the ESPIRIT kernel for coil sensitivity profile computation


mrsiReconParams.FiltParam.Water_minFreq = -obj.acq_params.resfreq*2; % hz %Min freq for water removal by HSVD (0Hz = 4.7ppm)
mrsiReconParams.FiltParam.Water_maxFreq = obj.acq_params.resfreq/2; % hz %Max freq for water removal by HSVD (0Hz = 4.7ppm)
mrsiReconParams.FiltParam.Comp = 16; % Nb of component for the HSVD water removal (advised: 16 at 3T and 32 at 7T)

mrsiReconParams.GaussianSigma = 3; % The kernel width for B0 fieldmap correction

mrsiReconParams.LRTGVModelParams.Maxit = 1500; % Maximum number of iteration
mrsiReconParams.LRTGVModelParams.Minit = 600; % Minimum number of iteration
mrsiReconParams.LRTGVModelParams.check_it = 25; % Step where the iterative recon check the different convergence parameters (output in log file)
mrsiReconParams.LRTGVModelParams.Plot_it = 100; % Step where the iterative recon make plot of the temporary results (output in Log_Files.../LowRankTGV_Recon/)
mrsiReconParams.LRTGVModelParams.CorrB0Map_it = 25; % Step where the iterative recon estimatea a correction to the B0 Field Map
mrsiReconParams.LRTGVModelParams.CorrB0Map_Maxcount = 100; % Max amount of  B0 Field Map correction
mrsiReconParams.LRTGVModelParams.CorrB0Map_MaxPPM = -0.5;
mrsiReconParams.LRTGVModelParams.CorrB0Map_MinPPM = -3.5;
mrsiReconParams.LRTGVModelParams.Orthogonalize_it = 25; % Step where the iterative recon forces orthogonality in the temporaal & spatial component
mrsiReconParams.LRTGVModelParams.SpecItFact = 2;
mrsiReconParams.LRTGVModelParams.reduction = 100; %100 invivo %1000 phantom % Nb of iterations for soft introduction of the regularization in the recon
mrsiReconParams.LRTGVModelParams.min_SpectStep = 1/32; % minimum spectral update step
mrsiReconParams.LRTGVModelParams.max_SpectStep = 1/2; %1/2 invivo, 1/4 if diverge 1/8 phantom % maximum spectral update step
mrsiReconParams.LRTGVModelParams.min_taup = 1/32; % minimal spatial update step
mrsiReconParams.LRTGVModelParams.max_taup = 1/8; % %1/8 invivo , 1/16 if diverge, 1/16 for Synthetic Data, 1/32 phantom  % maximal spatial update step

mrsiReconParams.ppm = (-obj.acq_params.ppm_ref + ((1:obj.acq_params.np_met)* ...
    obj.acq_params.spectralwidth/(obj.acq_params.np_met * obj.acq_params.resfreq)));

BrainMask = ones(MatSize);
ImMask = ones(MatSize);
SkMask = 0; % SkullMask

mrsiReconParams.BrainMask        = BrainMask;
mrsiReconParams.ImMask           = ImMask;
mrsiReconParams.SkMask           = SkMask;
mrsiReconParams.mrsiData(1,:,:,:) = mrsiData_tkk;
mrsiReconParams.SENSE = ones([1 MatSize]);

mrProt.samplerate = obj.acq_params.spectralwidth;
mrProt.DwellTime = 1.0/mrProt.samplerate;
mrsiReconParams.mrProt = mrProt;
mrsiReconParams.mrProt.VSize = obj.acq_params.np_met;
mrsiReconParams.mrProt.NMRFreq = obj.acq_params.resfreq;

Size_data = [1,obj.acq_params.np_met,MatSize(1),MatSize(2)];

mrsiReconParams.BrainMask2D = BrainMap;
mrsiReconParams.BrainMask2D = round(mrsiReconParams.BrainMask2D/max(mrsiReconParams.BrainMask2D(:)));

mrsiReconParams.ImMask2D = imresize(mean(mrsiReconParams.ImMask,3),[Size_data(3),Size_data(4)]);
mrsiReconParams.ImMask2D = round(mrsiReconParams.ImMask2D/max(mrsiReconParams.ImMask2D(:)));

mrsiReconParams.SkMask2D = mrsiReconParams.ImMask2D - mrsiReconParams.BrainMask2D;
mrsiReconParams.SkMask2D(mrsiReconParams.SkMask2D<0) = 0;
end