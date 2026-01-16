% Copyright All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, MRS4Brain research group @ CIBM MRI EPFL AIT, 2024
% See the LICENSE.TXT file for more details.

function  Nbasis = FindNBasisAllCoil_BA(~,Lipids_tkk,mrsiReconParams)
% mrsiReconParams.mrsiData dims: time-k-k

lipid_mask = mrsiReconParams.SkMask2D;
MaxNbasis = 64; %Must be a power of 2

[~,high_bnd_L] = min(abs(mrsiReconParams.LipidMaxPPM - mrsiReconParams.ppm));
[~,low_bnd_L] = min(abs(mrsiReconParams.LipidMinPPM - mrsiReconParams.ppm));

[NCoil,Nt,Nx,Ny] = size(mrsiReconParams.mrsiData);

% Lipid_SingSpect_fsc = zeros(Nt,MaxNbasis,NCoil);
% LipidMask_rrf = zeros(Nx,Ny,Nt);

Lipid_SingSpect_fsc = [];
for coil = 1:NCoil
    CoilLipids_rrf = permute(fft(mrsiReconParams.SENSE(coil,:,:).*ifft(ifft(Lipids_tkk,[],2),[],3),[],1),[2,3,1]);

    Lipid_Vol = squeeze(sum(abs(CoilLipids_rrf),3));
    
    Signal = Lipid_Vol(:);
    if mrsiReconParams.Threshold_LipMask >= 0
        LThr = mean(Signal) + mrsiReconParams.Threshold_LipMask*std(Signal);
    else
        LThr = 0;
    end
    
    CoilLipid_mask=((Lipid_Vol.*lipid_mask) > LThr);
    LipidMask_rrf = repmat(squeeze(CoilLipid_mask),[1 1 Nt]);
    Lipid_stack_rf = reshape(CoilLipids_rrf(LipidMask_rrf > 0),[],Nt);
    [~,~,Vorig] = svd(Lipid_stack_rf,'econ');

    while (nnz(CoilLipid_mask)<MaxNbasis)
        MaxNbasis = MaxNbasis/2;
    end

    Lipid_SingSpect_fsc(:,:,coil) = Vorig(:,1:MaxNbasis);
end

Nbasis = MaxNbasis/2;
Nstep = MaxNbasis/4;
IOP = diag(ones(Nt,1));
% LipidProj = 0*IOP;
% LipFree_rrf = zeros(Nx,Ny,Nt);
CombLipFreeData_rrf = zeros(Nx,Ny,Nt);
lipid_mask = mrsiReconParams.SkMask2D;
BrainMask = mrsiReconParams.BrainMask2D .* ~lipid_mask;

CoilData_crrf = permute(mrsiReconParams.mrsiData,[1,3,4,2]);%ckkt
CoilData_crrf = fft(ifft(ifft(CoilData_crrf,[],2),[],3),[],4);%crrf

LipSupStep = 1;
% figs = figure('Visible', 'off');
while Nstep >= 1
    CombLipFreeData_rrf = CombLipFreeData_rrf*0;
%     fprintf("coil = ");
    for coil = 1:NCoil       
        LipidProj = Lipid_SingSpect_fsc(:,1:Nbasis,coil)*Lipid_SingSpect_fsc(:,1:Nbasis,coil)';
    
        LipFree_rrf = reshape(reshape(CoilData_crrf(coil,:,:,:),[],Nt)*(IOP-LipidProj),[Nx Ny Nt]); 
        
        CombLipFreeData_rrf = CombLipFreeData_rrf + squeeze(conj(mrsiReconParams.SENSE(coil,:,:))).*LipFree_rrf;       
    end
   
    EnergyMap = sum(abs(CombLipFreeData_rrf(:,:,low_bnd_L:high_bnd_L).^2),3);
    SkullEnergy = quantile(EnergyMap(lipid_mask > 0),0.97);
    BrainEnergy = quantile(EnergyMap(BrainMask > 0),0.97);
    Ratio = (BrainEnergy/SkullEnergy);

    if(mrsiReconParams.L2SVDparams.PercentThres) > Ratio
		Nbasis = Nbasis + Nstep;
    else
		Nbasis = Nbasis - Nstep;
    end
	Nstep = Nstep/2;

    LipSupStep = LipSupStep + 1;
end

% close(figs)

if Nbasis > mrsiReconParams.L2SVDparams.NBasisMax
    Nbasis = mrsiReconParams.L2SVDparams.NBasisMax;
elseif Nbasis < mrsiReconParams.L2SVDparams.NBasisMin
    Nbasis = mrsiReconParams.L2SVDparams.NBasisMin;
end

end