% Copyright All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, MRS4Brain research group @ CIBM MRI EPFL AIT, 2024
% See the LICENSE.TXT file for more details.

function [mrsiDataLR_tkk,LipidFree_frr,Lipid_rrf,NBasis,...
    LipidProj] = ProjSVDLipidSuppression_BA(~,mrsiData_tkk,LipData_tkk,mrsiReconParams,NBasis,NameData)

 % mrsiReconParams.mrsiData dims: time-k-k
Data_kkf = fft(permute(mrsiData_tkk,[2,3,1]),[],3);
Data_rrf = ifft(ifft(Data_kkf,[],1),[],2);

LipData_rrt = ifft(ifft(permute(LipData_tkk,[2,3,1]),[],1),[],2);

N = size(LipData_rrt);
% HzpP = mrsiReconParams.mrProt.samplerate/N(3);
 
LipData_rrf = fft(LipData_rrt,[],3);

lipid_mask = mrsiReconParams.SkMask2D;

[~,high_bnd_L] = min(abs(mrsiReconParams.LipidMaxPPM - mrsiReconParams.ppm));
[~,low_bnd_L] = min(abs(mrsiReconParams.LipidMinPPM  - mrsiReconParams.ppm));

Lipid_Vol = squeeze(sum(abs(LipData_rrf(:,:,low_bnd_L:high_bnd_L)),3));
Lipid_Vol = Lipid_Vol/max(Lipid_Vol(:));
SortedLip = sort(Lipid_Vol(:),'descend');

Noise = Lipid_Vol(:);%=Lipid_Vol(mrsiReconParams.SkMask2D==1);
if mrsiReconParams.Threshold_LipMask >= 0
	LThr = mean(Noise) + mrsiReconParams.Threshold_LipMask*std(Noise);%+%1.0*std(Noise); %quantile(Lipid_Vol(:),0.5);%0; 
else
	LThr = 0;
end

% [~,NumVox] = min(abs(LThr-SortedLip(:)));
% fprintf(["NumVox="+num2str(NumVox)]);
lipid_mask = ((Lipid_Vol.*lipid_mask) > LThr);

% meta_mask = mrsiReconParams.BrainMask2D;
% Brainmask_rr1 = reshape(meta_mask,[size(meta_mask,1) size(meta_mask,2) 1]);
%% apply L2 lipid-basis recon to dual-density image

SVRatio = mrsiReconParams.L2SVDparams.PercentThres;
if (nargin == 4) % no Nbasis given
    [LipidProj,~,NBasis,~] = make_SVD_LipidBasis(LipData_rrf,Data_rrf,lipid_mask,SVRatio,mrsiReconParams); 
elseif (nargin > 4) % Nbasis given
    [LipidProj,~,NBasis,~] = make_SVD_LipidBasis(LipData_rrf,Data_rrf,lipid_mask,SVRatio,mrsiReconParams,NBasis,NameData);    
end

Lipid_kkf = reshape(Data_kkf,[],N(end)) * LipidProj;
Lipid_kkf = reshape(Lipid_kkf,[N(1),N(2),N(3)]);

Lipid_kkt = ifft(Lipid_kkf,[],3);
Lipid_rrf = ifft(ifft(Lipid_kkf,[],1),[],2);
mrsiDataLR_tkk = (mrsiData_tkk - permute(Lipid_kkt,[3,1,2]));
LipidFree_frr = fft(ifft(ifft(mrsiDataLR_tkk,[],2),[],3),[],1);

% [~,high_bnd_L] = min(abs(0 - mrsiReconParams.ppm));
% [~,low_bnd_L] = min(abs(-4.7 - mrsiReconParams.ppm));

end

function [LipidProj,Slipid,Nbasis,Lipids] = make_SVD_LipidBasis(lipid_rrf,data_rrf,lipid_mask,SVratio,ReconParams,Nbasis,NameData)
%GET_LIPIDBASIS Summary of this function goes here
%   Detailed explanation goes here

count = 0;
N = size(data_rrf);
IOP = diag(ones(N(end),1));
Lipid_Basis = zeros(N(3), sum(lipid_mask(:)));

[~,high_bnd_L] = min(abs(ReconParams.LipidMaxPPM - ReconParams.ppm));
[~,low_bnd_L] = min(abs(ReconParams.LipidMinPPM  - ReconParams.ppm));

% SVD basis of the skull lipids 

LipidMask_rrf = repmat(lipid_mask,[1 1 N(end)]);
[~,Sorig,Vorig] = svd(reshape(lipid_rrf(LipidMask_rrf > 0),[],N(end)),'econ');

if (nargin <6) % No Nbasis given
	Nbasis = 64;%32
	Nstep = 32;%16
	while Nstep >= 1	    
	    LipidProj = Vorig(:,1:Nbasis) * Vorig(:,1:Nbasis)';
	    
	    LipFree_rrf = reshape(reshape(data_rrf,[],N(end))*(IOP-LipidProj),[N(1) N(2) N(3)]);
	    EnergyMap=sum(abs(LipFree_rrf(:,:,low_bnd_L:high_bnd_L)).^2,3);
	    
	    SkullEnergy = mean(EnergyMap(lipid_mask > 0));
	    BrainEnergy = mean(EnergyMap(ReconParams.BrainMask2D > 0));
	    %BrainEnergy/SkullEnergy;
	    if SkullEnergy > (BrainEnergy/SVratio)
		    Nbasis = Nbasis + Nstep;
	    else
		    Nbasis = Nbasis - Nstep;
	    end
	    Nstep = Nstep/2;
	end

    if Nbasis > ReconParams.L2SVDparams.NBasisMax
        Nbasis = ReconParams.L2SVDparams.NBasisMax;
    elseif Nbasis < ReconParams.L2SVDparams.NBasisMin
        Nbasis = ReconParams.L2SVDparams.NBasisMin;
    end

else % Nbasis given
	LipidProj = Vorig(:,1:Nbasis) * Vorig(:,1:Nbasis)';
	    
	LipFree_rrf = reshape(reshape(data_rrf,[],N(end))*(IOP-LipidProj),[N(1) N(2) N(3)]);
	EnergyMap = sum(abs(LipFree_rrf(:,:,low_bnd_L:high_bnd_L)).^2,3);
	    
	SkullEnergy = mean(EnergyMap(lipid_mask > 0));
	BrainEnergy = mean(EnergyMap(ReconParams.BrainMask2D > 0));

end	
Lipids = Vorig(:,1:Nbasis);
Slipid = Sorig(1:Nbasis,1:Nbasis);
LipidProj = Lipids*Lipids';

if (nargin==7)
    s=[ReconParams.results_folder,'\LipidMaps\Lipid_Images.ps'];
    if exist(s); delete(s); end
    figs = figure('visible', 'off');
    subplot(211); plot(ReconParams.ppm,real(Lipids)); 
    legend(strcat(repmat('#',size(Lipids,2),1),num2str((1:size(Lipids,2))')),'Box','off','Location','eastoutside');
    grid on; xlabel 'ppm'; xlim tight; title 'Lipid spectra in full range';
    set(gca,'Box','off','XDir','reverse','TickDir','out','Color','none');
    subplot(212); plot(ReconParams.ppm,real(Lipids)); 
    legend(strcat(repmat('#',size(Lipids,2),1),num2str((1:size(Lipids,2))')),'Box','off','Location','eastoutside');
    grid on; xlabel 'ppm'; xlim([0 4]); title('Lipid spectra in 0 - 4 ppm');
    set(gca,'Box','off','XDir','reverse','TickDir','out','Color','none');
    print(figs,'-append','-dpsc2',s); warning('off');
    
    figs = figure('visible','off');
    imagesc(~lipid_mask.*sum(abs(lipid_rrf),3));%,[ 0,10*mean(image2plot(:))] )
    colormap default; colorbar;
    title('Filtered out Lipids in Brain & Outside head')
    print(figs,'-append','-dpsc2',s); warning('off');

    imagesc(~lipid_mask.*sum(abs(data_rrf),3));%,[ 0,10*mean(image2plot(:))] )
    colormap default; colorbar;
    title('Original Data in Brain & Outside head')
    print(figs,'-append','-dpsc2',s); warning('off');

    imagesc(~lipid_mask.*squeeze(sum(abs(LipFree_rrf),3))); %,[ 0,10*mean(image2plot(:))] )
    title('Lipid-Free Data in Brain & Outside head');
    colormap default; colorbar;
    print(figs,'-append','-dpsc2',s); warning('off');

    imagesc(sum(abs(lipid_rrf),3));%,[ 0,10*mean(image2plot(:))] )
    colormap default; colorbar;
    title('Filtered out Lipids in Image')
    print(figs,'-append','-dpsc2',s); warning('off');


    imagesc(sum(abs(data_rrf),3));%,[ 0,10*mean(image2plot(:))] )
    colormap default; colorbar;
    title('Original Data in Image')
    print(figs,'-append','-dpsc2',s); warning('off');

    imagesc(squeeze(sum(abs(LipFree_rrf),3)));%,[ 0,10*mean(image2plot(:))] )
    title('Lipid-Free Data in Image');
    colormap default; colorbar;
    print(figs,'-append','-dpsc2',s); warning('off');

    imagesc(lipid_mask);%,[ 0,10*mean(image2plot(:))] )
    title('Lipid masks');
    colormap default; colorbar;
    print(figs,'-append','-dpsc2',s); warning('off');
end
end
