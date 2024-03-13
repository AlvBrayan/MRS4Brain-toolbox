%%  preprocessing_SVS_FidA_isison.m 

% Copyright All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, MRS4Brain research group @ CIBM MRI EPFL AIT, 2024
% Copyright 2020 Jamie Near
% See the LICENSE.TXT file for more details.

% Jessie Mosso, CIBM - MRS4Brain group, LIFMET, 2021
% Guillaume Briand, CIBM - MRS4Brain group, 2023
% 
% USAGE : SVS Class public method 
% out = obj.preprocessing_SVS_FidA_isison(prog_dbox)
% 
% DESCRIPTION :
% Apply FidA preprocessing with ISIS on
%
% INPUTS :
% obj       = SVS Class object with properties and methods
% prog_dbox = MRS4Brain Toolbox progress dialog box
%
% OUTPUT :
% msg       = Error message
function msg = preprocessing_SVS_FidA_isison(obj,prog_dbox)
msg = {''};
lcmodelfac = 1; % PARAMETER
svs_param = obj.SVS_param;

for i = 1:length(obj.SVS_struct)
    prog_dbox.Message = ['Apply preprocessing for : ',obj.SVS_struct(i).exp_name];
    prog_dbox.Value = i/(length(obj.SVS_struct)+1);
    try
        study = obj.SVS_struct(i).raw_study;
        if study.multiplicity > 1
            %% step 1 - do preprocessing with FIDA
            %% 1-Convert Bruker study to FID A structure
            [out0a,out0b] = obj.convert2FidA_isison(i); % Convert study struct to FidA struct
            % apply Line Broadening (LB)
            dw = out0a.dwelltime;
            tt = 0:dw:dw*(out0a.n-1);
            % impair
            out0alb = out0a;
            fids0alb = out0alb.fids.*repmat(exp(-tt*pi*svs_param.LBall).',1,out0a.averages);
            out0alb.fids = fids0alb;
            out0alb.specs = fftshift(fft(out0alb.fids.',[],2),2).';
            % pair
            out0blb = out0b;
            fids0blb = out0blb.fids.*repmat(exp(-tt*pi*svs_param.LBall).',1,out0b.averages);
            out0blb.fids = fids0blb;
            out0blb.specs = fftshift(fft(out0blb.fids.',[],2),2).';
    
            %% 2-align av.
            if svs_param.FidA_align_avg
                [out1alb,fsa,phsa] = obj.alignAverages_FidA(out0alb, ...
                    svs_param.minppm,svs_param.maxppm,svs_param.t_max,'y');%5.2,9,0.5,'y');
                [out1blb,fsb,phsb] = obj.alignAverages_FidA(out0blb, ...
                    svs_param.minppm,svs_param.maxppm,svs_param.t_max,'y');%5.2,9,0.5,'y');
                % fs        = Vector of frequency shifts (in Hz) used for alignment.
                % phs       = Vector of phase shifts (in degrees) used for alignment.
    
                % remove Line Broadening (LB)
                out1a = out1alb;
                fids1a = out1alb.fids.*repmat(exp(tt*pi*svs_param.LBall).',1,out1alb.averages);
                out1a.fids = fids1a;
                out1a.specs = fftshift(fft(out1a.fids.',[],2),2).';
                out1b = out1blb;
                fids1b = out1blb.fids.*repmat(exp(tt*pi*svs_param.LBall).',1,out1blb.averages);
                out1b.fids = fids1b;
                out1b.specs = fftshift(fft(out1b.fids.',[],2),2).';
            else
                fsa = zeros(round(study.multiplicity/2),1); fsb = fsa; phsa = fsa; phsb = fsa;
                out1a = out0a;
                out1b = out0b;
            end
    
            %% 3-outlier removal
            if svs_param.FidA_rm_badavg
                [out2a,metrica,badAveragesa] = obj.rmbadaverages_FidA(out1a,svs_param.sd_thresh,'f'); %performs 10Hz LB inside
                [out2b,metricb,badAveragesb] = obj.rmbadaverages_FidA(out1b,svs_param.sd_thresh,'f'); %performs 10Hz LB inside
        
                %apply LB
                out2alb = out2a;
                fids2alb = out2a.fids.*repmat(exp(-tt*pi*svs_param.LBall).',1,out2a.averages);
                out2alb.fids = fids2alb;
                out2alb.specs = fftshift(fft(out2alb.fids.',[],2),2).';
                out2blb = out2b;
                fids2blb = out2b.fids.*repmat(exp(-tt*pi*svs_param.LBall).',1,out2b.averages);
                out2blb.fids = fids2blb;
                out2blb.specs = fftshift(fft(out2blb.fids.',[],2),2).';
            else
                metrica = zeros(round(study.multiplicity/2),1); metricb = metrica;
                badAveragesa = []; badAveragesb = [];
            end

            %combine on/off
            fida = out1a.fids.';
            fidb = out1b.fids.';
            fidtot(1:2:size(fida,1)*2,:) = fida;
            fidtot(2:2:size(fida,1)*2,:) = fidb;

            fid2sum = zeros(size(fidtot,1)/2,size(fidtot,2));
            for k = 1:size(fidtot,1)/2
                fid2sum(k,:) = sum(fidtot((k-1)*2+1:k*2,:));
            end

            ind = 1;
            fidmocor = zeros(size(fid2sum,1),size(fid2sum,2));
            for k = 1:size(fidtot,1)/2
                if ~ismember(k,badAveragesa)
                    if ~ismember(k,badAveragesb)
                        fidmocor(ind,:) = fid2sum(k,:);
                        ind = ind + 1;
                    end
                end
            end
            fidmocor(ind:end,:) = [];
            
            fidmocor = conj(fidmocor);
    
            %% 4-add all the info to the Matlab study structure

            processed_study = study;
            processed_study.fidaprocess.phsa = phsa;
            processed_study.fidaprocess.fsa = fsa;
            processed_study.fidaprocess.metrica = metrica;
            processed_study.fidaprocess.badAveragesa = badAveragesa;
            processed_study.fidaprocess.phsb = phsb;
            processed_study.fidaprocess.fsb = fsb;
            processed_study.fidaprocess.metricb = metricb;
            processed_study.fidaprocess.badAveragesb = badAveragesb;

            processed_study.params.nt = size(fidmocor,1);
            processed_study.multiplicity = size(fidmocor,1);

            processed_study.process.apodparam1 = zeros(1,size(fidmocor,1));
            processed_study.process.apodparam2 = zeros(1,size(fidmocor,1));
            processed_study.process.phasecorr0 = zeros(1,size(fidmocor,1));
            processed_study.process.phasecorr1 = zeros(1,size(fidmocor,1));

            processed_study.data.real = zeros(size(fidmocor,1),1,size(fidmocor,2));
            processed_study.data.real(:,1,:) = real(fidmocor);
            processed_study.data.imag = zeros(size(fidmocor,1),1,size(fidmocor,2));
            processed_study.data.imag(:,1,:) = imag(fidmocor);

            %% STEP 2 - APPLY PREPROCESSING, SUM AND PHASE THE SUM
    
            fidmocor = squeeze(processed_study.data.real) + 1j*squeeze(processed_study.data.imag);
            sumfid = sum(fidmocor); %./size(fidmocor,1);
            sumfid = sumfid .* lcmodelfac;
    
            %% save
            sum_processed_study = processed_study;
            sum_processed_study.data.real = zeros(1,1,study.np/2);
            sum_processed_study.data.imag = zeros(1,1,study.np/2);
    
            sum_processed_study.data.real(1,1,:) = real(sumfid);
            sum_processed_study.data.imag(1,1,:) = imag(sumfid);
    
            sum_processed_study.multiplicity = 1;
            sum_processed_study.process.lsfid = 0;
            sum_processed_study.process.apodparam1 = 0;
            sum_processed_study.process.apodparam2 = 0;
            sum_processed_study.process.phasecorr0 = 0;
            sum_processed_study.process.phasecorr1 = 0;
            sum_processed_study.process.B0 = zeros(1,study.np/2);
            
        else % Working with fid files
            fidmocor = squeeze(study.data.real).'+1i*squeeze(study.data.imag).';
            grpdly = round(study.params.grpdly) + 1;
            fidmocor = [fidmocor(grpdly:end),zeros(1,grpdly-1)];
            sumfid = fidmocor;

            processed_study = study;

            processed_study.data.real = zeros(1,1,processed_study.np/2);
            processed_study.data.imag = zeros(1,1,processed_study.np/2);

            processed_study.data.real(1,1,:) = real(sumfid);
            processed_study.data.imag(1,1,:) = imag(sumfid);

            processed_study.multiplicity = 1;
            processed_study.process.lsfid = 0;
            processed_study.process.apodparam1 = 0;
            processed_study.process.apodparam2 = 0;
            processed_study.process.phasecorr0 = 0;
            processed_study.process.phasecorr1 = 0;
            processed_study.process.B0 = zeros(1,processed_study.np/2);

            sum_processed_study = processed_study;
        end

        filename = study.filename(1:end-4);
        processed_study.filename = [filename '_processed.mat'];
        processed_study.liststring = fullfile(obj.result_dir,obj.foldername, ...
            'processed',processed_study.filename);
        
        obj.SVS_struct(i).processed_study = processed_study;

        if(~exist(fullfile(obj.result_dir,obj.foldername,'processed'),"dir"))
            mkdir(fullfile(obj.result_dir,obj.foldername,'processed'));
        end
        save(processed_study.liststring,'processed_study');

        filename = study.filename(1:end-4);
        sum_processed_study.filename = ['SUM_' filename '_processed.mat'];
        sum_processed_study.liststring = fullfile(obj.result_dir,obj.foldername, ...
            'processed','sum',sum_processed_study.filename);

        obj.SVS_struct(i).sum_processed_study = sum_processed_study;

        if(~exist(fullfile(obj.result_dir,obj.foldername,'processed','sum'),"dir"))
            mkdir(fullfile(obj.result_dir,obj.foldername,'processed','sum'));
        end
        save(sum_processed_study.liststring,'sum_processed_study');

    catch ME
        msg = {'Error while doing the processing on the experiment : ',obj.SVS_struct(i).exp_name,  ...
            'Error message : ',ME.message};
    end
end
