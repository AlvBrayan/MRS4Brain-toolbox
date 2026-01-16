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
function msg = preprocessing_SVS_FidA_isison_smallvoxel(obj,prog_dbox)
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

            % Combine shots --> specifially when the voxel is small - it is then
            % difficult to align individual shots because the ones with slice
            % inversion have almost 0 signal (same inner and outer volume - due to poor OVS transition bands) 
            out0c=out0b;
            out0c.fids=out0a.fids+out0b.fids; 
            out0c.specs=fftshift(fft(out0c.fids.',[],2),2).';
            
            %apply LB
            dw=out0c.dwelltime;
            tt=[0:dw:dw*(out0c.n-1)];
            out0clb=out0c;
            fids0clb=out0clb.fids.*repmat(exp(-tt*pi*svs_param.LBall).',1,out0c.averages);
            out0clb.fids=fids0clb;
            out0clb.specs=fftshift(fft(out0clb.fids.',[],2),2).';
    
            %% 2-align av.

            if svs_param.FidA_align_avg
                [out1clb,fsc,phsc] = obj.alignAverages_FidA(out0clb, ...
                    svs_param.minppm,svs_param.maxppm,svs_param.t_max,'y');%5.2,9,0.5,'y');
                % fs        = Vector of frequency shifts (in Hz) used for alignment.
                % phs       = Vector of phase shifts (in degrees) used for alignment.

                %remove LB
                out1c=out1clb;
                fids1c=out1clb.fids.*repmat(exp(tt*pi*svs_param.LBall).',1,out1clb.averages);
                out1c.fids=fids1c;
                out1c.specs=fftshift(fft(out1c.fids.',[],2),2).';
            else
                fsc = zeros(round(study.multiplicity/2),1); phsc = fsc;
                out1c = out0c;
            end
    
            %% 3-outlier removal
            if svs_param.FidA_rm_badavg

                [out2c,metricc,badAveragesc] = obj.rmbadaverages_FidA(out1c,svs_param.sd_thresh,'f'); %performs 10Hz LB inside
        
                %apply LB

                out2clb=out2c;
                fids2clb=out2c.fids.*repmat(exp(-tt*pi*svs_param.LBall).',1,out2c.averages);
                out2clb.fids=fids2clb;
                out2clb.specs=fftshift(fft(out2clb.fids.',[],2),2).';
                if ~isempty(badAveragesc)
                    fprintf(2,['!!! ATTENTION !!! Noisy shots: ',num2str(badAveragesc(:)'),' ---> REMOVED ☹☹ \n'])
                    txt  = vertcat({obj.SVS_struct(i).exp_name},{''},{['!!! ATTENTION !!! Noisy shots: ',num2str(badAveragesc(:)'),' ---> REMOVED ☹☹']});
                    tdir = fullfile(obj.result_dir,obj.foldername,['Removed_shots_in_',obj.SVS_struct(i).exp_name,'.txt']); writecell(txt,tdir);   
                else
                    fprintf(2,'!!! ATTENTION !!! No noisy shots 😊😊 \n')
                    txt  = vertcat({obj.SVS_struct(i).exp_name},{''},{'!!! ATTENTION !!! No noisy shots 😊😊'});
                    tdir = fullfile(obj.result_dir,obj.foldername,['No_removed_shots_in_',obj.SVS_struct(i).exp_name,'.txt']); writecell(txt,tdir);    
                end
            else
                out2c=out1c;
                out2c.fids=conj(out2c.fids);
                out2c.specs=fftshift(ifft(out2c.fids,[],out2c.dims.t),out2c.dims.t);
                metricc = zeros(round(study.multiplicity/2),1);
                badAveragesc = [];
            end

            %combine on/off
            fidc=out1c.fids.'; 
            fidtot(1:1:size(fidc,1),:)=fidc;
            fid2sum=fidtot; % contains corrected data but not filtered for outliers

            ind = 1;
            fidmocor = zeros(size(fid2sum,1),size(fid2sum,2));
            for k = 1:size(fidtot,1)/2
                if ~ismember(k,badAveragesc)
                    fidmocor(ind,:) = fid2sum(k,:);
                    ind = ind + 1;
                end
            end
            fidmocor(ind:end,:) = [];
            
            fidmocor = conj(fidmocor);
    
            %% 4-add all the info to the Matlab study structure

            processed_study = study;
            processed_study.fidaprocess.phsc = phsc;
            processed_study.fidaprocess.fsc = fsc;
            processed_study.fidaprocess.metricc = metricc;
            processed_study.fidaprocess.badAveragesc = badAveragesc;
            % processed_study.params.nt = size(fidmocor,1);
            %JM 18032024
            processed_study.params.nt = size(fidmocor,1)*2;
            processed_study.multiplicity = size(fidmocor,1)*2; % changed by BA 30/07/2024

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
        
        clear fid2sum fidtot

    catch ME
        msg = {'Error while doing the processing on the experiment : ',obj.SVS_struct(i).exp_name,  ...
            'Error message : ',ME.message};
    end
end
