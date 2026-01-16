%%  preprocessing_DWS_FidA_isisoff.m 

% Copyright All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, MRS4Brain research group @ CIBM MRI EPFL AIT, 2024
% Copyright 2020 Jamie Near
% See the LICENSE.TXT file for more details.

% Jessie Mosso, CIBM - MRS4Brain group, LIFMET, 2021
% Guillaume Briand, CIBM - MRS4Brain group, 2023
% 
% USAGE : DWS Class public method 
% out = obj.preprocessing_DWS_FidA_isisoff(prog_dbox)
% 
% DESCRIPTION :
% Apply FidA preprocessing with ISIS on
%
% INPUTS :
% obj       = DWS Class object with properties and methods
% prog_dbox = MRS4Brain Toolbox progress dialog box
%
% OUTPUT :
% msg       = Error message
function msg = preprocessing_DWS_FidA_isisoff(obj,prog_dbox)
msg = {''};
lcmodelfac = 1; % PARAMETER
dws_param = obj.DWS_param;

for i = 1:length(obj.DWS_struct)
    prog_dbox.Message = ['Apply preprocessing for : ',obj.DWS_struct(i).exp_name];
    prog_dbox.Value = i/(length(obj.DWS_struct)+1);
    try
        study = obj.DWS_struct(i).raw_study;
        if study.multiplicity > 1
            %% step 1 - do preprocessing with FIDA
            %% 1-Convert Bruker study to FID A structure
            out0 = obj.convert2FidA_isisoff(i); % Convert study struct to FidA struct
            % apply Line Broadening (LB)
            dw = out0.dwelltime;
            tt = 0:dw:dw*(out0.n-1);
            out0lb = out0;
            fids0lb = out0lb.fids.*repmat(exp(-tt*pi*dws_param.LBall).',1,out0.averages);
            out0lb.fids = fids0lb;
            out0lb.specs = fftshift(fft(out0lb.fids.',[],2),2).';
    
            %% 2-align av.
            if dws_param.FidA_align_avg
                [out1lb,fs,phs] = obj.alignAverages_FidA(out0lb, ...
                    dws_param.minppm,dws_param.maxppm,dws_param.t_max,'y');%5.2,9,0.5,'y');
        
                % remove Line Broadening (LB)
                out1 = out1lb;
                fids1 = out1lb.fids.*repmat(exp(tt*pi*dws_param.LBall).',1,out1lb.averages);
                out1.fids = fids1;
                out1.specs = fftshift(fft(out1.fids.',[],2),2).';
            else
                fs = zeros(study.multiplicity,1); phs = fs;
                out1 = out0;
            end
    
    
            %% 3-outlier removal
            if dws_param.FidA_rm_badavg
                [out2,metric,badAverages] = obj.rmbadaverages_FidA(out1,dws_param.sd_thresh,'f'); %performs 10Hz LB inside
    
                %apply LB
                out2lb = out2;
                fids2lb = out2.fids.*repmat(exp(-tt*pi*dws_param.LBall).',1,out2.averages);
                out2lb.fids = fids2lb;
                out2lb.specs = fftshift(fft(out2lb.fids.',[],2),2).';
                if ~isempty(badAverages)
                    fprintf(2,['!!! ATTENTION !!! Noisy shots: ',num2str(badAverages(:)'),' ---> REMOVED ☹☹ \n'])
                    txt  = vertcat({obj.DWS_struct(i).exp_name},{''},{['!!! ATTENTION !!! Noisy shots: ',num2str(badAverages(:)'),' ---> REMOVED ☹☹']});
                    tdir = fullfile(obj.result_dir,obj.foldername,['Removed_shots_in_',obj.DWS_struct(i).exp_name,'.txt']); writecell(txt,tdir);    
                else
                    fprintf(2,'!!! ATTENTION !!! No noisy shots 😊😊 \n')
                    txt  = vertcat({obj.DWS_struct(i).exp_name},{''},{'!!! ATTENTION !!! No noisy shots 😊😊'});
                    tdir = fullfile(obj.result_dir,obj.foldername,['No_removed_shots_in_',obj.DWS_struct(i).exp_name,'.txt']); writecell(txt,tdir);  
                end
            else
                metric = zeros(study.multiplicity,1);
                badAverages = [];
                out2 = out1;
            end
    
            %% 4-add all the info to the Matlab study structure
            processed_study = study;
            processed_study.fidaprocess.phs = phs;
            processed_study.fidaprocess.fs = fs;
            processed_study.fidaprocess.metric = metric;
            processed_study.fidaprocess.badAverages = badAverages;
            processed_study.params.nt = size(out2.fids,2);
            processed_study.multiplicity = size(out2.fids,2);
            processed_study.process.apodparam1 = zeros(1,size(out2.fids,2));
            processed_study.process.apodparam2 = zeros(1,size(out2.fids,2));
            processed_study.process.phasecorr0 = zeros(1,size(out2.fids,2));
            processed_study.process.phasecorr1 = zeros(1,size(out2.fids,2));
            processed_study.data.real = zeros(size(out2.fids,2),1,size(out2.fids,1));
            processed_study.data.real(:,1,:) = real(out2.fids.');
            processed_study.data.imag = zeros(size(out2.fids,2),1,size(out2.fids,1));
            processed_study.data.imag(:,1,:) = imag(out2.fids.');

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
            grpdly = round(study.params.grpdly)+1;
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
        
        obj.DWS_struct(i).processed_study = processed_study;

        if(~exist(fullfile(obj.result_dir,obj.foldername,'processed'),"dir"))
            mkdir(fullfile(obj.result_dir,obj.foldername,'processed'));
        end
        save(processed_study.liststring,'processed_study');

        filename = study.filename(1:end-4);
        sum_processed_study.filename = ['SUM_' filename '_processed.mat'];
        sum_processed_study.liststring = fullfile(obj.result_dir,obj.foldername, ...
            'processed','sum',sum_processed_study.filename);

        obj.DWS_struct(i).sum_processed_study = sum_processed_study;

        if(~exist(fullfile(obj.result_dir,obj.foldername,'processed','sum'),"dir"))
            mkdir(fullfile(obj.result_dir,obj.foldername,'processed','sum'));
        end
        save(sum_processed_study.liststring,'sum_processed_study');

    catch ME
        msg = {'Error while doing the processing on the experiment : ',obj.DWS_struct(i).exp_name,  ...
            'Error message : ',ME.message};
    end
end
