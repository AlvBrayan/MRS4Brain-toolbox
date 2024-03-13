%%  data2RAW.m 

% Copyright All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, MRS4Brain research group @ CIBM MRI EPFL AIT, 2024
% See the LICENSE.TXT file for more details.

% Guillaume Briand, CIBM - MRS4Brain group, 2023
% 
% USAGE : SVS Class public method 
% msg = obj.data2RAW()
% 
% DESCRIPTION :
% Convert processed data from study structure to RAW files in order to be
% quantified by LCModel
%
% INPUTS :
% obj       = SVS Class object with properties and methods
%
% OUTPUT :
% msg       = Error message
function msg = data2RAW(obj)
% Function to convert data to RAW files in order to be quantify with
% LCModel
msg = {''};

folder_name = fullfile(obj.result_dir,obj.foldername,'quantified');
if ~exist(folder_name,"dir")
    mkdir(folder_name);
end
lcmodelraw_name = fullfile(folder_name,'raw_lcmodel');
if ~exist(lcmodelraw_name,"dir")
    mkdir(lcmodelraw_name);
end

%% Data save
% Save metabolite or water data
try
    for i = 1:length(obj.SVS_struct)
        tmp = obj.SVS_struct(i).sum_processed_study;
        filename = tmp.filename(1:end-4);
        filename = fullfile(lcmodelraw_name,[filename,'.RAW']);
    
        if size(tmp.params.nt,1) == 1
            norm_factor = tmp.params.nt;
        else
            norm_factor = 1;
        end
        
        data.real = squeeze(tmp.data.real) ./ norm_factor;
        data.imag = squeeze(tmp.data.imag) ./ norm_factor;
    
        sdata = length(data.real);
        NFFT = 2^nextpow2((sdata));
    
        RAW.real = zeros([NFFT 1]);
        RAW.imag = zeros([NFFT 1]);
        RAW.real = data.real; % fill in your data
        RAW.imag = data.imag; % fill in your data
        
        ID = '';
        rawfile = filename;
        tramp = 1;
        volume = tmp.params.vox1 .* tmp.params.vox2 .* tmp.params.vox3;
        
        %header raw file
        fileid = fopen(rawfile, 'w','b');
        fprintf(fileid,' \n $NMID\n');
        fprintf(fileid,' ID=\''%s\''\n',ID);
        fprintf(fileid,' FMTDAT=\''(2E13.5)\''\n');
        datastring=strrep(sprintf(' TRAMP= % 13.5E\n', tramp), 'E+0', 'E+');
        datastring=strrep(datastring, 'E-0', 'E-');
        fprintf(fileid,datastring);
        datastring=strrep(sprintf(' VOLUME= % 13.5E\n', volume), 'E+0', 'E+');
        datastring=strrep(datastring, 'E-0', 'E-');
        fprintf(fileid,datastring);
        
        
        fprintf(fileid,' $END\n');
        %data
        for k = 1:1:length(RAW.real)
            datastring = sprintf(' % 13.5E % 13.5E\n',RAW.real(k), RAW.imag(k));
            datastring = strrep(datastring, 'E+0', 'E+');
            datastring = strrep(datastring, 'E-0', 'E-');
            fprintf(fileid,datastring);
        
        end
        fclose(fileid);
    end
catch ME
    msg = {'Error while transforming processed data to RAW file :',ME.message};
end
end