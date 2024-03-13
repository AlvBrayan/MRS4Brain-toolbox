%%  data2RAW.m 

% Copyright All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, MRS4Brain research group @ CIBM MRI EPFL AIT, 2024
% See the LICENSE.TXT file for more details.

% Brayan Alves, CIBM - MRS4Brain group, 2022
% Guillaume Briand, CIBM - MRS4Brain group, 2023
% 
% USAGE : MRSI Class public method 
% msg = obj.data2RAW(data_signal,reference_signal,output_name,slice)
% 
% DESCRIPTION :
% Convert data to RAW files in order to be quantify with LCModel
%
% INPUTS :
% obj           = MRSI Class object with properties and methods
% data_signal   = data fid signal
% ref_signal    = reference fid signal
% output_name   = MRSI slice number (Slice_Nk)
% slice         = slice number
%
% OUTPUT :
% msg       = Error message
function msg = data2RAW(obj,data_signal,reference_signal,output_name,slice)
msg = {''};
folder_name = fullfile(obj.results_folder,output_name);
data_folder = fullfile(folder_name,'Data');

bool_position = squeeze(obj.Final_mask(slice,:,:));
Nx = size(bool_position,1);
Ny = size(bool_position,2);

try
    %% Data save  
    % Save reference data
    for i = 1:Nx
        for j = 1:Ny
            if(bool_position(i,j))
                filename = [append(output_name+"@"+i+"_"+j+"w.RAW")]; % give it a name
                filename = fullfile(data_folder,filename);
                sdata = length(reference_signal(:,i,j));   %size
                NFFT = 2^nextpow2((sdata));
                RAW.real = zeros([NFFT 1]);
                RAW.imag = zeros([NFFT 1]);
                RAW.real(1:sdata) = real(reference_signal(:,i,j)); % fill in your data
                RAW.imag(1:sdata) = imag(reference_signal(:,i,j)); % fill in your data
    
                ID = '';
                rawfile = filename;
                tramp = 1;
                volume = 1;
    
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
        end
    end
    
    % Save metabolite data
    for i = 1:Nx
        for j=1:Ny
            if(bool_position(i,j))
                filename = [append(output_name+"@"+i+"_"+j+".RAW")]; % give it a name
                filename = fullfile(data_folder,filename);
                sdata = length(data_signal(:,i,j));   %size
                NFFT = 2^nextpow2((sdata));
                RAW.real = zeros([NFFT 1]);
                RAW.imag = zeros([NFFT 1]);
                RAW.real(1:sdata) = real(data_signal(:,i,j)); % fill in your data
                RAW.imag(1:sdata) = imag(data_signal(:,i,j)); % fill in your data
    
                ID = '';
                rawfile = filename;
                tramp = 1;
                volume = 1;
    
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
        end
    end
catch ME
    msg = {'Error while transforming data into RAW files',ME.message};
end
end