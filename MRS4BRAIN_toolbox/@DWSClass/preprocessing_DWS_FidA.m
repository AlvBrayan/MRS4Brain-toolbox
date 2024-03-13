%%  preprocessing_DWS_FidA.m 

% Copyright All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, MRS4Brain research group @ CIBM MRI EPFL AIT, 2024
% Copyright 2020 Jamie Near
% See the LICENSE.TXT file for more details.

% Jessie Mosso, CIBM - MRS4Brain group, LIFMET, 2021
% Guillaume Briand, CIBM - MRS4Brain group, 2023
% 
% USAGE : DWS Class public method 
% out = obj.preprocessing_DWS_FidA(prog_dbox)
% 
% DESCRIPTION :
% Apply FidA preprocessing with or without ISIS on
%
% INPUTS :
% obj       = DWS Class object with properties and methods
% prog_dbox = MRS4Brain Toolbox progress dialog box
%
% OUTPUT :
% msg       = Error message
function msg = preprocessing_DWS_FidA(obj,prog_dbox)
% Preprocessing function using FidA from Jamie Near
if obj.DWS_param.ISIS_bool
    msg = obj.preprocessing_DWS_FidA_isison(prog_dbox);
else
    msg = obj.preprocessing_DWS_FidA_isisoff(prog_dbox);
end
end
