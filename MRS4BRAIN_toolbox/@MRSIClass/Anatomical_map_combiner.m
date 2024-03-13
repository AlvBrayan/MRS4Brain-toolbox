% Copyright All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, MRS4Brain research group @ CIBM MRI EPFL AIT, 2024
% See the LICENSE.TXT file for more details.

function f = Anatomical_map_combiner(obj,niiimage,N_slice,mask,interpolation,color,name)

% niiimage = numerical matrix containing the anatomical image
% N_slice = number of the slice to display
% mask = Map to superimpose with the image
% interpolation = boolean value to confirm if there is an interpolation
% color = colormap for the mask (look at
% https://ch.mathworks.com/help/matlab/ref/colormap.html if there is an
% error)
% title = title of the plot

test = rot90(niiimage(:,:,N_slice),-1);

map = mask(end:-1:1,:);

if(interpolation)
    [X,Y]=meshgrid(1:Mat_size(1),1:Mat_size(2));
    [XI,YI]=meshgrid(1:(30/245):Mat_size(1),1:(30/245):Mat_size(2));
    lin=interp2(X,Y,map,XI,YI,'cubic');
else
    lin = kron(map,ones(8,8));
end


BUpperLeftCorner = [3,3];

BB = padarray(lin,BUpperLeftCorner-1,NaN,'pre');
BB = padarray(BB,size(niiimage(:,:,N_slice))-size(BB),NaN,'post');

f = figure('Name',name);
%plot first data 
% title(name)
ax1 = axes; 
im = imagesc(ax1,test); 
im.AlphaData = 1; % change this value to change the background image transparency 
axis square; 
hold all; 
%plot second data 
ax2 = axes; 
im1 = imagesc(ax2,BB); 
im1.AlphaData = 0.45; % change this value to change the foreground image transparency 
axis square;
%link axes 
linkaxes([ax1,ax2]) 
%%Hide the top axes 
ax2.Visible = 'off'; 
ax2.XTick = []; 
ax2.YTick = []; 
%add differenct colormap to different data if you wish 
colormap(ax1,'gray') 
colormap(ax2,color)



%set the axes and colorbar position 
cb2 = colorbar(ax2,'Position',[.88 .11 .0675 .815]); 
set(findall(gcf,'-property','FontSize'),'FontSize',18)
set(findall(gcf,'-property','FontWeight'),'FontWeight','bold')
