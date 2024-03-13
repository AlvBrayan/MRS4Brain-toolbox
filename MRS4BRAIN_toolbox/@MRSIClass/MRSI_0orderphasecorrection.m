% Copyright All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, MRS4Brain research group @ CIBM MRI EPFL AIT, 2024
% See the LICENSE.TXT file for more details.

function [corrected_MRSI_frr, phase_map_rr] = MRSI_0orderphasecorrection(obj,MRSI_frr,ppm,ppm_max,ppm_min,border,err,eta)

N_f = size(MRSI_frr,1);
N_x = size(MRSI_frr,2);
N_y = size(MRSI_frr,3);

corrected_MRSI_frr = zeros(N_f,N_x,N_y);
phase_map_rr = zeros(N_x,N_y);
integ_map_rr = zeros(N_x,N_y);


for x=1:N_x
    for y=1:N_y
        if(border(x,y))
            [corr_spect,phase0,integ_val] = obj.order0phasecorrection(squeeze(MRSI_frr(:,x,y)),ppm,ppm_max,ppm_min,err,eta);
            corrected_MRSI_frr(:,x,y) = corr_spect;
            phase_map_rr(x,y) = phase0;
            integ_map_rr(x,y) = integ_val;

        end
    end
end


%Find the 5 maximum values using integ_map
N_find = 5;
ind_x = [];
ind_y = [];

temp_map = integ_map_rr;
while(N_find > 0)
    [test_x,test_y] = find(max(max(temp_map))==temp_map);
    ind_x = [ind_x,test_y];
    ind_y = [ind_y,test_x];

    temp_map(test_x,test_y) = 0;
    N_find = N_find-1;
end

for ii=1:length(ind_x)
    old_data = squeeze(MRSI_frr(:,ind_y(ii),ind_x(ii)));
    new_data = squeeze(corrected_MRSI_frr(:,ind_y(ii),ind_x(ii)));

    old_plot = figure();
    plot(ppm,real(old_data))
    xlim([0 4])
    set(gca,'Xdir','reverse')
    title(['Original : Column = ',num2str(ind_x(ii)),' / Line = ',num2str(ind_y(ii))])
    xlabel(['ppm'])
    set(findall(gcf,'-property','FontSize'),'FontSize',14)
    set(findall(gcf,'-property','FontWeight'),'FontWeight','bold')

    new_plot = figure();
    plot(ppm,real(new_data))
    xlim([0 4])
    set(gca,'Xdir','reverse')
    title(['Phase Corrected ( ' num2str(phase_map_rr(ind_y(ii),ind_x(ii))) ' deg) : Column = ',num2str(ind_x(ii)),' / Line = ',num2str(ind_y(ii))])
    xlabel(['ppm'])
    set(findall(gcf,'-property','FontSize'),'FontSize',14)
    set(findall(gcf,'-property','FontWeight'),'FontWeight','bold')

    old_plot.Position(1:3) = [500 ,525, 900];
    new_plot.Position(1:3) = [500 , 50, 900];

    pause(0.5);

    close(old_plot);
    close(new_plot);
end
    
end