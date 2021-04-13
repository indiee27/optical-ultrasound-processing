%% ***** Read data:
filepath = cd;
filebase = '200400_50umstep_61x61';

[data,taxis,~,~,Nt,Npos] = read_grid_scan_data(filepath,filebase);
data = squeeze(data(:,:,2));

% return
data(:,taxis<=1E-6) = 0;



%% ***** Define / read grid parameters:
Nx = 61;
Ny = 61;
dx = 50E-6;
dy = dx;
xaxis = (0:(Nx-1))*dx;
yaxis = (0:(Ny-1))*dy;

dt = taxis(2)-taxis(1);
c = 1480;

fmin = 0E6;
fmax = 125E6;
F = generate_filter(taxis,fmin,fmax,128);
F = ones(Nx*Ny,1)*F;



%% ***** Post-process data:
% Apply hydrophone calibration:
faxis = 1/dt/Nt * ([0:ceil(Nt/2)-1,-floor(Nt/2):-1]);
cal_data_int = get_calibration_data(faxis,'2672');
cal_data_int(abs(faxis)<=1E6) = Inf;
cal_data_int = ones(Nx*Ny,1)*cal_data_int;
data = real(ifft(  fft(data,[],2)./cal_data_int.*F  ,[],2));

% Envelope detection:
data_env = abs(hilbert(data'))';

% Reshape to 3D matrix:
data     = reshape(data,Nx,Ny,Nt);
data_env = reshape(data_env,Nx,Ny,Nt);

% Extract maximum observed in time and space:
max_t    = squeeze(max(data_env,[],3));
[maxx,maxy] = find(max_t==max(max_t(:)));
Ascan    = squeeze(data(maxx,maxy,:));
arrtime  = find(  squeeze(data_env(maxx,maxy,:))==max(squeeze(data_env(maxx,maxy,:))),1);
z0       = (arrtime-1)*dt*c;

mask     = ((1:Nt)>arrtime-50 & (1:Nt)<arrtime+50)';
SPECT    = 20*log10(abs(fft(  Ascan .* mask  )));
faxis    = 1/dt/Nt * [0:ceil(Nt/2)-1,-floor(Nt/2):-1];



%% ***** Display maximum pressure:
figure;colormap hot;
subplot(1,3,[1 2]);
hold on;
    imagesc(xaxis*1000,yaxis*1000,max_t);
    contour(xaxis*1000,yaxis*1000,max_t, [1 1]*max(max_t(:))/2, 'color','g','linewidth',2,'linestyle',':');
hold off;
axis equal tight;
xlabel('x [mm]');
ylabel('y [mm]');
caxis([0 max(max_t(:))]);
title('Maximum observed pressure [MPa]');
colorbar;

subplot(2,3,3);
plot(taxis*1E6,Ascan.*mask);
xlabel('Time [\mus]');
ylabel('Pressure [MPa]');
title('Pressure in location exhibiting peak amplitude')
axis([0 2 -1.5 1.5]);

subplot(2,3,6);
plot(faxis(1:Nt/2)/1E6,SPECT(1:Nt/2),'linewidth',2);
line([0 50],(max(SPECT)-6)*[1 1],'linestyle','--','color','r','linewidth',2);
xlabel('Frequency [MHz]');
ylabel('Power [dB]');
title('Corresponding power spectrum')
axis([0 50 max(SPECT)-[50 0]]);





%% ***** Propagate data along z:
dz = (linspace(0,16,9)-1000*z0)/1000;
pASA_tot = propagate_ASA_multi(data,dz,dx,dy,dt,c,[2 2 1]);

%%
figure;
colormap hot;

for aa = 1:min([length(dz),9])

    toplot = squeeze(pASA_tot(:,:,:,aa));
    toplot = reshape(toplot,Nx*Ny,Nt);
    toplot = abs(hilbert(toplot'))';
    toplot = reshape(toplot,Nx,Ny,Nt);
    toplot = squeeze(max(toplot,[],3));

    subplot(3,3,aa);
    hold on;
    imagesc(xaxis*1000,yaxis*1000,toplot);
    contour(xaxis*1000,yaxis*1000,toplot, [1 1]*max(toplot(:))/2, 'color','g','linewidth',2,'linestyle',':');
    hold off;
    axis equal tight;
    caxis([0 max(toplot(:))]);
    xlabel('x [mm]');
    ylabel('y [mm]');
    title({['Maximum observed pressure [MPa] - distance = ',num2str((dz(aa)+z0)*1000),' mm'] ; ...
           ['FWHM extent x: ',num2str(sum(  max(toplot,[],1) >= max(toplot(:))/2 ) * dx * 1000), ' mm, y: ',num2str(sum(  max(toplot,[],2) >= max(toplot(:))/2 ) * dy * 1000),' mm']});
    colorbar;
end