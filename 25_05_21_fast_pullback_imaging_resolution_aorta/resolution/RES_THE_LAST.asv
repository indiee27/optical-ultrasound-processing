%% Script for FWHM of resolution wire phantom
% side viewing fibre 2021-05-21
%% Things to do for each scan set
% first run plot_image on one of the scans to find the coords of the wires
% and a noise image, then change the following variables:
% - scan_list
% - f_min and f_max
% - number of wires
% - wire coordinates
% - noise image location
% - ax_dim and lat_dim
% - dy

%% List all the fnames
% NB I've left out the 30mm scans
scan_list = {'speed10ms_reprate100hz_length20mm', 'speed10ms_reprate200hz_length20mm', 'speed15ms_reprate150hz_length20mm', 'speed15ms_reprate300hz_length20mm', 'speed20ms_reprate200hz_length20mm', 'speed20ms_reprate400hz_length20mm', 'speed25ms_reprate250hz_length20mm', 'speed25ms_reprate500hz_length20mm', 'speed30ms_reprate300hz_length20mm', 'speed30ms_reprate600hz_length20mm', 'speed35ms_reprate350hz_length20mm', 'speed35ms_reprate700hz_length20mm', 'speed40ms_reprate400hz_length20mm', 'speed40ms_reprate800hz_length20mm', 'speed45ms_reprate450hz_length20mm', 'speed45ms_reprate900hz_length20mm', 'speed50ms_reprate500hz_length20mm', 'speed50ms_reprate1000hz_length20mm', 'speed60ms_reprate600hz_length20mm', 'speed70ms_reprate700hz_length20mm', 'speed80ms_reprate800hz_length20mm', 'speed90ms_reprate900hz_length20mm', 'speed100ms_reprate1000hz_length20mm'};
s = length(scan_list);

%% constants for all the images
f_min = 0.05;
f_max = 0.5;

% wire constants
% NB the vectors for wires are [column number, row number]
w = 5; % number of wires in the scan
% wire coords
wire_1 = [1153,63]; 
wire_2 = [1146,76];
wire_3 = [534,146];
wire_4 = [535,162];
wire_5 = [544,176];

% list to iterate through
wires = {wire_1, wire_2, wire_3, wire_4, wire_5};

% dimensions
ax_dim = 40;
lat_dim = 10;

% graphing constants
t_axis = linspace(10e-9,26.67e-6,2667);
dt = t_axis(2) - t_axis(1);
dy = 100e-6; % double check this value
c = 1500;
d_ax = dt * 1e6 * c/2;

%% axial and lateral cross section shapes
lat_res_total = zeros(s,w); 
ax_res_total = zeros(s,w); 
ref_amps = zeros(s,w);
res_depth = zeros(s,w);
snrs = zeros(s,w);

%% import image specific data
for n = 1:s

f_name = scan_list{n};
t_data = importdata([f_name '_data.txt']);
hdr_info = importdata([f_name '_hdr.txt'], '\t'); % this needs to be increase
pullback_speed = hdr_info.data(5);
rep_rate = hdr_info.data(3);

%% Step size check

if pullback_speed/rep_rate ~= dy
   msgbox(['iteration incorrect'])
   n
end

%% Low pass filter
[b,a]=butter(4,f_min,'high');
ts=filtfilt(b,a,t_data);
[b,a]=butter(4,f_max,'low');
ts=filtfilt(b,a,ts);

%% cross-talk removal and downsample in time axis
nwin = 50;
if ~mod(nwin,2)
    nwin = nwin + 1;
end

Y_new=zeros(size(ts));
Y_means=zeros(size(ts));
Y_mod=zeros(size(ts));
T2s=zeros(3,size(ts,2));
for i=1:size(ts,2)
    
    % indices with which to calculate average
    %     if i<=size(Y,2)-nwin+1;
    %         avg_ind_min=i;
    %         avg_ind_max=i+nwin-1;
    %     else
    %         avg_ind_min=size(Y,2)-nwin+1;
    %         avg_ind_max=size(Y,2);
    %     end
    avg_ind_min=max(1,i-(nwin-1)/2);
    avg_ind_max=min(size(ts,2),i+(nwin-1)/2);
    
    % calculate average
    Y_mean=mean(ts(:,avg_ind_min:avg_ind_max),2);
    
    % derivative for temporal shift
    dY_mean=[diff(Y_mean);0];
    %dY2_mean=[diff(dY_mean);0];
    
    %X2=[ones(size(Y,1),1),Y_mean,dY_mean,dY2_mean];
    X2=[ones(size(ts,1),1),Y_mean,dY_mean];
    %X2=[ones(size(Y,1),1),Y_mean];
    
    % least-squares fit
    T2=pinv(X2'*X2)*X2'*ts(:,i);
    
    % subtract fit
    Y_new(:,i)=ts(:,i)-X2*T2;
    
    % store variables
    Y_means(:,i)=Y_mean;
    Y_mod(:,i)=X2*T2;
    T2s(:,i)=T2;
    
end

%% recon
forrecon = zeros(size(Y_new,1),size(Y_new,2)+400);
forrecon(:,201:size(Y_new,2)+200) = Y_new;

p_xy = kspaceLineRecon(forrecon,dy,dt/2,c);

img = abs(hilbert(p_xy(:,201:end-200)));
%catim = medfilt2(catim);
figure;
imagesc(img);

% define noise image for SNR
noise_img = img(1050:1450, 65:100);
noise_floor = mean(mean(noise_img));

%% locations and resolutions

% iteration loop
for x = 1:w
    wire = wires{x};
    
    % location coords
    column_num = wire(1);
    row_num = wire(2);
    
    % res img and peak
    res_img = img(column_num:column_num+10, row_num:row_num+10);
    [res_peak, coords] = maxND(res_img);
    
    % updating values for the loop
    % note 1 should become n when the overall loop is added
    snrs(n,x) = 20*log10(res_peak/noise_floor);
    ref_amps(n,x) = res_peak;
    res_depth(n,x) = (coords(1) + column_num) * d_ax/1000;
    
    % axial
    ax_left = squeeze(img(coords(1) + column_num - ax_dim:coords(1) + column_num - 1, coords(2) + row_num - 1));
    ax_right = squeeze(img(coords(1) + column_num:coords(1) + column_num + ax_dim, coords(2) + row_num - 1));

    % lateral
    lat_left = squeeze(img(coords(1) + column_num - 1, coords(2) + row_num - lat_dim:coords(2) + row_num - 1));
    lat_right = squeeze(img(coords(1) + column_num - 1, coords(2) + row_num:coords(2) + row_num + lat_dim));
    
    % ax interp
    ax_left = interp(ax_left,2);
    ax_right = interp(ax_right,2);
    
    % lat interp
    lat_left = interp(lat_left,2);
    lat_right = interp(lat_right,2);
    
    % axial max and min value of FWHM 
    ax_first = find(ax_right < res_peak/2, 1); % first value above half max 
    ax_last = find(ax_left < res_peak/2,1,'last');
    
    % lateral max and min of FWHM
    lat_first = find(lat_right < res_peak/2,1);
    lat_last = find(lat_left < res_peak/2,1,'last');
    
    % resolutions
    lat_res = ((lat_dim*2) - lat_last + lat_first) * 25;
    ax_res = ((ax_dim*2) - ax_last + ax_first) * d_ax/2;
    
    % add to matrix
    lat_res_total(n,x) = lat_res;
    ax_res_total(n,x) = ax_res;
    
    % figure
    figure;
    subplot(1,3,1);
    imagesc(linspace(dy,dy*15,15)*1000,linspace(d_ax,d_ax*30,30)/1000,res_img);
    xlabel('x [mm]');
    ylabel('Depth [mm]');
    title({'Wire image to find', 'resolution'});
    
    ax_plot = [ax_left;ax_right];
    
    subplot(1,3,2);
    plot(linspace(d_ax/2,d_ax*81,162)/1000,ax_plot);
    xlabel('Depth [mm]');
    ylabel('Signal strength [a.u.]');
    title({'Axial cross-section', 'through PSF'});

    lat_plot = [lat_left,lat_right];

    subplot(1,3,3);
    plot(linspace(25,25*21,42)/1000,lat_plot);
    xlabel('x [mm]');
    ylabel('Signal strength [a.u.]');
    title({'Lateral cross-section', 'through PSF'});

end
end

%% plot the final results

figure;
subplot(2,2,1);
plot(ax_res_total');
xlabel('Depth [mm]');
ylabel('Axial Resolution [\mu m]');
title('Axial Resolution');
xticklabels(squeeze(res_depth(1,:)));

subplot(2,2,2);
plot(lat_res_total');
xlabel('Depth [mm]');
ylabel('Lateral Resolution [\mu m]');
title('Lateral Resolution');
xticklabels(squeeze(res_depth(1,:)));

subplot(2,2,3);
plot(ref_amps');
xlabel('Depth [mm]');
ylabel('Signal Peak [a.u.]');
title('Peak Signal Strength');
xticklabels(squeeze(res_depth(1,:)));

subplot(2,2,4);
plot(snrs');
xlabel('Depth [mm]');
ylabel('SNR [dB]');
title('Signal-to-noise Ratio');
xticklabels(squeeze(res_depth(1,:)));


