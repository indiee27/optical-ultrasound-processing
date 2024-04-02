%% Script to plot the resolution of reconstructed images 
% - tungsten wire phantom
% - fast pullback imaging
%% 

%mkdir res_figs 

q = figure;
set(gcf, 'color', 'white');

% filtering
nwin = 40;
filt_high = 0.05;
filt_low = 0.5;

list = ['1' '2' '4' '8'];

c = 1500;
points = 10;
dB = -35;
   
%% PLOT FIRST FIG
f_path = 'scan_3';
    
% load files
t_axis = importdata([f_path '_taxis.txt']);
hdr_info = importdata([f_path '_hdr.txt'], '\t');  
tdata = importdata([f_path '_data.txt']);

samples = hdr_info.data(1);
a_lines = hdr_info.data(2);
rep_rate = hdr_info.data(3);
scan_length = hdr_info.data(4); 
pullback_speed = hdr_info.data(5);

% params
dt = t_axis(2) - t_axis(1);
scan_depth = round(((samples*1500*1e-8)/2)*1000,0);
step = scan_length/a_lines; dy = step*1e-3;
ratio = round(step*(1/(scan_depth/samples)));

[b,a] = butter(4, filt_high, 'high');
ts = filtfilt(b, a, tdata);
[b,a] = butter(4, filt_low, 'low');
ts = filtfilt(b, a, ts);  

if ~mod(nwin,2)
   nwin = nwin + 1;
end

% set matrices
Y_new = zeros(size(ts));
Y_means = zeros(size(ts));
Y_mod = zeros(size(ts));
T2s = zeros(4,size(ts,2));

for i = 1:size(ts,2)
    if i <= size(ts,2) - nwin+1
        avg_ind_min = i;
        avg_ind_max = i + nwin - 1;
    else
        avg_ind_min = size(ts,2) - nwin+1;
        avg_ind_max = size(ts,2);
    end

    % calculate average
    Y_mean = mean(ts(:,avg_ind_min:avg_ind_max),2);

    % derivative for temporal shift
    dY_mean = [diff(Y_mean);0];
    dY2_mean = [diff(dY_mean);0];

    X2 = [ones(size(ts,1),1), Y_mean, dY_mean, dY2_mean];
    %X2 = [ones(size(ts,1),1), Y_mean, dY_mean];

    % least-squares fit
    T2 = pinv(X2'*X2)*X2'*ts(:,i);

    % subtract fit
    Y_new(:,i)=ts(:,i)-X2*T2;

    % store variables
    Y_means(:,i)=Y_mean;
    Y_mod(:,i)=X2*T2;
    T2s(:,i)=T2; 
end

% create image
forrecon = zeros(size(Y_new,1), size(Y_new,2)+400);
forrecon(:,201:size(Y_new,2)+200) = Y_new;

p_xy = kspaceLineRecon(forrecon, dy*2, dt, c);
img = 20*log10(abs(hilbert(p_xy(:,201:end-200))));

set(0, 'CurrentFigure', q)
subplot(4,2,[1 4])
set(gca, 'box', 'off', 'fontsize', 12, 'fontname', 'arial')
imagesc(img(1:round(3*samples/4),:));
CA = caxis;
caxis([dB 0] + max(CA));
xticks([1 a_lines/2 a_lines-1])
xticklabels({0 scan_length/2 scan_length})
xlabel('Length [mm]')
yticks([1 samples/4 samples/2 3*samples/4 samples-1])
yticklabels({0 scan_depth/4 scan_depth/2 3*scan_depth/4 scan_depth})
ylabel('Depth [mm]')
colormap hot
colorbar
set(gca, 'DataAspectRatio', [1 ratio 1]) 
   
%% resolution measurements
results = cell(1,length(list));
for s = 1:4
       
   d = str2double(list(s));

   result = zeros(4, points);
   a_lines = hdr_info.data(2)/d;
   step = scan_length/a_lines; dy = step*1e-3;
   
   tseries = zeros(samples, a_lines-1);
   for g = 1:(a_lines-1)
       h = 1 + (g - 1)*d;
       tseries(:,g) = tdata(:,h);
   end

   % params
   [b,a] = butter(4, filt_high, 'high');
   ts = filtfilt(b, a, tseries);
   [b,a] = butter(4, filt_low, 'low');
   ts = filtfilt(b, a, ts);

   if ~mod(nwin,2)
       nwin = nwin + 1;
   end

   % set matrices
   Y_new = zeros(size(ts));
   Y_means=zeros(size(ts));
   Y_mod=zeros(size(ts));
   T2s=zeros(4,size(ts,2));

   for i = 1:size(ts,2)
        if i <= size(ts,2) - nwin+1
            avg_ind_min = i;
            avg_ind_max = i + nwin - 1;
        else
            avg_ind_min = size(ts,2) - nwin+1;
            avg_ind_max = size(ts,2);
        end

        % calculate average
        Y_mean = mean(ts(:,avg_ind_min:avg_ind_max),2);

        % derivative for temporal shift
        dY_mean = [diff(Y_mean);0];
        dY2_mean = [diff(dY_mean);0];

        X2 = [ones(size(ts,1),1), Y_mean, dY_mean, dY2_mean];
        %X2 = [ones(size(ts,1),1), Y_mean, dY_mean];

        % least-squares fit
        T2 = pinv(X2'*X2)*X2'*ts(:,i);

        % subtract fit
        Y_new(:,i)=ts(:,i)-X2*T2;

        % store variables
        Y_means(:,i)=Y_mean;
        Y_mod(:,i)=X2*T2;
        T2s(:,i)=T2; 
   end

   % create image
   forrecon = zeros(size(Y_new,1), size(Y_new,2)+400);
   forrecon(:,201:size(Y_new,2)+200) = Y_new;

   p_xy = kspaceLineRecon(forrecon, dy*2, dt, c);
   img = 20*log10(abs(hilbert(p_xy(:,201:end-200))));

   % find noise floor
   noise_img = abs(hilbert(p_xy));
   noise_floor = noise_img(100:200, 200:300);
   noise = mean(mean(noise_img));

   for n = 1:points
       hFig = figure;
       set(gcf, 'Position', get(0, 'Screensize'));
       imagesc(img)
       title(['Depth ' num2str(n)]);

       % find position
       h2 = drawrectangle;
       pos_recon = h2.Position;
       pos_recon = round(pos_recon);

       rawdata_recon = 20*log(abs(hilbert(p_xy(:,201:end-200))));
       ROI = rawdata_recon(pos_recon(2):(pos_recon(2)+pos_recon(4)),pos_recon(1):(pos_recon(1)+pos_recon(3)));

       signal_axial = max(ROI');
       signal_lateral = max(ROI', [], 2);

       [val, peak_loc] = max(signal_lateral);
       signal_norm_lateral = signal_lateral - val;

       [val_2, peak_loc_2] = max(signal_axial);
       signal_norm_axial = signal_axial - val_2;

       res_depth = (peak_loc + pos_recon(2)) * ((dt * 1e6 * c/2)/1000);
       result(1,n) = res_depth;

       ROI_SNR = noise_img(pos_recon(2):(pos_recon(2)+pos_recon(4)),pos_recon(1):(pos_recon(1)+pos_recon(3)));
       S = maxND(ROI_SNR);
       SNR = abs(20*log10(S/noise));

       % find the lateral FWHM
       dB_point_pos1 = find(signal_norm_lateral(peak_loc:end) <=-6, 1, 'first') - 1;
       dB_point_neg1 = find(fliplr(signal_norm_lateral(1:peak_loc)) <=-6, 1, 'first') - 1;

       % find the axial FWHM
       dB_point_pos2 = find(signal_norm_axial(peak_loc_2:end) <= -6, 1, 'first') - 1;
       dB_point_neg2 = find(fliplr(signal_norm_axial(1:peak_loc_2)) <= -6, 1, 'first') -1 ;

       FWHM_lateral = dB_point_pos1 + dB_point_neg1;
       FWHM_axial = dB_point_pos2 + dB_point_neg2;

       R_axial = 1e6*FWHM_axial*(dt*c/2);
       R_lateral = 1e6*FWHM_lateral*dy;

       % add to results matrix
       result(2,n) = R_axial;
       result(3,n) = R_lateral;
       result(4,n) = SNR;

       results{s} = result;
       close(hFig)
   end

   clear tseries
   clear Y_new
   clear Y_means
   clear Y_mod
   clear T2s
   clear forrecon
  
end
   
%% SNR subplot
color_list = {'r', 'b', 'k', 'g'};

% SNR figure
set(0, 'CurrentFigure', q)
subplot(4,2,[5 6])
set(gca,'box', 'off','fontsize',12,'fontname','arial')
for z = 1:4
    SNR_plot = results{z};
    scatter(SNR_plot(1,1:n), smoothdata(smoothdata(SNR_plot(4,1:n))), color_list{z})
    hold on
end
xlabel('Depth [mm]')
ylabel('SNR [dB]')


%% Res subplots
set(0, 'CurrentFigure', q)

% subplot axial
subplot(4,2,7)
set(gca,'box', 'off','fontsize',12,'fontname','arial')
for z = 1:4
    ax_plot = results{z};
    scatter(ax_plot(1,1:n), smoothdata(smoothdata(smoothdata(ax_plot(2,1:n)))), color_list{z})
    hold on
end
xlabel('Depth [mm]')
ylabel('Axial Resolution [\mum]')
%ylim([0 200])

% subplot lateral
subplot(4,2,8)
set(gca,'box', 'off','fontsize',12,'fontname','arial')
for z = 1:4
    lat_plot = results{z};
    scatter(lat_plot(1,1:n), smoothdata(smoothdata(smoothdata(lat_plot(3,1:n)))), color_list{z})
    hold on
end
xlabel('Depth [mm]')
ylabel('Lateral Resolution [\mum]')
%ylim([0 500])

figure_name = (['/res_figs/' f_path '_spacing_' num2str(step*1e3) '.png']);
saveas(gcf, [pwd figure_name]);
figure_name = (['/res_figs/' f_path '_spacing_' num2str(step*1e3) '.fig']);
saveas(gcf, [pwd figure_name]);