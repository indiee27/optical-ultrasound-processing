%% Load in Data
fname = 'speed60ms_reprate600hz_length20mm';

t_data = importdata([fname '_data.txt']);
t_axis = importdata([fname '_taxis.txt']);
%t_axis = linspace(10e-9,26.67e-6,2667);
hdr_inf = importdata([fname '_hdr.txt'], '\t');
hdr_inf.data(1)
t_data = t_data';
dt = t_axis(2)-t_axis(1);
dy = 50e-6;
c = 1500;
%%
% high pass filter
    
[b,a]=butter(4,0.2,'high');     % high pass filter
ts=filtfilt(b,a,t_data);

nwin=40;%20-60
if ~mod(nwin,2)
    nwin=nwin+1;
end
Y_new = ts;Y_new=zeros(size(ts));
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


forrecon = zeros(size(Y_new,1),size(Y_new,2)+400);
forrecon(:,201:size(Y_new,2)+200) = Y_new;

p_xy = kspaceLineRecon(forrecon,dy,dt/2,c);

dB = -30;

%% Recon data ROI selection
result=zeros(2,3);
for n=1:3
%% Raw data ROI selection
figure;
set(gcf,'Position');
imagesc((abs(hilbert(Y_new))));
% CA = caxis;
% caxis([dB 0]+max(CA));
title('Raw data (log+hilbert)');
% colormap(hot);
h1 = imrect;
pos_raw = getPosition(h1);
pos_raw=round(pos_raw);
close gcf;
%%

rawdata=(abs(hilbert(Y_new)));
%rawdata_recon=20*log10(abs(hilbert(p_xy)));
% The 'rawdata' is in dB, i.e. not voltage values, this means the value at
% which the signal drops to half is not 0.5 the value of the peak. Instead
% is it -6 dB, this correspond to half the max value. Also, because of this
% we can't use the Gaussian fit either. Instead you should fit to the
% actual raw data -> abs(hilbert(p_xy)) without the log10+



ROI_raw=rawdata(pos_raw(2):(pos_raw(2)+pos_raw(4)),pos_raw(1):(pos_raw(1)+pos_raw(3)));
%ROI_recon_1=rawdata_recon(pos_recon(2):(pos_recon(2)+pos_recon(4)),pos_recon(1):(pos_recon(1)+pos_recon(3)));


%ROI_raw1=filtfilt(y,x,ROI_raw);
%%

signal_axial = max(ROI_raw');
signal_lateral=max(ROI_raw',[],2)';


[y,x]=butter(5,0.2,'low');   %low pass filter to eliminate high f components
signal_lateral=filtfilt(y,x,signal_lateral);
signal_axial=filtfilt(y,x,signal_axial);


signal_lateral_med = medfilt1(signal_lateral,5); %median filter 
signal_axial_med = medfilt1(signal_axial,5);

[val, peak_loc] = max(signal_lateral);	% location of the mainlobe
signal_norm_lateral = signal_lateral - val;     % normalized signal



[val2, peak_loc2] = max(signal_axial);

signal_norm_axial=signal_axial-val2;

%plot(signal);



%Gaussian fit

%{
f_lateral= fit([1:1:size(ROI,1)]',signal_norm_lateral','gauss4');
[peak_f, maxloc]=max(f_lateral);

bottom_lateral_f=min(f_lateral);
dB_point_pos3 = find(f_lateral(peak_loc:end) <= (bottom_lateral/2), 1, 'first') - 1;
dB_point_neg3 = find(fliplr(f_lateral(1:peak_loc)) <= (bottom_lateral/2), 1, 'first') - 1;
%}

figure(2);
subplot(1,2,1);

plot(signal_norm_lateral,'Color',[0,0.7,0.9]);

title('Lateral signal');
xlabel('Location');
ylabel('Amplitude');

subplot(1,2,2);
plot(signal_norm_axial,'Color',[0,0.7,0.9]);
title('Axial signal');
xlabel('location');
ylabel('Amplitude');

%%

 %find the lateral resolution

dB_point_pos1 = find(signal_norm_lateral(peak_loc:end) <=-6, 1, 'first') - 1;
dB_point_neg1 = find(fliplr(signal_norm_lateral(1:peak_loc)) <=-6, 1, 'first') - 1;

% find the axial resolution
dB_point_pos2 = find(signal_norm_axial(peak_loc2:end) <= -6, 1, 'first') - 1;
dB_point_neg2 = find(fliplr(signal_norm_axial(1:peak_loc2)) <= -6, 1, 'first') -1 ;
%suptitle(['lateral Bottom=',num2str(bottom_lateral),',','axial bottom=',num2str(bottom_axial)])

%%
%  Linear interpolation to find the exact lateral FWHM
dB_point_pos_interp1 = dB_point_pos1 + ...
    (((-6) - signal_norm_lateral(peak_loc+dB_point_pos1)) / (signal_norm_lateral(peak_loc+dB_point_pos1) - signal_norm_lateral(peak_loc+dB_point_pos1-1)));
dB_point_neg_interp1 = dB_point_neg1 + ...
    (((-6) - signal_norm_lateral(peak_loc-dB_point_neg1)) / (signal_norm_lateral(peak_loc-dB_point_neg1) - signal_norm_lateral(peak_loc-dB_point_neg1+1)));

if (numel(dB_point_pos_interp1) == 0 || numel(dB_point_neg1) == 0)
    FWHM_lateral = [];
else
    FWHM_lateral = dB_point_pos_interp1 + dB_point_neg_interp1;
end


%  Linear interpolation to find the exact axial FWHM
dB_point_pos_interp2 = dB_point_pos2 + ...
    (((-6) - signal_norm_axial(peak_loc2+dB_point_pos2)) / (signal_norm_axial(peak_loc2+dB_point_pos2) - signal_norm_axial(peak_loc2+dB_point_pos2-1)));
dB_point_neg_interp2 = dB_point_neg2 + ...
    (((-6) - signal_norm_axial(peak_loc2-dB_point_neg2)) / (signal_norm_axial(peak_loc2-dB_point_neg2) - signal_norm_axial(peak_loc2-dB_point_neg2+1)));

if (numel(dB_point_pos_interp2) == 0 || numel(dB_point_neg2) == 0)
   FWHM_axial = [];
else
   FWHM_axial = dB_point_pos_interp2 + dB_point_neg_interp2;
end
%%




R_axial=1e6*FWHM_axial* (dt*c/2);
R_lateral=1e6*FWHM_lateral*dy;

result(1,n)=R_axial;
result(2,n)=R_lateral;




end