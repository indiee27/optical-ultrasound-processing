
mkdir images
list = ['1', '2', '4', '8'];
c = 1500;
dB = -40;
nwin=40;

for x = 3
% load the data
t = (['scan_' num2str(x)]);
tdata = importdata([t '_data.txt']);

taxis = importdata([t '_taxis.txt']);
dt = (taxis(2)-taxis(1));

% header info
hdr_info = importdata([t '_hdr.txt'], '\t'); 
samples = hdr_info.data(1);
a_lines = hdr_info.data(2);
rep_rate = hdr_info.data(3);
scan_length = hdr_info.data(4); 
pullback_speed = hdr_info.data(5);

for s = 1:4

d = str2double(list(s));

% sizes
a_lines = hdr_info.data(2)/d;
step = scan_length/a_lines; dy = step*1e-3;
scan_depth = round(((samples*1500*1e-8)/2)*1000,0);
ratio = round(step*(1/(scan_depth/samples)));

% data acquisition
tseries = zeros(samples, a_lines-1);
for g = 1:(a_lines-1)
    h = 1 + (g - 1)*d;
    tseries(:,g) = tdata(:,h);
end

% low pass filter
[b,a]=butter(4,0.05,'high');
ts=filtfilt(b,a,tseries);
%ts = smooth(ts);
[b,a]=butter(4,0.5,'low');
ts=filtfilt(b,a,ts);

if ~mod(nwin,2)
    nwin=nwin+1;
end
%Y_new = ts;
Y_new=zeros(size(ts));
Y_means=zeros(size(ts));
Y_mod=zeros(size(ts));
T2s=zeros(4,size(ts,2));
for i=1:size(ts,2)
    
    % indices with which to calculate average
    if i<=size(ts,2)-nwin+1
        avg_ind_min=i;
        avg_ind_max=i+nwin-1;
    else
        avg_ind_min=size(ts,2)-nwin+1;
        avg_ind_max=size(ts,2);
    end
    
    % calculate average
    Y_mean=mean(ts(:,avg_ind_min:avg_ind_max),2);
    
    % derivative for temporal shift
    dY_mean=[diff(Y_mean);0];
    dY2_mean=[diff(dY_mean);0];
    
    X2=[ones(size(ts,1),1),Y_mean,dY_mean,dY2_mean];
    %X2=[ones(size(ts,1),1),Y_mean,dY_mean];
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

p_xy = kspaceLineRecon(forrecon,dy*2,dt,c);
img = 20*log10(abs(hilbert(p_xy)));

figure(1);
%subplot(2,1,1);
%imagesc(20*log10(abs(hilbert(forrecon))));
%CA = caxis;
%caxis([dB 0]+max(CA));
%subtitle('Raw data (log+hilbert)');
%xticks([200 200+a_lines/2 200+a_lines])
%xticklabels({0 scan_length/2 scan_length})
%xlabel('Length [mm]')
%yticks([1 samples/2 samples])
%yticklabels({0 scan_depth/2 scan_depth})
%ylabel('Depth [mm]')
%colorbar
%subplot(2,1,2);
imagesc(img(1:3*samples/4,201:end-200));
CA = caxis;
caxis([dB 0]+max(CA));
xticks([1 a_lines/2 a_lines])
xticklabels({0 scan_length/2 scan_length})
xlabel('Length [mm]')
yticks([1 samples/4 samples/2 samples])
yticklabels({0 scan_depth/4 scan_depth/2 scan_depth})
ylabel('Depth [mm]')
title(['Scan ' num2str(x) ': speed: ' num2str(pullback_speed) ' mm/s, spacing: ' num2str(dy*1e6) ' \mu m']);
subtitle('Reconstructed data (log+hilbert)');
colormap(hot);
colorbar
set(gca, 'DataAspectRatio', [1 d*ratio 1])
figure_name = (['/images/' t ' spacing ' num2str(step*1e3) '.png']);
saveas(gcf,[pwd figure_name]);

%close(1)

clear tseries
clear Y_means
clear Y_mod
clear Y_new
clear T2s
clear forrecon
end
end
