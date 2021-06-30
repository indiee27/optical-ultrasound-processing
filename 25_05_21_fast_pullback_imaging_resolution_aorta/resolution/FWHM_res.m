% script for FWHM of resolution phantom

%% Import image and data
fname = 'speed60ms_reprate600hz_length20mm';

t_data = importdata([fname '_data.txt']);
t_axis = importdata([fname '_taxis.txt']);
%t_axis = linspace(10e-9,26.67e-6,2667);
hdr_inf = importdata([fname '_hdr.txt'], '\t');

%% Set variables
dt = t_axis(2) - t_axis(1);
Nt = numel(t_axis);
dx = ((hdr_inf.data(5))/(hdr_inf.data(3)))*1e-3;
Nx = hdr_inf.data(2);
dy = 100e-6;
c = 1500;

%% Frequency Filter
% butter function has (n, Wn, type) where n = nth order butterworth, Wn =
% normalized cutoff frequency and type is high or low

% might be interesting to see if you can write a function to calculate the
% best cutoff frequency??

[b,a] = butter(4,0.02,'high');
ts = filtfilt(b,a,t_data);
[b,a] = butter(4,0.7,'low');
ts = filtfilt(b,a,ts);

%% WHAT IS THIS?
nwin = 5; %20-60
if ~mod(nwin,2)
    nwin = nwin+1;
end 

%% Cross Talk Removal
Y_new = ts; %this is pointless?
Y_new = zeros(size(ts));
Y_means = zeros(size(ts));
Y_mod = zeros(size(ts));

T2s = zeros(3, size(ts,2)); %rename this variable it sucks

for i = 1:size(ts,2)
    % indices with which to calculate average
    
%         if i<=size(Y,2)-nwin+1;
%             avg_ind_min=i;
%             avg_ind_max=i+nwin-1;
%         else
%             avg_ind_min=size(Y,2)-nwin+1;
%             avg_ind_max=size(Y,2);
%         end
    avg_ind_min = max(1,i-(nwin-1)/2);
    avg_ind_max = min(size(ts,2),i+(nwin-1)/2);
     
    % calculate average
    Y_mean = mean(ts(:,avg_ind_min:avg_ind_max),2);
    
    % derivative for temporal shift
    dY_mean = [diff(Y_mean);0];
    %dY2_mean=[diff(dY_mean);0];
    
    %X2=[ones(size(Y,1),1),Y_mean,dY_mean,dY2_mean];
    X2 = [ones(size(ts,1),1),Y_mean,dY_mean];
    %X2=[ones(size(Y,1),1),Y_mean];
    
    % least-squares fit
    T2 = pinv(X2'*X2)*X2'*ts(:,i);
    
    % subtract fit
    Y_new(:,i) = ts(:,i)-X2*T2;
    
    % store variables
    Y_means(:,i) = Y_mean;
    Y_mod(:,i) = X2*T2;
    T2s(:,i) = T2;
    
end

%% Reconstruction of Data
% size of recon image

forrecon = zeros(size(Y_new,1),size(Y_new,2)+400);
forrecon(:,201:size(Y_new,2)+200) = Y_new;

%% Image Construction with Buffer
p_xy = kspaceLineRecon(forrecon,dy,dt/2,c);

%% Image Construction with no Buffer

data_rec = p_xy(:,201:end-200);

%% FWHM
% this needs to be on the hilbert transformed data but before it is log
% transformed
dB = -35;

%% Plot Hilbert transformed data
% figure;
% hilbert_image = abs(hilbert(Y_new));
% imagesc(hilbert_image);
% 
figure;
hilbert_image_recon = abs(hilbert(p_xy));
imagesc(hilbert_image_recon);

%% Wire subsections


%% Find the max value

%{
for n=1:5
    signal_axial = max(['wire_' 'n']);
    signal_lateral = max(['wire_' 'n'], [], 2);

%% Plot cross sections - axial and lateral

% figure(2);
% subplot(1,2,1);
% 
% plot(signal_norm_lateral,'Color',[0,0.7,0.9]);
% 
% title('Lateral signal');
% xlabel('Location');
% ylabel('Amplitude');
% 
% subplot(1,2,2);
% plot(signal_norm_axial,'Color',[0,0.7,0.9]);
% title('Axial signal');
% xlabel('location');
% ylabel('Amplitude');

%% Measure FWHM of cross section


%% Save the values for each motor speed

end
%}
