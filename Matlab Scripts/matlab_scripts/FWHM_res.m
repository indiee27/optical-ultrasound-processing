% script for FWHM of resolution phantom

%% Import image and data
fname = '';

t_data = importdata([fname '_data.txt']);
t_axis = importdata([fname '_taxis.txt']);
%t_axis = linspace(10e-9,26.67e-6,2667);
hdr_inf = importdata([fname '_hdr.txt'], '\t');

%% Set variables
dt = t_axis(2) - t_axis(1);
dy = 50e-6;
c = 1500;

%% Filter
% butter function has (n, Wn, type) where n = nth order butterworth, Wn =
% normalized cutoff frequency and type is high or low

% might be interesting to see if you can write a function to calculate the
% best cutoff frequency??

[b,a] = butter(4,0.05,'high');
ts = filtfilt(b,a,t_data);
[b,a] = butter(4,0.5,'low');
ts = filtfilt(b,a,ts);

%% WHAT IS THIS?
nwin = 50;
if ~mod(nwin,2)
    nwin = nwin+1;
end 

%% Data Mapping
Y_new = zeros(size(ts));
Y_means = Y_new;
Y_mod = Y_new;

T2s = zeros(3, size(ts,2)); %rename this variable it sucks

for i = 1:size(ts,2)
    % indices with which to calculate average
    
    %     if i<=size(Y,2)-nwin+1;
    %         avg_ind_min=i;
    %         avg_ind_max=i+nwin-1;
    %     else
    %         avg_ind_min=size(Y,2)-nwin+1;
    %         avg_ind_max=size(Y,2);
    %     end
    avg_ind_min = max(1,i-(nwin-1)/2);
    avg_ind_max = min(size(ts,2),i+(nwin-1)/2);
    
    % calculate average
    Y_mean = mean(ts(:,avg_ind_min:avg_ind_max),2);
    
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

%% for recon??

forrecon = zeros(size(Y_new,1),size(Y_new,2)+400);
forrecon(:,201:size(Y_new,2)+200) = Y_new;

%% image construction
p_xy = kspaceLineRecon(forrecon,dy,dt/2,c);


%% FWHM
% this needs to be on the hilbert transformed data but before it is log
% transformed
FWHM_val = fwhm(x,y);