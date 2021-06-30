%% this is for the research in art competition
t_data = importdata('speed50ms_reprate1000hz_length30mm_data.txt');

t_axis = linspace(10e-9,26.67e-6,2667);
%% its a kind of floral shape so try reflecting it and reflecting again?
% change color gradient throughout?
%% 3D image??

dt = t_axis(2) - t_axis(1);
dy = 100e-6;
c = 1500;

%ts = t_data;
[b,a]=butter(4,0.05,'high');
ts=filtfilt(b,a,t_data);
[b,a]=butter(4,0.7,'low');
ts=filtfilt(b,a,ts);

nwin=10;
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

img = (abs(hilbert(p_xy(:,201:end-200))));
dB = -25;
%figure;
%imagesc(img)
I2 = flip(img ,2);           %# horizontal flip
I3 = flip(img ,1);           %# vertical flip
I4 = flip(I3,2);    %# horizontal+vertical flip
% figure;
% imagesc(img)
% set(gca, 'visible', 'off')
% figure;
% imagesc(I2)
% set(gca, 'visible', 'off')
% figure;
% imagesc(I3)
% set(gca, 'visible', 'off')
% figure;
% imagesc(I4)
% set(gca, 'visible', 'off')
%CA = caxis;
%caxis([dB 0]+max(CA));
%colormap copper
%shading interp

added_img = [I3 I4; img I2];
imagesc(added_img)
set(gca, 'visible', 'off')
colormap((pink))
shading interp