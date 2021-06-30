%% Script for FWHM of resolution wire phantom
% side viewing fibre 2021-05-21

%% List all the fnames
% NB I've left out the 30mm scans
scan_list = {'speed10ms_reprate100hz_length20mm', 'speed10ms_reprate200hz_length20mm', 'speed15ms_reprate150hz_length20mm', 'speed15ms_reprate300hz_length20mm', 'speed20ms_reprate200hz_length20mm', 'speed20ms_reprate400hz_length20mm', 'speed25ms_reprate250hz_length20mm', 'speed25ms_reprate500hz_length20mm', 'speed30ms_reprate300hz_length20mm', 'speed30ms_reprate600hz_length20mm', 'speed35ms_reprate350hz_length20mm', 'speed35ms_reprate700hz_length20mm', 'speed40ms_reprate400hz_length20mm', 'speed40ms_reprate800hz_length20mm', 'speed45ms_reprate450hz_length20mm', 'speed45ms_reprate900hz_length20mm', 'speed50ms_reprate500hz_length20mm', 'speed50ms_reprate1000hz_length20mm', 'speed60ms_reprate600hz_length20mm', 'speed70ms_reprate700hz_length20mm', 'speed80ms_reprate800hz_length20mm', 'speed90ms_reprate900hz_length20mm', 'speed100ms_reprate1000hz_length20mm'};

%% Loop through all the images
% matrix initialisation
lateral_resolution = zeros(2,2);
axial_resolution = zeros(2,2);
speed_reprate = zeros(2,2);

for n = 1:2%length(scan_list)
    % filename assign
    f_name = scan_list{n};
    
   
    % import data
    t_data = importdata([f_name '_data.txt']);
    t_axis = importdata([f_name '_taxis.txt']);
    %t_axis = linspace(10e-9,26.67e-6,2667);
    hdr_inf = importdata([f_name '_hdr.txt'], '\t');
    speed_reprate(n, 1) = hdr_inf.data(5);
    speed_reprate(n, 2) = hdr_inf.data(3);
    
    % set constants
    dt = t_axis(2) - t_axis(1);
    dy = 100e-6;
    c = 1500;
    dx = hdr_inf.data(5)/hdr_inf.data(3);
    
    %{
if pullback_speed/rep_rate = dx
    continue
else
    n = n+1
end
%}
    % high pass filter
    [b,a] = butter(4,0.02,'high');
    ts = filtfilt(b,a,t_data);
    [b,a] = butter(4,0.7,'low');
    ts = filtfilt(b,a,ts);
    
    % window
    nwin = 5;%20-60
    if ~mod(nwin,2)
        nwin = nwin+1;
    end 
    
    % matrix initialising
    Y_new = zeros(size(ts));
    Y_means = zeros(size(ts));
    Y_mod = zeros(size(ts));
    T2s = zeros(3,size(ts,2));
    
    %% cross-talk removal
    for i = 1:size(ts,2)
    
        % indices with which to calculate average
        % if i<=size(Y,2)-nwin+1;
        %   avg_ind_min=i;
        %   avg_ind_max=i+nwin-1;
        % else
        %   avg_ind_min=size(Y,2)-nwin+1;
        %   avg_ind_max=size(Y,2);
        % end
        
        avg_ind_min = max(1, i-(nwin-1)/2);
        avg_ind_max = min(size(ts,2), i+(nwin-1)/2);
    
        % calculate average
        Y_mean = mean(ts(:,avg_ind_min:avg_ind_max),2);
    
        % derivative for temporal shift
        dY_mean = [diff(Y_mean);0];
        %dY2_mean = [diff(dY_mean);0];
    
        %X2 = [ones(size(Y,1),1), Y_mean,dY_mean,dY2_mean];
        X2 = [ones(size(ts,1),1), Y_mean,dY_mean];
        %X2 = [ones(size(Y,1),1), Y_mean];
    
        % least-squares fit
        T2 = pinv(X2'*X2)*X2'*ts(:,i);
    
        % subtract fit
        Y_new(:,i) = ts(:,i) - X2*T2;
    
        % store variables
        Y_means(:,i) = Y_mean;
        Y_mod(:,i) = X2*T2;
        T2s(:,i) = T2;
    
    end
    
    %% image reconstruction
    forrecon = zeros(size(Y_new,1),size(Y_new,2)+400);
    forrecon(:,201:size(Y_new,2)+200) = Y_new;

    p_xy = kspaceLineRecon(forrecon,dy,dt/2,c);
    
    img = abs(hilbert(p_xy(:,201:end-200)));
    % THIS LINE WAS MISSING
    
    fig1 = figure(1);
    imagesc(img)
    set(gcf);
    
    %% Loop through the five wires in each image
        
    % results matrix
    lateral_results = [];
    axial_results = [];
    
    % wire list
    wires = {wire_1, wire_2, wire_3, wire_4, wire_5};
    
    wire_1 = imrect;
    wire_2 = imrect;
    wire_3 = imrect;
    wire_4 = imrect;
    wire_5 = imrect;
    
    for x = 1:2%length(wires)
        %wire = wires{x};
        
        pos_raw = getPosition(wire_1);
        pos_raw = round(pos_raw);
        
        % data
        raw_data = abs(img);
        
        % data for wire section
        ROI_raw = raw_data(pos_raw(2):(pos_raw(2)+pos_raw(4)),pos_raw(1):(pos_raw(1)+pos_raw(3)));
        
        %% max signal
        signal_axial = max(ROI_raw);
        signal_lateral = max(ROI_raw,[],2);
        
        % normalize central lobe
        [val, peak_loc] = max(signal_lateral);	% location of the mainlobe
        signal_norm_lateral = signal_lateral - val;     % normalized signal
        
        [val2, peak_loc2] = max(signal_axial);
        signal_norm_axial = signal_axial - val2;
        
        %% plot signal
        fig2 = figure;
        
        % laterial
        subplot(1,2,1);
        plot(signal_norm_lateral,'Color',[0,0.7,0.9]);
        title('Lateral signal');
        xlabel('Location');
        ylabel('Amplitude');
        
        % axial
        subplot(1,2,2);
        plot(signal_norm_axial,'Color',[0,0.7,0.9]);
        title('Axial signal');
        xlabel('location');
        ylabel('Amplitude');
        
        % add to results matrix
        lateral_results(:,x) = signal_norm_lateral;
        axial_results(:,x) = signal_norm_axial;
        
        % close 
        msgbox('Save the fucking graph Indie');
                      
    end 
    
    %% average of 5 wires
    lateral_average = sum(lateral_results,1)/5;
    axial_average = sum(axial_results,1)/5;
    
%     % find peak
%     [val, peak_loc] = max(lateral_average);
%     [val2, peak_loc2] = max(axial_average);
%     
%     % find lateral resolution
%     dB_point_pos1 = find(lateral_average(peak_loc:end) <=-6, 1, 'first') - 1;
%     dB_point_neg1 = find(fliplr(lateral_average(1:peak_loc)) <=-6, 1, 'first') - 1;
%     
%     % find the axial resolution
%     dB_point_pos2 = find(axial_average(peak_loc2:end) <= -6, 1, 'first') - 1;
%     dB_point_neg2 = find(fliplr(axial_average(1:peak_loc2)) <= -6, 1, 'first') -1 ;
% 
%     %  Linear interpolation to find the exact lateral FWHM
%     dB_point_pos_interp1 = dB_point_pos1 + (((-6) - lateral_average(peak_loc + dB_point_pos1)) / (signal_norm_lateral(peak_loc + dB_point_pos1) - signal_norm_lateral(peak_loc + dB_point_pos1 - 1)));
%     dB_point_neg_interp1 = dB_point_neg1 + (((-6) - lateral_average(peak_loc - dB_point_neg1)) / (signal_norm_lateral(peak_loc - dB_point_neg1) - signal_norm_lateral(peak_loc - dB_point_neg1 + 1)));
% 
%     if (numel(dB_point_pos_interp1) == 0 || numel(dB_point_neg1) == 0)
%         FWHM_lateral = [];
%     else
%         FWHM_lateral = dB_point_pos_interp1 + dB_point_neg_interp1;
%     end
% 
%     %  Linear interpolation to find the exact axial FWHM
%     dB_point_pos_interp2 = dB_point_pos2 + (((-6) - signal_norm_axial(peak_loc2 + dB_point_pos2)) / (signal_norm_axial(peak_loc2+dB_point_pos2) - signal_norm_axial(peak_loc2+dB_point_pos2-1)));
%     dB_point_neg_interp2 = dB_point_neg2 + (((-6) - signal_norm_axial(peak_loc2 - dB_point_neg2)) / (signal_norm_axial(peak_loc2-dB_point_neg2) - signal_norm_axial(peak_loc2-dB_point_neg2+1)));
% 
%     if (numel(dB_point_pos_interp2) == 0 || numel(dB_point_neg2) == 0)
%         FWHM_axial = [];
%     else
%         FWHM_axial = dB_point_pos_interp2 + dB_point_neg_interp2;
%     end
%     
%     close gcf;
%   
    FWHM_lateral = fwhm((lateral_average), dx);
    FWHM_axial = fwhm((axial_average), dx);
    
    % add to results
    lateral_resolution(:,n) = FWHM_lateral;
    axial_resolution(:,n) = FWHM_axial;
    
end