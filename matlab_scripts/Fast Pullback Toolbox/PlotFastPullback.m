function [img] = PlotFastPullback(file_num, x, varargin)

% INPUTS:
%   file_name       - string of filepath
%   x               - reconstruction spacing - every x a-lines
%
% OPTIONAL INPUTS:
%   Optional 'string' value pairs that may be used to modify the default
%   settings
%   'Interp'        - interpolation method
%   'TGC'           - Boolean controlling whether time gain compensation is
%                     applied
%   'xTalk'         - Boolean controlling whether the cross talk algorithm is applied
%   'Interp'        - string input controlling the interpolation method
%                     used by interp2 in the reconstruction (default = '*nearest').
%   'nwin'          - nwindow
%   'Plot raw data' - plot raw data as well as reconstructed
%   'SoundSpeed'    - speed of sound in medium
%   'dB'            - dB
%
% OUTPUTS:
%   img             - matrix result

% define defaults
num_req_inputs = 2;
TGC = false;
xTalk = true;
nwin = 40;
plot_raw = false;
filt_high = 0.05;
filt_low = 0.5;
c = 1500;
dB = -35;
interp_method = '*nearest';

if nargin < num_req_inputs
    error('Not enough input variables');
elseif ~isempty(varargin)
    for input_index = 1:2:length(varargin)
        switch varargin{input_index}
            case 'Interp'
                interp_method = varargin{input_index + 1};
            case 'TGC'
                TGC = varargin{input_index + 1};
            case 'xTalk'
                xTalk = varargin{input_index + 1};
            case 'nwin'
                nwin = varargin{input_index + 1};
            case 'PlotRaw'
                plot_raw = varargin{input_index + 1};
            case 'FiltHigh'
                filt_high = varargin{input_index + 1};
            case 'FiltLow'
                filt_low = varargin{input_index + 1};
            case 'SoundSpeed'
                c = varargin{input_index + 1};
            case 'dB'
                dB = varargin{input_index + 1};
            otherwise
                error('Unknown optional input');
        end
    end
end

% create saving folder
if ~isfolder('images')
    mkdir images
end
file_num = num2str(file_num);

% import the three files
tdata = importdata(['scan_' file_num '_data.txt']);
taxis = importdata(['scan_' file_num '_taxis.txt']); dt = (taxis(2) - taxis(1));
hdr_info = importdata(['scan_' file_num '_hdr.txt'], '\t');

% header info
samples = hdr_info.data(1);
scan_length = hdr_info.data(4);
pullback_speed = hdr_info.data(5);
a_lines = hdr_info.data(2)/x;

% sizes
step = scan_length/a_lines; dy = step*1e-3;
scan_depth = round(((samples*1500*1e-8)/2)*1000,0);
ratio = round((scan_length/a_lines)*(1/(scan_depth/samples)));

% data acquisition
tseries = zeros(samples, a_lines);
for g = 1:(a_lines - 1)
    h = 1 + (g - 1)*x;
    tseries(:,g) = tdata(:,h);
end
clear g;
clear h;

[b, a] = butter(4, filt_high, 'high');
ts = filtfilt(b, a, tseries);
[b, a] = butter(4, filt_low, 'low');
ts = filtfilt(b, a, ts);

if ~mod(nwin,2)
    nwin = nwin+1;
end

% CROSSTALK
if xTalk
    Y_new = zeros(size(ts));
    Y_means = Y_new; Y_mod = Y_new;
    T2s = zeros(4, size(ts,2));
    for i = 1:size(ts,2)
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
end

% TIME GAIN COMPENSATION
if TGC 
    i_max = round(0.1 * samples,0);
    lamda = 2.5;
    for l = 1:a_lines
        for k = 1:samples
            gain_factor = (min(k, i_max)/ i_max)^lamda;
            Y_new(k,l) = Y_new(k,l)*gain_factor;
        end
    end
end
    
% RECONSTRUCTION
forrecon = zeros(size(Y_new,1), size(Y_new,2)+400);
forrecon(:,201:size(Y_new,2)+200) = Y_new;

p_xy = kspaceLineRecon(forrecon, dy*2, dt, c, 'Interp', interp_method);
img = 20*log10(abs(hilbert(p_xy)));

if plot_raw == false
    imagesc(img(:,201:end-200));
    CA = caxis;
    caxis([dB 0]+max(CA));
    xticks([1 a_lines/2 a_lines-1])
    xticklabels({0 scan_length/2 scan_length})
    xlabel('Length [mm]')
    yticks([1 samples/2 samples-1])
    yticklabels({0 scan_depth/2 scan_depth})
    ylabel('Depth [mm]')
    title(['Scan ' num2str(x) ': speed: ' num2str(pullback_speed) ' mm/s, spacing: ' num2str(dy*1e6) ' \mu m']);
    subtitle('Reconstructed data (log+hilbert)');
    colormap(hot);
    colorbar
    set(gca, 'DataAspectRatio', [1 ratio 1])
    figure_name = (['/images/' file_num ' spacing ' num2str(step*1e3) '.png']);
    saveas(gcf,[pwd figure_name]);
    
elseif plot_raw == true
    figure(1);
    subplot(2,1,1);
    imagesc(20*log10(abs(hilbert(forrecon))));
    CA = caxis;
    caxis([dB 0]+max(CA));
    subtitle('Raw data (log+hilbert)');
    xticks([200 200+a_lines/2 200+a_lines])
    xticklabels({0 scan_length/2 scan_length})
    xlabel('Length [mm]')
    yticks([1 samples/2 samples-1])
    yticklabels({0 scan_depth/2 scan_depth})
    ylabel('Depth [mm]')
    colorbar
    subplot(2,1,2);
    imagesc(img(:,201:end-200));
    CA = caxis;
    caxis([dB 0]+max(CA));
    xticks([1 a_lines/2 a_lines-1])
    xticklabels({0 scan_length/2 scan_length})
    xlabel('Length [mm]')
    yticks([1 samples/2 samples-1])
    yticklabels({0 scan_depth/2 scan_depth})
    ylabel('Depth [mm]')
    title(['Scan ' num2str(x) ': speed: ' num2str(pullback_speed) ' mm/s, spacing: ' num2str(dy*1e6) ' \mu m']);
    subtitle('Reconstructed data (log+hilbert)');
    colormap(hot);
    colorbar
    set(gca, 'DataAspectRatio', [1 ratio 1])
    figure_name = (['/images/' file_num ' spacing ' num2str(step*1e3) '.png']);
    saveas(gcf,[pwd figure_name]);
else
    error('Cannot plot')
end
    

