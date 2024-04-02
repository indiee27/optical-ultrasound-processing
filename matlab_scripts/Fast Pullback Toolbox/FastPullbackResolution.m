function results = FastPullbackResolution(file_num, x, pts, varargin)

% INPUTS:
%   file_name       - scan number
%   x               - reconstruction spacing - every x a-lines
%   pts             - number of assessment points
%
% OPTIONAL INPUTS:
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
% OUPUTS:
%   'SNR'           - matrix of SNR vs depth
%   'lat'           - matrix of lateral vs depth
%   'ax'            - matrix of axial vs depth

% default values
num_req_inputs = 3;
TGC = false;
xTalk = true;
nwin = 40;
filt_high = 0.05;
filt_low = 0.5;
c = 1500;
dB = -35;
interp_method = '*nearest';
dt = 1e-8;
dy = 25e-6;

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
            case 'FiltHigh'
                filt_high = varargin{input_index + 1};
            case 'FiltLow'
                filt_low = varargin{input_index + 1};
            case 'SoundSpeed'
                c = varargin{input_index + 1};
            case 'dB'
                dB = varargin{input_index + 1};
            case 'dt'
                dt = varargin{input_index + 1};
            case 'dy'
                dy = varargin{input_index} + 1;
            otherwise
                error('Unknown optional input');
        end
    end
end

img = PlotFastPullback(file_num, x, ...
    'Interp', interp_method, 'TGC', TGC, 'xTalk', xTalk, 'nwin', nwin, 'FiltHigh', filt_high, 'FiltLow', filt_low, 'SoundSpeed', c, 'dB', dB);

snr = zeros(2,pts);
lat = snr; ax = snr;

rawdata = (img.^10)/20;
noise_floor = rawdata(100:200, 200:300);
noise = mean(mean(noise_floor));

for n = 1:pts
    hFig = figure;
    set(gcf, 'Position', get(0, 'Screensize'));
    imagesc(img)
    title(['Depth ' num2str(n)]);
    colormap hot
    
    % find position
    h2 = drawrectangle;
    pos_recon = h2.Position;
    pos_recon = round(pos_recon);
    
    rawdata_recon = rawdata(:,201:end-200);
    ROI = rawdata_recon(pos_recon(2):(pos_recon(2)+pos_recon(4)),pos_recon(1):(pos_recon(1)+pos_recon(3)));
    
    signal_axial = max(ROI');
    signal_lateral = max(ROI', [], 2);
    
    [val, peak_loc] = max(signal_lateral);
    signal_norm_lateral = signal_lateral - val;
    
    [val2, peak_loc2] = max(signal_axial);
    signal_norm_axial = signal_axial - val2;
    
    res_depth = (peak_loc + pos_recon(2)) * ((dt * 1e6 * c/2)/1000);

    ROI_SNR = rawdata(pos_recon(2):(pos_recon(2)+pos_recon(4)),pos_recon(1):(pos_recon(1)+pos_recon(3)));
    S = maxND(ROI_SNR);
    SNR = abs(20*log10(S/noise));

    % find the lateral FWHM
    dB_point_pos1 = find(signal_norm_lateral(peak_loc:end) <=-6, 1, 'first') - 1;
    dB_point_neg1 = find(fliplr(signal_norm_lateral(1:peak_loc)) <=-6, 1, 'first') - 1;
    
    % find the axial FWHM
    dB_point_pos2 = find(signal_norm_axial(peak_loc2:end) <= -6, 1, 'first') - 1;
    dB_point_neg2 = find(fliplr(signal_norm_axial(1:peak_loc2)) <= -6, 1, 'first') -1 ;
    FWHM_lateral = dB_point_pos1 + dB_point_neg1;
    FWHM_axial = dB_point_pos2 + dB_point_neg2;

    R_axial = 1e6*FWHM_axial*(dt*c/2);
    R_lateral = 1e6*FWHM_lateral*dy*x;
    
    % add to results
    snr(1,n) = res_depth; snr(2,n) = SNR;
    lat(1,n) = res_depth; lat(2,n) = R_lateral;
    ax(1,n) = res_depth; ax(2,n) = R_axial;
    
    close(hFig)
end

results = [snr; lat; ax];