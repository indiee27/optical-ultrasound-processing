function res = PointRes(ROI, varargin)

% INPUTS:
%   ROI         - Section of img matrix
%
% OUTPUTS:
%   'res'       - [axial, lateral] FWHM

signal_axial = max(ROI);
signal_lateral = max(ROI);

[val, peak_loc] = max(signal_lateral);
signal_norm_lateral = signal_lateral - val;
    
[val2, peak_loc2] = max(signal_axial);
signal_norm_axial = signal_axial - val2;

% find the lateral FWHM
dB_point_pos1 = find(signal_norm_lateral(peak_loc:end) <=-6, 1, 'first') - 1;
dB_point_neg1 = find(fliplr(signal_norm_lateral(1:peak_loc)) <=-6, 1, 'first') - 1;

% find the axial FWHM
dB_point_pos2 = find(signal_norm_axial(peak_loc2:end) <= -6, 1, 'first') - 1;
dB_point_neg2 = find(fliplr(signal_norm_axial(1:peak_loc2)) <= -6, 1, 'first') -1 ;
FWHM_lateral = dB_point_pos1 + dB_point_neg1;
FWHM_axial = dB_point_pos2 + dB_point_neg2;

res = [FWHM_axial FWHM_lateral];

