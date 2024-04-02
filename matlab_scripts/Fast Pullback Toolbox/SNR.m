function value = SNR(ROI, noise)

% INPUTS:
%   ROI         - Section of img matrix
%
% OUTPUTS:
%   'res'       - [axial, lateral] FWHM

S = maxND(ROI);
value = abs(20*log10(S/noise));