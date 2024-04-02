function [data, ratio] = H5toUS(fname)
dt = 1e-8;
%dt = h5readatt(fname, '/', 'dt');
dx = h5readatt(fname, '/', 'dx')*1e-3;
%samples = h5readatt(fname, '/', 'samples');
%Nx = h5readatt(fname, '/', 'Nx');
%length = h5readatt(fname, '/', 'length_mm');
%depth = h5readatt(fname, '/', 'depth_mm');

US = h5read(fname, '/US');
US = butterfilt(US);
US = crosstalk(US,40);
US = TGC(US);
USmat = zeros(size(US,1), size(US,2)+400);
USmat(:,201:size(US,2)+200) = US;
USmat = kspaceLineRecon(USmat, dx, dt/2, 1500);

data = 20*log10(abs(hilbert(USmat(:,201:end-200))));
f = figure('Visible','off');
imagesc(data);
CA = clim;
close(f);
data(data<(-35+max(CA))) = -35+max(CA);
data(data>(0+max(CA))) = 0+max(CA);

%ratio_y = round((length/Nx)*(1/(depth/samples)));
%ratio = [1 ratio_y 1];

data = uint8(255 * mat2gray(data));