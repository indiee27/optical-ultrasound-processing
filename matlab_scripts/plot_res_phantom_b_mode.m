%% Script to plot the resolution of reconstructed images
% - tungsten wire phantom
% - b mode

%mkdir res_figs

q = figure;
set(gcf, 'color', 'white');

% params
c = 1500;
dB = -35;
nwin = 40;

filt_high = 0.05;
filt_low = 0.5;

list = ['1' '2' '4' '8'];
points = 10;

%% Plot fig
f_path = 'scan_3';

hdr_info = importdata([f_path '_hdr.txt']);
Nx = hdr_info.data(1);
Ny = hdr_info.data(2);
dx = hdr_info.data(3);
dy = hdr_info.data(4);

% load files
f1 = importdata([f_path '_1.txt']);
samples = length(f1);
clear f1

tdata = zeros(Ny, samples);
for m = 1:Ny
    fid = fopen([t '_' num2str(m) '.txt']);
    aa = fread(fid);
    fclose(fid);
       
    tdata(m,:) = str2double(char(aa'));
    clear aa;
end
tdata = tdata'; 

t_axis = importdata([f_path '_taxis.txt']);
dt = (taxis(2) - taxis(1));

% sizes
scan_length = Ny*dy;
scan_depth = round(((samples*c*1e-8)/2)*1000,0);
