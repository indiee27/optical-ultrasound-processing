% load the data

load('10mm');

c_s = 1500;
dx = 10/100; %[mm] i.e. 40 mm travelled with 100 Hz acquisition - might need adjusting
dx = dx/1000; %[m]
dt = 1e-8;

Yrecon = kspaceLineRecon(Y_new,dx,dt/2,c_s);
imagesc(xdim,ydim,20*log10(abs(hilbert(Yrecon))));
colormap('hot');
caxis([-30 -10]);

