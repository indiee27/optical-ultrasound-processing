function filter = generate_filter(taxis,fmin,fmax,NH)
% filter = generate_filter(taxis,fmin,fmax,NH)
% 
% NH = number of points in hann filter - EVEN NUMBER!

NH = NH - mod(NH,2);

H = hann(2*NH);
H = H(NH+1:end);

Nt = length(taxis);
dt = taxis(2)-taxis(1);
faxis = 1/dt/Nt * [0:ceil(Nt/2)-1,-floor(Nt/2):-1];

filter = ones(1,Nt);

cut_off_low  = find(abs(faxis) >= fmin, 1) + floor(NH/2)-1;
cut_off_high = find(faxis-fmax >= 0,    1) - floor(NH/2);
filter(1:cut_off_low-1) = 0;
filter(cut_off_high+1:end) = 0;

if isempty(cut_off_high);cut_off_high=ceil(Nt/2)-1; end

% High-pass (> fmin):
cnt = 0;
while cut_off_low-cnt >= 1 && cnt < NH
    filter(cut_off_low-cnt) = H(cnt+1);
    cnt = cnt+1;
end

% Low-pass (<fmax):
cnt = 0;
while cut_off_high+cnt < ceil(Nt/2)-1 && cnt < NH
    filter(cut_off_high+cnt) = H(cnt+1);
    cnt = cnt+1;
end

% Copy to negative frequencies:
filter(end:-1:end-floor(Nt/2)+2) = filter(2:floor(Nt/2));

% %%
% figure;plot(faxis/1E6,filter);
% axis([fmin/1E6-5 fmax/1E6+5 0 1]);