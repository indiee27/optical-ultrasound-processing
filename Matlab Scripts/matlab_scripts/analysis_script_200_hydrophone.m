%%-------------------------------------------------------------------------
%
% Ultrasound transmitter single point analysis script. This is edited from
% a script written by Dr. Erwin Alles and is the quick go to script for
% analysing data recorded using needle hydrophones when testing ultrasound
% transmitters.
%
%%-------------------------------------------------------------------------

%% Open the file and load the data from it

% This section simply loads the data file and assigns the data to vectors
% in matlab so that it is easy to manipulate.

[file,path] = uigetfile('*.txt');
fid = fopen([path,file]);
aa = fread(fid); %reading file
fclose(fid); %closing file
 
data = str2num(char(aa'));  %changes strings to numbers.
Nt = size(data,2);          %finding the length of time series (nt = size 
                                %of data), 2 is dimension i.e. 1st dimension is 
                                %length, 2nd is time.
clear aa; 
 
%window = 1:Nt;
window = 51:500;    %changes which bit of the time series you are looking at 
                        %"1" is the beginning and Nt is the end (length of 
                        %time series).
 
taxis = data(1,window);     %defining time axis. t axis = data so if data 
                                %is M x N matrix, taking first point of M 
                                %(entire column), then specify points i.e. 
                                %1 point on M and entire window on N
% (5,5:10) means point 5 on M and points 5-10 inclusive of 5 and 10 on N
 
dt = taxis(2)-taxis(1);     %time sep. dt is any two points
faxis = fftshift([0:(length(window)/2-1), -(length(window)/2):-1]/length(window)/dt);
%frequency axis is the Fourier transform of time axis

ch0 = data(2,window);   %now essentially reading the channel data out of 
                            %original data... ch is input on DAQ card
                            %ch0, ch1 are inputs on DAQ... 1. US info 
                            %2. pulse (note, these channels might be
                            %reversed depending on which PC you are using
ch1 = data(3,window);
 
%filter the low frequency noise
[b,a]=butter(4,0.3/50,'high');      %this creates a butterworth filter with
                                        %an order of 4, and a low frequency
                                        %cut-off of 0.1 MHz 
ch0=filtfilt(b,a,ch0); %this applies the filter created above
 
%% Set up the hydrophone calibration data and distance

%The data recorded on the hydrophone is in Volts and so must be converted
%to Pressure (MPa). We have a hydrophone calibration to do this, which
%gives us values of pressure per volt for discrete frequencies. This
%section interpolates over these discrete values to cover the measurement
%range. Then since it is for Frequency values, we Fourier transform the
%recorded data, apply the calibration, then reverse Fourier transform to
%get the calibrated time series.

% added 0 and Nyquist (50 MHz) frequencies
sen = [31 31 45 49 47 37 44 45 36 40 38 34 37 37 36 38 40 41 42 42 43 43 44 44 44 41 40 38 38 38 38 38]*10^(-3); %Sensitivity of hydrophone
freq = [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 50]; %corresponding freq
% These are taken from the factory calibration data of the needle
% hydrophone used. The hydrophones are typically calibrated at 1 MHz
% intervals from 1 to 30 MHz, so values are added at 0 MHz and 50 MHz (the
% Nyquist limit of our acquisition card to account for the uncalibrated
% region
 
N = size(ch0,2); 
fft_freqs=[0:(N-1)]/N*100;      %creates a frequency axis for the Fourier
                                    %transform in MHz
 
halfsen=interp1(freq,sen,fft_freqs(1:(N/2+1))); 
%halfsen = covers sensitivity for half of freq space. i.e. this covers the 
%positive frequency space since the Fourier Transform is symmetrical about
%zero. The interp1 function is used to interpolate from the sensitivity
%data above
fullsen=[halfsen fliplr(halfsen(2:end-1))];  %mirrors the positive sensitivity
                                                %to give the full range

fft_ch0 = fft(ch0);     %Fourier transform of the time series to get the
                            %data in frequency space.
scale_ch0 = fft_ch0./fullsen;   %This is the application of the calibration
                                    %to the data. i.e. it converts the
                                    %voltage values to pressure
 
ch2 = real(ifft(scale_ch0)); %ch2 is now scaled version of ch1. This is
                                %a reassignment of the data (pretty much
                                %unnecessary as I don't think we need the
                                %unscaled data after this point but oh
                                %well. We take the real part as the FFT is
                                %complex in nature.
 
CH2 = 20*log10(abs(fftshift(scale_ch0))); %power spectra

%This is a Fourier transform of the scaled data, the data is shifted using
%the fftshift to correct the phase. We take the absolute value as we only
%care about the magnitude of the frequency and apply 20*log10 as this
%converts from the pressure value to a power value (dB) which is the
%customary way to plot and quote this data.
 
 
%% In this section we plot the data

%what good is data unless we can see it! This is a simple and
%unsophisticated way of plotting the data. For publication we would tidy
%this up a lot. But this is better for quick displays of the data.

figure;
set(gcf,'position',[9    49   944   948]);
clear screen_size;
set(gcf, 'color', 'white');
%this opens a figure in MATLAB, sets the size and the figure color

%Here we create a subplot to display the time series.
subplot(2,1,1);
set(gca,'fontsize',14,'fontweight','demi'); %set the font size and weight
plot(taxis*1E6,[ch2],'linewidth',1);        %plot the data against time in us
xlabel('Time [\mus]');                      %label the time axis
ylabel('Pressure [MPa]');                   %label the pressure axis
title('Ultrasound time series');            %Give a title to the figure
pkpk = max(ch2) - min(ch2);                 %Find the peak-to-peak pressure
                                                %this essentially
                                                %subracts the minimum
                                                %ultrasound pressure from
                                                %the maximum pressure. It
                                                %is a value that is often
                                                %quoted to summarise the
                                                %ultrasound characteristics
set(gca,'xtick',[0:.5:5],'YMinorTick','off','XMinorTick','off','tickdir', 'out','box', 'off','fontsize',12,'fontname','arial','linewidth',1.5)%,'ticklength',tick_length)
%tidies up the labelling of the x axis

%Now we work on plotting the power spectrum

CH2 = CH2-max(CH2);     %here we normalise the power spectrum so that the
                            %value is zero. Again this is the customary way
                            %to display this as it makes it easier to see
                            %values that are used to quote the bandwidth,
                            %such as the -6 dB limit.

% select positive frequency space part of the spectrum
data6dB = CH2(end/2+1:end);

%upsample the spectrum for calculating the -6 dB bandwidth. This is another
%important value. It basically defines the frequencies at which the
%pressure falls to half the maximum value. We upsampled to increase the
%resolution of this bandwidth measurement
fband = faxis(end/2+1:end)/1E6;
upfband = linspace(0,125,1000);
updata6dB = interp1(fband,data6dB,upfband);

[maxdB, maxloc] = max(updata6dB);     %find the location of the maximum power

%find the bandwidth by finding the first point either side of the maximum
%to drop below the -6 dB limit
lowband = maxloc-find(fliplr(updata6dB(1:maxloc))<-6,1);    %find the low
                                                            %frequency
                                                            %limit
highband = find(updata6dB(maxloc:end)<-6,1)+maxloc;     %find the high
                                                        %frequency limit
bandwidth = upfband(highband) - upfband(lowband);   %find the bandwidth
                                                    %from the difference

%Here we create a subplot to display the power spectrum.
subplot(2,1,2);
set(gca,'fontsize',14,'fontweight','demi'); %set the font size and weight
plot(faxis/1E6,[CH2],'linewidth',1);        %plot the data against frequency in MHz
axis([0 50 -40 5]);                         %set the axis limits
xlabel('Frequency [MHz]');
ylabel('Power [dB]');
title('Power spectra');
%annotate the plot with the -6 dB bandwidth and the peak to peak pressure
%so that it is easy to see
text(5,-30,['Bandwidth = ' num2str(bandwidth) ' MHz']);
text(5,-35,['pk-pk Pressure = ' num2str(pkpk) ' MPa']);

