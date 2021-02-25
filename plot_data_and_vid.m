%%-------------------------------------------------------------------------
%
% This script is made to load the image data, reshape it, then process it
% and plot it as an image as well as make a video of the data which is a
% reconstruction of how it was acquired.
%
%%-------------------------------------------------------------------------
% would be good to turn this whole script into a function to do this for
% every distance
% maybe work out how to make some ROC curves from there?
% where is the video?

%% First load and reshape the data

probe_loc = '20 mm per second 20 mm travel';     %this defines the location of the US probe when the
                    %measurement was taken. i.e. the different data folders
                    %labelled X1 - 9
filename= ['.\' probe_loc '\depth_scan_00000.dat'];   %define the filename using the 
                                            %probe location
headername = ['.\' probe_loc '\depth_scan_00000.txt'];%define the header file using the 
                                            %probe location

hdrinf = importfile(headername);    %uses a matlab autogenerated function 
                                    %to load the header file info. These
                                    %matlab autogenerated functions are a
                                    %great way to be lazy when you need a
                                    %script like this
                                    
samp_per_scan=hdrinf{6,2};     %this is the numbers of samples per each a-line,
                                %i.e. the number of time points

% finds the number of a-lines that were acquired during the image
% acquisition. We need to know this so that we can reshape the data
fid=fopen(filename,'r');

Y0 = fread(fid,'double');   %reads the data from the file into a matrix
%then we reshape the data using the number of samples per scan. So that it
%is arrange as a-lines acquired side by side.

Y0 = reshape(Y0,samp_per_scan,size(Y0,1)/samp_per_scan);

%% Apply some frequency filters to the data

%here we apply frequency filters to the data. These are used to improve the
%signal to noise and the resolution of the image. You can have a play with
%the cutoffs to see what the effect is.

%define the cutoffs and nyquist frequency
nyquistfreq = hdrinf{3,2}/(2e3); %the nyquist frequency is half the sample
                                    %rate. Note that the sample rate is
                                    %recorded in kHz so we convert to MHz
                                    %for the nyquist
lowcutoff = 3e3;    %define the cut off frequencies in MHz
highcutoff = 40e3;

%Now we apply the filters to the data (this is the same as in the other
%imaging scripts)
% high pass filter
[b,a]=butter(4,lowcutoff/nyquistfreq,'high');
Y0=filtfilt(b,a,Y0);


% low pass filter
[b,a]=butter(4,highcutoff/nyquistfreq,'low');
Y0=filtfilt(b,a,Y0);

%% Apply the cross-talk removal algorithm
% this cross-talk removal is described in detail in the 2016 BOE vascular
% tissue imaging paper. Essentially, when we are imaging we have ultrasound
% that travels directly from the transmitter to the receiver. This signal
% is quite strong, but also quite uniform over time. In order to remove it
% and clean up the image we model the signal using the previous A-lines and
% then subtract it from the current image line. By using this type of
% moving average, the signal from tissue, which is not stationary over time
% is not removed. One obvious flaw with this method is that if there is a
% real signal that does not vary over time then it may be removed too.

nwin=50;        %sets the window size for the moving average
if ~mod(nwin,2)     %makes window size odd if it is even for consistency
    nwin=nwin+1;
end

% create new matrices for the general linear model cross talk removal. Read
% the paper to get more details on this
Y_new=zeros(size(Y0));        %matrix for the new Y0
Y_means=zeros(size(Y0));      %matrix for the means
Y_mod=zeros(size(Y0));        %matrix of the modulus
T2s=zeros(3,size(Y0,2));
for i=1:size(Y0,2)

    avg_ind_min=max(1,i-(nwin-1)/2);            %find the lower limit of
                                                    %index for the window
                                                    %for the moving
                                                    %average. We need this
                                                    %since the moving
                                                    %average is based on
                                                    %the previous a-lines
                                                    %and for the first
                                                    %samples, up to the
                                                    %size of the window,
                                                    %there will be less
                                                    %samples before the
                                                    %a-lines than the
                                                    %number defined by the
                                                    %window
    avg_ind_max=min(size(Y0,2),i+(nwin-1)/2); %find the upper limit of 
                                                    %index for the window

    % calculate average - this accounts for variations in the amplitude of
    % the cross-talk from line to line
    Y_mean=mean(Y0(:,avg_ind_min:avg_ind_max),2);

    % derivative for temporal shift - this accounts for slight temporal
    % shifts in the cross-talk from line to line
    dY_mean=[diff(Y_mean);0];

    X2=[ones(size(Y0,1),1),Y_mean,dY_mean]; %create matrix to fit to the data

    % least-squares fit of the cross-talk to the data
    T2=pinv(X2'*X2)*X2'*Y0(:,i);

    % subtract fit - this removes the modelled cross-talk from the data,
    % leaving the new image data with the cross-talk removed
    Y_new(:,i)=Y0(:,i)-X2*T2;

    % store variables
    Y_means(:,i)=Y_mean;
    Y_mod(:,i)=X2*T2;
    T2s(:,i)=T2;

end


%% Apply signal suppression to the beginning of the image
%here some depth dependent signal suppression is applied to the image to
%remove any residual cross-talk and to tidy up the final image

% make the x and y dimensions for the image
c_s = 1500;     %define the speed of sound
rep_rate = 100;
samp_rate = hdrinf{3,2};
%use the sample rate and the speed of sounds to find the total image depth
ydim = linspace(0,(size(Y_new,1)*(1/samp_rate))*c_s/2,size(Y_new,1));
xdim = linspace(0,size(Y_new,2)/rep_rate,size(Y_new,2));

%take the hilbert transform to get the signal envelope and then apply a log
%transform to display the image
lY_new = 20*log10(abs(hilbert(Y_new)));

%apply signal suppression to the first few time points. This removes any
%residual cross-talk. The numbers were chosen empirically and follow a
%similar method to that used in the manuscript in this folder.
nummax = 400;
gamma = 0.3;
susp = ones(size(lY_new,1),1);

%creates a vector to apply as signal suppresion to each a-line. Has a depth
%dependent magnitude
for num = 1:nummax
    susp(num) = (num/nummax)^gamma;
end
%apply the suppresion to the data
slY_new = lY_new./susp;

%% Plot a figure with the US image

% plot figure with set size
ax1 = figure('pos',[10,50,1000,500]);
imagesc(xdim,ydim,slY_new);
colormap hot

% define the dB color scale
cmax = maxND(slY_new);
dBrange = 30;
caxis([cmax-dBrange cmax]);

% label the plot
xlabel('Time (s)');
ylabel('Depth (mm)');

print(ax1,probe_loc,'-dsvg');         %this line is used to save the
% figure as an svg file (image)


return

%% Make a video of the data

writerObj = VideoWriter(probe_loc);     %define a video object in matlab
writerObj.FrameRate = rep_rate;               %set the video frame rate

figure;         %open a figure to put the video frames in

open(writerObj);    %open the video object to write data to it

% now we loop through the data, showing 100 a-lines at a time. The image is
% cycled to the left, adding a new a-line on the right hand side with each
% frame in the video.
for fnum = 1:(size(Y_new,2)-101)
    imagesc(xdim(fnum:fnum+100),ydim,slY_new(:,fnum:fnum+100));
    caxis([cmax-dBrange cmax]);
    colormap hot
    xlabel('Time (s)');
    ylabel('Depth (mm)');
    F = getframe(gcf);
    writeVideo(writerObj, F);
    
    pause(0.1);
end

close(writerObj);   %close the video object now that all the frames are written