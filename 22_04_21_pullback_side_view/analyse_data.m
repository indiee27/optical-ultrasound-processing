%%-------------------------------------------------------------------------
%
% A script for analysing image data. This combines cross-talk removal and
% image reconstruction
%
%%-------------------------------------------------------------------------

%% Import the data and set the parameters
%assign files names
fileprefix = 'scan_16';
filehdr = [fileprefix '_hdr.txt'];
hdrinf = importdata(filehdr,'\t');
%Nx = hdrinf.data(2);
Nx = 200;
Nt = hdrinf.data(1);
max_file_number = 4;
number_of_images = (max_file_number+1)*2;
all_images = zeros(number_of_images,Nt,Nx);

dB = -20;
lowcutoff = 2;
highcutoff = 35;

for fnum = 0:max_file_number

    
    filehdr = [fileprefix '_hdr.txt'];
    filedata = [fileprefix '_image_0_f_' num2str(fnum) '.txt'];
    filetaxis = [fileprefix '_taxis.txt'];

    data = importdata(filedata);   %load the csv data from the .txt file
    data = data';                       %transpose the data to get the
                                            %orientation correct
    taxis = importdata(filetaxis)'; %import the time axis
    dt = taxis(2)-taxis(1);             %find the time step (this is needed for
                                            %the image reconstruction

    hdrinf = importdata(filehdr,'\t');    %load the header file to find the 
                                                %scan size and step size for
                                                %imagin
    dx = ((hdrinf.data(5))/(hdrinf.data(3)))*1e-3;       %x dimenstion step size [m]
    Nx = hdrinf.data(2);            %number of steps in the x dimension

    c = 1500;       %define the speed of sound [m/s]
    data = data';
    % we high pass filter the data to reduce the noise from the low frequencies
    % where the signal is weak. Further we can sacrifice the SNR by increasing
    % the low cut off, in doing so we can improve the resolution by removing
    % low frequency signal. There is a trade off.
    nyquistfreq = 50;       %nyquist frequency of the acquisition card [MHz]
    %lowcutoff = 5;          %low frequency cutoff - this is chosen empirically
                                %and can be changed to improve the image
                                %quality [MHz]
    [b,a]=butter(4,lowcutoff/nyquistfreq,'high');    %A butterworth filter is used. This filter
                                        %type can be changed and played with.
                                        %We have typically used butterworth as
                                        %it offers good performance
    data=filtfilt(b,a,data);        %this line applies the butterworth filter
                                        %to the data
    % we can also apply a low pass filter for similar reasons, this can be used
    % to improve the SNR by removing the low signal high frequencies.
    %highcutoff = 35;        %high frequency cutoff [MHz]
    [b,a]=butter(4,highcutoff/nyquistfreq,'low');
    data=filtfilt(b,a,data);

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

    nwin=25;        %sets the window size for the moving average
    if ~mod(nwin,2)     %makes window size odd if it is even for consistency
        nwin=nwin+1;
    end

    % create new matrices for the general linear model cross talk removal. Read
    % the paper to get more details on this
    Y_new=zeros(size(data));        %matrix for the new data
    Y_means=zeros(size(data));      %matrix for the means
    Y_mod=zeros(size(data));        %matrix of the modulus
    T2s=zeros(3,size(data,2));
    for i=1:size(data,2)

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
        avg_ind_max=min(size(data,2),i+(nwin-1)/2); %find the upper limit of 
                                                        %index for the window

        % calculate average - this accounts for variations in the amplitude of
        % the cross-talk from line to line
        Y_mean=mean(data(:,avg_ind_min:avg_ind_max),2);

        % derivative for temporal shift - this accounts for slight temporal
        % shifts in the cross-talk from line to line
        dY_mean=[diff(Y_mean);0];

        X2=[ones(size(data,1),1),Y_mean,dY_mean]; %create matrix to fit to the data

        % least-squares fit of the cross-talk to the data
        T2=pinv(X2'*X2)*X2'*data(:,i);

        % subtract fit - this removes the modelled cross-talk from the data,
        % leaving the new image data with the cross-talk removed
        Y_new(:,i)=data(:,i)-X2*T2;

        % store variables
        Y_means(:,i)=Y_mean;
        Y_mod(:,i)=X2*T2;
        T2s(:,i)=T2;

    end

    %% Reconstruct the image data

    %firstly we make a new matrix to put the data into before we reconstruct
    %it. The x-dimension is increased by 200 lines on either side of the data,
    %this gives a buffer of empty (zero signal a-lines) on either side of the
    %image. The reason we do this is that the reconstruction algorithm uses a
    %fast Fourier transform method. This assumes that the data is infinite,
    %i.e. that it repeats either side of the matrix. If we didn't have a buffer
    %of zeros, then when reconstructing, the right hand side of the image would
    %produce artefacts and features in the left hand side and vice versa.
    data_buf = zeros(size(Y_new,1),size(Y_new,2)+400);

    data_buf(:,201:size(Y_new,2)+200) = Y_new;  %insert the data into the empty
                                                    %matrix
    p_xy = kspaceLineRecon(data_buf,dx,dt/2,c); %this uses the kspaceLineRecon
                                                %algorithm from k-wave to
                                                %reconstruct the data. Look
                                                %this up and read the notes on
                                                %it. Figure out what it does
                                                %and how it works. Originally
                                                %it was designed for PA
                                                %imaging. This is why we half
                                                %the time step. Try to think
                                                %about why this is and why it
                                                %still works.
    data_rec = p_xy(:,201:end-200); %simply makes an image without the empty
                                    %buffer regions on the side

    %% Plot the final reconstructed image

    %data_rec(1:150,:) = data_rec(1:150,:)*0.05; %this is just an empirical
                                                %suppresion of any remaining
                                                %cross-talk. This can be done
                                                %in a more elegant way -> for
                                                %example by using a depth
                                                %verying gain for the signal.
                                                %this is something that you
                                                %could try implementing

    figure;
    %here we take the data and apply a hilbert transform. This gives us the
    %signal envelope for displaying the image. Then we apply a log transform
    %(log10) to improve the dynamic range of the plotted image. Ultrasound
    %images are typically plotted on a logarithmic scale.
    cur_im = 20*log10(abs(hilbert(squeeze(data_rec))));
    %We plot the image here, defining the x and z axes using the a-line step
    %spacing and the time axis
    imagesc(linspace(dx,dx*Nx,Nx)*1e3,c*1e3*taxis/2,cur_im);
    xlabel('x [mm]');       %label the x axis
    ylabel('depth [mm]');   %label the y axis (or what we think of as z here)
    CA = caxis;     %find the values of the colour axis
    %dB = -40;       %set the dynamic range for the image [dB]
    caxis([dB 0]+max(CA));  %define the new limits of the colour axis
    colormap(hot);          %set the colour map to hot (this is our preferred
                                %colour map
    set(gcf, 'Position',  [100, 100, 1000, 300])    %set the figure size
    axis tight equal    %make the axes tight to the image and to scale

    all_images((fnum*2)+1,:,:) = cur_im;
    %% Repeat for backward travelling images
    
    filehdr = [fileprefix '_hdr.txt'];
    filedata = [fileprefix '_image_0_b_' num2str(fnum) '.txt'];
    filetaxis = [fileprefix '_taxis.txt'];

    data = importdata(filedata);   %load the csv data from the .txt file
    data = data';                       %transpose the data to get the
                                            %orientation correct
    taxis = importdata(filetaxis)'; %import the time axis
    dt = taxis(2)-taxis(1);             %find the time step (this is needed for
                                            %the image reconstruction

    hdrinf = importdata(filehdr,'\t');    %load the header file to find the 
                                                %scan size and step size for
                                                %imagin
    dx = ((hdrinf.data(5))/(hdrinf.data(3)))*1e-3;;       %x dimenstion step size [m]
    Nx = hdrinf.data(2);            %number of steps in the x dimension

    c = 1500;       %define the speed of sound [m/s]
    data = data';
    % we high pass filter the data to reduce the noise from the low frequencies
    % where the signal is weak. Further we can sacrifice the SNR by increasing
    % the low cut off, in doing so we can improve the resolution by removing
    % low frequency signal. There is a trade off.
    nyquistfreq = 50;       %nyquist frequency of the acquisition card [MHz]
    %lowcutoff = 5;          %low frequency cutoff - this is chosen empirically
                                %and can be changed to improve the image
                                %quality [MHz]
    [b,a]=butter(4,lowcutoff/nyquistfreq,'high');    %A butterworth filter is used. This filter
                                        %type can be changed and played with.
                                        %We have typically used butterworth as
                                        %it offers good performance
    data=filtfilt(b,a,data);        %this line applies the butterworth filter
                                        %to the data
    % we can also apply a low pass filter for similar reasons, this can be used
    % to improve the SNR by removing the low signal high frequencies.
    %highcutoff = 35;        %high frequency cutoff [MHz]
    [b,a]=butter(4,highcutoff/nyquistfreq,'low');
    data=filtfilt(b,a,data);

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

    nwin=25;        %sets the window size for the moving average
    if ~mod(nwin,2)     %makes window size odd if it is even for consistency
        nwin=nwin+1;
    end

    % create new matrices for the general linear model cross talk removal. Read
    % the paper to get more details on this
    Y_new=zeros(size(data));        %matrix for the new data
    Y_means=zeros(size(data));      %matrix for the means
    Y_mod=zeros(size(data));        %matrix of the modulus
    T2s=zeros(3,size(data,2));
    for i=1:size(data,2)

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
        avg_ind_max=min(size(data,2),i+(nwin-1)/2); %find the upper limit of 
                                                        %index for the window

        % calculate average - this accounts for variations in the amplitude of
        % the cross-talk from line to line
        Y_mean=mean(data(:,avg_ind_min:avg_ind_max),2);

        % derivative for temporal shift - this accounts for slight temporal
        % shifts in the cross-talk from line to line
        dY_mean=[diff(Y_mean);0];

        X2=[ones(size(data,1),1),Y_mean,dY_mean]; %create matrix to fit to the data

        % least-squares fit of the cross-talk to the data
        T2=pinv(X2'*X2)*X2'*data(:,i);

        % subtract fit - this removes the modelled cross-talk from the data,
        % leaving the new image data with the cross-talk removed
        Y_new(:,i)=data(:,i)-X2*T2;

        % store variables
        Y_means(:,i)=Y_mean;
        Y_mod(:,i)=X2*T2;
        T2s(:,i)=T2;

    end

    %% Reconstruct the image data

    %firstly we make a new matrix to put the data into before we reconstruct
    %it. The x-dimension is increased by 200 lines on either side of the data,
    %this gives a buffer of empty (zero signal a-lines) on either side of the
    %image. The reason we do this is that the reconstruction algorithm uses a
    %fast Fourier transform method. This assumes that the data is infinite,
    %i.e. that it repeats either side of the matrix. If we didn't have a buffer
    %of zeros, then when reconstructing, the right hand side of the image would
    %produce artefacts and features in the left hand side and vice versa.
    data_buf = zeros(size(Y_new,1),size(Y_new,2)+400);

    data_buf(:,201:size(Y_new,2)+200) = Y_new;  %insert the data into the empty
                                                    %matrix
    p_xy = kspaceLineRecon(data_buf,dx,dt/2,c); %this uses the kspaceLineRecon
                                                %algorithm from k-wave to
                                                %reconstruct the data. Look
                                                %this up and read the notes on
                                                %it. Figure out what it does
                                                %and how it works. Originally
                                                %it was designed for PA
                                                %imaging. This is why we half
                                                %the time step. Try to think
                                                %about why this is and why it
                                                %still works.
    data_rec = p_xy(:,201:end-200); %simply makes an image without the empty
                                    %buffer regions on the side

    %% Plot the final reconstructed image

    %data_rec(1:150,:) = data_rec(1:150,:)*0.05; %this is just an empirical
                                                %suppresion of any remaining
                                                %cross-talk. This can be done
                                                %in a more elegant way -> for
                                                %example by using a depth
                                                %verying gain for the signal.
                                                %this is something that you
                                                %could try implementing

    figure;
    %here we take the data and apply a hilbert transform. This gives us the
    %signal envelope for displaying the image. Then we apply a log transform
    %(log10) to improve the dynamic range of the plotted image. Ultrasound
    %images are typically plotted on a logarithmic scale.
    cur_im = 20*log10(abs(hilbert(squeeze(data_rec))));
    %We plot the image here, defining the x and z axes using the a-line step
    %spacing and the time axis
    imagesc(linspace(dx,dx*Nx,Nx)*1e3,c*1e3*taxis/2,cur_im);
    xlabel('x [mm]');       %label the x axis
    ylabel('depth [mm]');   %label the y axis (or what we think of as z here)
    CA = caxis;     %find the values of the colour axis
    %dB = -40;       %set the dynamic range for the image [dB]
    caxis([dB 0]+max(CA));  %define the new limits of the colour axis
    colormap(hot);          %set the colour map to hot (this is our preferred
                                %colour map
    set(gcf, 'Position',  [100, 100, 1000, 300])    %set the figure size
    axis tight equal    %make the axes tight to the image and to scale
    
    all_images((fnum+1)*2,:,:) = cur_im;
end

savename = [fileprefix 'image'];
save(savename,'all_images');

%% Make a video of the results

writerObj = VideoWriter(savename);     %define a video object in matlab
writerObj.FrameRate = 2;               %set the video frame rate

figure;         %open a figure to put the video frames in

open(writerObj);    %open the video object to write data to it

cmax = maxND(all_images);
dBrange = -1*dB;

% now we loop through the data, showing 100 a-lines at a time. The image is
% cycled to the left, adding a new a-line on the right hand side with each
% frame in the video.
for fnum = 1:size(all_images,1)
    imagesc(linspace(dx,dx*Nx,Nx)*1e3,c*1e3*taxis/2,squeeze(all_images(fnum,:,:)));
    caxis([cmax-dBrange cmax]);
    colormap hot
    xlabel('Width (mm)');
    ylabel('Depth (mm)');
    axis tight equal
    F = getframe(gcf);
    writeVideo(writerObj, F);
    
    pause(0.1);
end

close(writerObj);   %close the video object now that all the frames are written