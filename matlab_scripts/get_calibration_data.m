function cal_data_int = get_calibration_data(faxis,probe)


if strcmp(probe,'75')           % SN2320
    cal_freq = (1:30)*1E6;
    cal_data = [5 6 8 10 10 8 8 8 8 9 12 13 12 11 11 10 10 10 10 10 10 10 10 9 9 8 8 8 9 10]*1E-3;              % V/MPa for 75 um needle hydrophone S/N23201
    CORR     = [0.021288 0.37849 0.56944 0.66922 0.73899 0.67209 0.76585 1.1408 1.1382 1.3663 ...
                1.6015 1.8295 2.7685 1.4455 1.3617 1.3753 1.4511 1.5126 1.4853 1.5161 1.5107 ...
                1.5101 1.5672 1.5844 1.5012 1.3717 1.2998 1.2656 1.4026 1.6339];                                % Richard worked this out on 18 Nov 2016
    cal_data = cal_data./CORR;
elseif strcmp(probe,'200')      % SN2393
    cal_freq = (1:30)*1E6;
    cal_data = [53 77 89 88 82 81 82 73 68 62 62 70 67 64 65 68 73 71 70 69 70 71 72 73 73 71 70 68 66 65]*1E-3;% V/MPa for 200 um needle hydrophone S/N2393
elseif strcmp(probe,'2672')     % 200 um NEW
    cal_freq = (1:30)*1E6;
    cal_data = [31 45 49 47 37 44 45 36 40 38 34 37 37 36 38 40 41 42 42 43 43 44 44 44 41 40 38 38 38 38]*1E-3;% V/MPa for 200 um needle hydrophone S/N2672
elseif strcmp(probe,'2698')     % 75 um NEW
    cal_freq = (1:30)*1E6;
    cal_data = [5 6 8 10 12 13 12 11 10 10 10 9 8 7 7 8 8 9 9 9 9 9 9 9 9 9 8 8 9 9]*1E-3;                      % V/MPa for 75 um needle hydrophone S/N2698
end

cal_data_int = interp1(cal_freq,cal_data,abs(faxis),'pchip');   
cal_data_int(abs(faxis)>max(cal_freq)) = cal_data(end);
cal_data_int(abs(faxis)<min(cal_freq)) = cal_data(1);
