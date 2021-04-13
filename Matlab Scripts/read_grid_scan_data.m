function [data,taxis,faxis,pos,Nt,Npos] = read_grid_scan_data(filepath,filebase)
% 
% function [data,taxis,faxis,pos,Nt,Npos] = read_grid_scan_data(filepath,filebase)
%
% Input parameters:
% filepath  - string pointing to the directory containing the grid scan data
% filebase  - the filename of the main data file. Template: NAME.bin
% 
% Output variables:
% data      - the matrix (Npos x Nt x 2) containing the measurement data.
%             contains data for both channel 0 (:,:,1) and 1 (:,:,2) [V]
% taxis     - vector containing the time corresponding to each sample [s]
% faxis     - vector containing the positive frequency components [Hz]
% pos       - Matrix (Npos x 3) containing the (x,y,z) coordinates of each
%             measurent location [m]
% Nt        - Number of time samples
% Npos      - Number of measurement locations


if ~strcmp(filepath(end),'\')
    filepath = [filepath,'\'];
end

if ~isempty(strfind(filebase,'.bin'))
    filebase = filebase(1:end-4);
end




filename = [filebase,'.bin'];
axesname = [filebase,'_axes.bin'];
possname = [filebase,'_positions.bin'];

s = dir([filepath,axesname]);
try
    Nt = s.bytes/16;
catch
    disp(' ');
    disp(['Warning: file ',[filepath,axesname],' not found.']);
    data = [];taxis = [];faxis = [];pos = [];Nt = [];Npos = [];
    return;
end
s = dir([filepath,filename]);
Npos = s.bytes/16/Nt;
clear s;

fid = fopen([filepath,axesname]);
axes = fread(fid,2*Nt,'double','b');
fclose(fid);
taxis = axes(1:Nt);
faxis = axes(Nt+1:floor(Nt*1.5))*1E6;
clear fid ans axes axesname;

fid = fopen([filepath,filename]);
data_read = fread(fid,2*Nt*Npos,'double','b');
fclose(fid);
data = reshape(data_read,Nt,2,Npos);
clear fid ans data_read filename;
data = permute(data,[3,1,2]);

fid = fopen([filepath,possname]);
pos_read = fread(fid,3*Npos,'double','b');
fclose(fid);
pos = reshape(pos_read,3,Npos)';
clear pos_read;

clear ans fid *name filepath filename;