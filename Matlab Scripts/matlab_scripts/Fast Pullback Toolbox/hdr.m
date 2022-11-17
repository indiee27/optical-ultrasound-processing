function info = hdr(file_num)

try 
    hdr_info = importdata(['scan_' num2str(file_num) '_hdr.txt'], '\t');
catch
    hdr_info = importdata(['scan_' num2str(file_num) '_hdr.txt'], '\t');
end

if length(hdr_info.data) == 5
    samples = hdr_info.data(1);
    a_lines = hdr_info.data(2);
    rep_rate = hdr_info.data(3);
    scan_length = hdr_info.data(4);
    pullback_speed = num2str(hdr_info.data(5));
    pullback_speed = ([pullback_speed ' mm/s']);
    
    % calcs
    step = scan_length/a_lines; 
    
elseif length(hdr_info.data) == 2
    a_lines = hdr_info.data(1);
    step = hdr_info.data(2);
    
    % calcs
    scan_length = a_lines * step;
    pullback_speed = 'B-mode';
    
    f1 = importdata(['scan_' num2str(file_num) '_data.txt']);
    samples = size(f1,2);
    clear f1
else
    error('Incorrect hdr file format')
end

scan_depth = round(((samples*1500*1e-8)/2)*1000,0);

% add to struct
info = struct;
info.samples = samples;
info.a_lines = a_lines;
info.step = step;
info.pullback_speed = pullback_speed;
info.scan_length = scan_length;
info.scan_depth = scan_depth;
