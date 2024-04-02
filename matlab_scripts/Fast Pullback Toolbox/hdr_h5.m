function metadata = hdr_h5(file)

hdr_info = h5info(file).Attributes;

old = {' ', '-'};
new = {'_'};
old2 = {'(', ')', '/'};
new2 = '';

metadata = struct;
for i = 1:length(hdr_info)
    fieldname = hdr_info(i).Name;
    fieldname = replace(fieldname, old, new);
    fieldname = replace(fieldname, old2, new2);
    metadata.(fieldname) = hdr_info(i).Value;
end

samples = metadata.Samples;
scan_depth = round(((samples*1500*1e-8)/2)*1000,0);
metadata.Scan_depth_mm = scan_depth;

step = (metadata.Scan_length_mm)/(metadata.A_lines);
metadata.Step = step;
