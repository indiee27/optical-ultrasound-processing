%% Function to create data of the right size
function tseries = SpacingCompress(tdata, step)
tseries = zeros(size(tdata,1),((size(tdata,2)/step)-1));
for g = 1:size(tseries,2)
    h = 1 + (g - 1)*step;
    tseries(:,g) = tdata(:,h);
end