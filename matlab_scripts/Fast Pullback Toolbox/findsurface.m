function [depth,peakPos] = findsurface(data, threshold)

[data_rows, data_cols] = size(data);
peakPos = zeros(data_rows,data_cols);
depth = zeros(1,data_cols);

for col = 1:data_cols
    val = find(data(:,col)>threshold);
    if isempty(val)
        peakPos(1:end,col) = 0;
        %peakPos(end,col) = 1;
        depth(col) = data_cols-1;
    else
        depth(col) = val(1,1);
        peakPos(val(1:end,1),col) = 1;
    end
end

%depth = accumarray(Cx, Rx, size(data));
%depth(depth == 0) = 1;

%peakPos = 0;