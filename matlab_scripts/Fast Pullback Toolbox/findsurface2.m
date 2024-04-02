function [depth, peakPos] = findsurface2(data, threshold, visibility)

[data_rows, data_cols] = size(data);
peakPos = zeros(data_rows, data_cols);

for col = 1:data_cols
    peakFound = false;
    for row = threshold + 1:data_rows
        if data(row, col) > threshold
            peakPos(row:row + visibility, col) = 1;
            peakFound = true;
        elseif peakFound
            break;
        end
    end
end

depth = zeros(1, size(peakPos, 2));
for col = 1:size(peakPos, 2)
    row_index = find(peakPos(:, col), 1, 'first');
    if ~isempty(row_index)
        depth(col) = row_index;
    end
end