function [surface, depth] = findsurfacefinal(US,varargin)

[rows,cols] = size(US);
surface = zeros(rows,cols);
depth = zeros(1,cols); depth(1:end) = rows;

band_size = 5;
threshold = 85;

num_req_inputs = 1;
if nargin < num_req_inputs
    error('Not enough input variables')
elseif ~isempty(varargin)
    for input_index = 1:2:length(varargin)
        switch varargin{input_index}
            case 'band_size'
                band_size = varargin{input_index + 1};
            case 'threshold'
                threshold = varargin{input_index + 1};
            otherwise 
                error('Unknown input')
        end
    end
end

for col = 1:cols
    column = US(:,col);
    area = find(column(:,1) >= threshold);

    if isempty(area)
        surface(1:end,col) = 0;
        if col <= band_size
            depth(col) = mean(depth(1:band_size));
        else
            depth(col) = mean(depth(col-band_size:col));
        end

    elseif size(area,1) < band_size
        surface(1:end,col) = 0;
        if col <= band_size
            depth(col) = mean(depth(1:band_size));
        else
            depth(col) = mean(depth(col-band_size:col));
        end
        
    else
        area_end = find(column(area(1,1):end) < threshold,1);
        if area_end >= band_size
            surface(area(1:area_end-1,1),col) = 1;
            depth(col) = area(1,1);

        elseif area_end < band_size
            x = 0;
            while (area_end < band_size) && (x < size(area,1)-band_size)
                x = x + 1;
                area_end = find(column(area(x,1):end) < threshold ,1);
            end

            if area_end-x < band_size
                surface(1:end,col) = 0;
                if col <= band_size
                    depth(col) = mean(depth(1:band_size));
                else
                    depth(col) = mean(depth(col-band_size:col));
                end

            elseif area_end-x >= band_size
                surface(area(x:area_end,1),col) = 1;
                depth(col) = area(x,1);
            end
        end
    end
end
