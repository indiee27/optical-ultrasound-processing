function filtered = butterfilt(data, varargin)

high = 0.05;
low = 0.5;

num_req_inputs = 1;
if nargin < num_req_inputs
    error('Not enough input variables')
elseif ~isempty(varargin)
    for input_index = 1:2:length(varargin)
        switch varargin{input_index}
            case 'high'
                high = varargin{input_index + 1};
            case 'highcutoff'
                high = varargin{input_index + 1};
            case 'low'
                low = varargin{input_index + 1};
            case 'lowcutoff'
                low = varargin{input_index + 1};
            otherwise 
                error('Unknown input')
        end
    end
end

[b,a] = butter(4, high, 'high');
ts = filtfilt(b, a, data);
[b,a] = butter(4, low, 'low');
filtered = filtfilt(b, a, ts);
