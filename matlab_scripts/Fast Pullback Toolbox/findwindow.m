function window = findwindow(data)

%ref = min(data');

%for x = 1:size(data,2)
%[~,loc] = find(data(:,x) > ref*1.2);
%end

%window = [min(loc) max(loc)];
window = zeros(2,size(data,2));
for x = 1:size(data,2)
    [loc,~] = find(data(:,x) > 0.5);
    window(1,x) = loc(1);
    window(2,x) = loc(end);
    [~,peak] = max(data(:,x));
    window(3,x) = peak;
end