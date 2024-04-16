function [values] = depth_correction(depth, step_size)

theta = rad2deg(tan(max(depth)/length(depth)*step_size));
%B = 90 - theta;
fD = 0.4;
NA = 12.71;
a = fD/2 + tan(deg2rad(NA));
values = zeros(1,length(depth));

for n = 1:length(depth)
    %d = max(depth) - depth(n);
    b = depth(n)/sin(deg2rad(theta)) * sin(deg2rad(NA));
    area = a * (b) * pi;
    I = 1/area;
    values(n) = I;
end
