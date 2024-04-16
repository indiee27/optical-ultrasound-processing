function [model_values] = FL_depth_correction(depth)
alpha = 12.71;  %in degrees
fdiameter = 0.4;    %fibre core diameter in mm
model_values = 1./(pi*(fdiameter+(depth*tan(deg2rad(alpha)))).^2);
%model_values = normalize(model_values, 'range');
