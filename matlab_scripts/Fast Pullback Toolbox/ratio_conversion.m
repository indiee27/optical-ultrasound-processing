function [x_axis, y_axis] = ratio_conversion(length, alines, depth, samples, data)

pixel_dimensions(1) = length/alines;
pixel_dimensions(2) = depth/samples;

image_depth = size(data,1)*pixel_dimensions(2);

image_width = size(data,2)*pixel_dimensions(1);

x_axis = linspace(0.0, image_width, size(data,2));

y_axis = linspace(0.0, image_depth, size(data,1));
