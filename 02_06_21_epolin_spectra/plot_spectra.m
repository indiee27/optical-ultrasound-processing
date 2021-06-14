%% LOAD FLAME DATA

f_ref_light = (importdata('flame_light.txt')); f_ref_dark = (importdata('flame_dark.txt'));
f_x = f_ref_light(1,:);
f_y = f_ref_light(2,:) - f_ref_dark(2,:);

% EPOLIN 2MG

f2a_light = (importdata('flame_2_1_light.txt')); f2a_dark = (importdata('flame_2_1_dark.txt'));
f2a = ((f2a_light(2,:))-(f2a_dark(2,:)));
f2b_light = (importdata('flame_2_2_light.txt')); f2b_dark = (importdata('flame_2_2_dark.txt'));
f2b = ((f2b_light(2,:))-(f2b_dark(2,:)));
f2c_light = (importdata('flame_2_3_light.txt')); f2c_dark = (importdata('flame_2_3_dark.txt'));
f2c = ((f2c_light(2,:))-(f2c_dark(2,:)));
f2d_light = (importdata('flame_2_4_light.txt')); f2d_dark = (importdata('flame_2_4_dark.txt'));
f2d = ((f2d_light(2,:))-(f2d_dark(2,:)));
f2e_light = (importdata('flame_2_5_light.txt')); f2e_dark = (importdata('flame_2_5_dark.txt'));
f2e = ((f2e_light(2,:))-(f2e_dark(2,:)));

figure;

title Flame 2mg
plot(f_x, f2a)
hold on
plot(f_x, f2b)
hold on
plot(f_x, f2c)
hold on
plot(f_x, f2d)
hold on
plot(f_x, f2e)
hold on
plot(f_x, f_y)
legend('2_1', '2_2', '2_3', '2_4', '2_5', 'Ref')

% EPOLIN 5MG

f5a_light = (importdata('flame_5_1_light.txt')); f5a_dark = (importdata('flame_5_1_dark.txt'));
f5a = ((f5a_light(2,:))-(f5a_dark(2,:)));
f5b_light = (importdata('flame_5_2_light.txt')); f5b_dark = (importdata('flame_5_2_dark.txt'));
f5b = ((f5b_light(2,:))-(f5b_dark(2,:)));
f5c_light = (importdata('flame_5_3_light.txt')); f5c_dark = (importdata('flame_5_3_dark.txt'));
f5c = ((f5c_light(2,:))-(f5c_dark(2,:)));
f5d_light = (importdata('flame_5_4_light.txt')); f5d_dark = (importdata('flame_5_4_dark.txt'));
f5d = ((f5d_light(2,:))-(f5d_dark(2,:)));
f5e_light = (importdata('flame_5_5_light.txt')); f5e_dark = (importdata('flame_5_5_dark.txt'));
f5e = ((f5e_light(2,:))-(f5e_dark(2,:)));

figure;
title Flame 5mg
plot(f_x, f5a)
hold on
plot(f_x, f5b)
hold on
plot(f_x, f5c)
hold on
plot(f_x, f5d)
hold on
plot(f_x, f5e)
hold on
plot(f_x, f_y)
legend('5_1', '5_2', '5_3', '5_4', '5_5', 'Ref')

%% LOAD NIR DATA

n_ref_light = (importdata('nir_light.txt')); n_ref_dark = (importdata('nir_dark.txt'));
n_x = n_ref_light(1,:);
n_y = n_ref_light(2,:) - n_ref_dark(2,:);

% EPOLIN 2MG

n2a_light = (importdata('nir_2_1_light.txt')); n2a_dark = (importdata('nir_2_1_dark.txt'));
n2a = ((n2a_light(2,:))-(n2a_dark(2,:)));
n2b_light = (importdata('nir_2_2_light.txt')); n2b_dark = (importdata('nir_2_2_dark.txt'));
n2b = ((n2b_light(2,:))-(n2b_dark(2,:)));
n2c_light = (importdata('nir_2_3_light.txt')); n2c_dark = (importdata('nir_2_3_dark.txt'));
n2c = ((n2c_light(2,:))-(n2c_dark(2,:)));
n2d_light = (importdata('nir_2_4_light.txt')); n2d_dark = (importdata('nir_2_4_dark.txt'));
n2d = ((n2d_light(2,:))-(n2d_dark(2,:)));
n2e_light = (importdata('nir_2_5_light.txt')); n2e_dark = (importdata('nir_2_5_dark.txt'));
n2e = ((n2e_light(2,:))-(n2e_dark(2,:)));

figure;
title NIR 2mg
plot(n_x, n2a)
hold on
plot(n_x, n2b)
hold on
plot(n_x, n2c)
hold on
plot(n_x, n2d)
hold on
plot(n_x, n2e)
hold on
plot(n_x, n_y)
legend('2_1', '2_2', '2_3', '2_4', '2_5', 'Ref')

% EPOLIN 5MG

n5a_light = (importdata('nir_5_1_light.txt')); n5a_dark = (importdata('nir_5_1_dark.txt'));
n5a = ((n5a_light(2,:))-(n5a_dark(2,:)));
n5b_light = (importdata('nir_5_2_light.txt')); n5b_dark = (importdata('nir_5_2_dark.txt'));
n5b = ((n5b_light(2,:))-(n5b_dark(2,:)));
n5c_light = (importdata('nir_5_3_light.txt')); n5c_dark = (importdata('nir_5_3_dark.txt'));
n5c = ((n5c_light(2,:))-(n5c_dark(2,:)));
n5d_light = (importdata('nir_5_4_light.txt')); n5d_dark = (importdata('nir_5_4_dark.txt'));
n5d = ((n5d_light(2,:))-(n5d_dark(2,:)));
n5e_light = (importdata('nir_5_5_light.txt')); n5e_dark = (importdata('nir_5_5_dark.txt'));
n5e = ((n5e_light(2,:))-(n5e_dark(2,:)));

figure;
title NIR 5mg
plot(n_x, n5a)
hold on
plot(n_x, n5b)
hold on
plot(n_x, n5c)
hold on
plot(n_x, n5d)
hold on
plot(n_x, n5e)
hold on
plot(n_x, n_y)
legend('5_1', '5_2', '5_3', '5_4', '5_5', 'Ref')

%% PERCENTAGES
% Replace negative values from the flame results
f2a(f2a < 0) = 0;
f2b(f2b < 0) = 0;
f2c(f2c < 0) = 0;
f2d(f2d < 0) = 0;
f2e(f2e < 0) = 0;
f5a(f5a < 0) = 0;
f5b(f5b < 0) = 0;
f5c(f5c < 0) = 0;
f5d(f5d < 0) = 0;
f5e(f5e < 0) = 0;
f_y(f_y < 0) = 0;

% Percentage calculations
f2a_y = f2a./f_y;
f2a_y(f2a_y > 1) = 1;
f2b_y = f2b./f_y;
f2b_y(f2b_y > 1) = 1;
f2c_y = f2c./f_y;
f2c_y(f2c_y > 1) = 1;
f2d_y = f2d./f_y;
f2d_y(f2d_y > 1) = 1;
f2e_y = f2e./f_y;
f2e_y(f2e_y > 1) = 1;
f5a_y = f5a./f_y;
f5a_y(f5a_y > 1) = 1;
f5b_y = f5b./f_y;
f5b_y(f5b_y > 1) = 1;
f5c_y = f5c./f_y;
f5c_y(f5c_y > 1) = 1;
f5d_y = f5d./f_y;
f5d_y(f5d_y > 1) = 1;
f5e_y = f5e./f_y;
f5e_y(f5e_y > 1) = 1;
n2a_y = n2a./n_y;
n2b_y = n2b./n_y;
n2c_y = n2c./n_y;
n2d_y = n2d./n_y;
n2e_y = n2e./n_y;
n5a_y = n5a./n_y;
n5b_y = n5b./n_y;
n5c_y = n5c./n_y;
n5d_y = n5d./n_y;
n5e_y = n5e./n_y;

figure
title Flame 2mg
xlim([300 1100])
ylim([0 1])
plot(f_x, f2a_y)
hold on
plot(f_x, f2b_y)
hold on
plot(f_x, f2c_y)
hold on
plot(f_x, f2d_y)
hold on
plot(f_x, f2e_y)
hold off
legend('2_1', '2_2', '2_3', '2_4', '2_5')

figure
title Flame 5mg
xlim([300 1100])
ylim([0 1])
plot(f_x, f5a_y)
hold on
plot(f_x, f5b_y)
hold on
plot(f_x, f5c_y)
hold on
plot(f_x, f5d_y)
hold on
plot(f_x, f5e_y)
hold off
legend('5_1', '5_2', '5_3', '5_4', '5_5')

figure;
title NIR 2mg

plot(n_x, n2a_y)
hold on
plot(n_x, n2b_y)
hold on
plot(n_x, n2c_y)
hold on
plot(n_x, n2d_y)
hold on
plot(n_x, n2e_y)
hold off
legend('2_1', '2_2', '2_3', '2_4', '2_5')

figure;
title NIR 5mg

plot(n_x, n5a_y)
hold on
plot(n_x, n5b_y)
hold on
plot(n_x, n5c_y)
hold on
plot(n_x, n5d_y)
hold on
plot(n_x, n5e_y)
hold off
legend('5_1', '5_2', '5_3', '5_4', '5_5')

