%% Epolin 1
% reference values
fopen('/02_06_21_epolin_spectra/');
f_ref_light_a = importdata('flame_light.txt'); f_ref_dark_a = (importdata('flame_dark.txt'));
f_x_a = f_ref_light_a(1,:);
f_y_a = f_ref_light_a(2,:) - f_ref_dark_a(2,:);
nir_ref_light_a = importdata('nir_light.txt'); nir_ref_dark_a = importdata('nir_dark.txt');
n_x_a = nir_ref_light_a(1,:);
n_y_a = nir_ref_light_a(2,:) - nir_ref_dark_a(2,:);

%2mg
%f2a_light_a = (importdata('flame_2_1_light.txt')); f2a_dark = (importdata('flame_2_1_dark.txt'));
%f2a_a = ((f2a_light_a(2,:))-(f2a_dark(2,:)));
f2b_light_a = (importdata('flame_2_2_light.txt')); f2b_dark_a = (importdata('flame_2_2_dark.txt'));
f2b_a = ((f2b_light_a(2,:))-(f2b_dark_a(2,:)));
f2c_light_a = (importdata('flame_2_3_light.txt')); f2c_dark_a = (importdata('flame_2_3_dark.txt'));
f2c_a = ((f2c_light_a(2,:))-(f2c_dark_a(2,:)));
f2d_light_a = (importdata('flame_2_4_light.txt')); f2d_dark_a = (importdata('flame_2_4_dark.txt'));
f2d_a = ((f2d_light_a(2,:))-(f2d_dark_a(2,:)));
f2e_light_a = (importdata('flame_2_5_light.txt')); f2e_dark_a = (importdata('flame_2_5_dark.txt'));
f2e_a = ((f2e_light_a(2,:))-(f2e_dark_a(2,:)));

f5a_light_a = (importdata('flame_5_1_light.txt')); f5a_dark_a = (importdata('flame_5_1_dark.txt'));
f5a_a = ((f5a_light(2,:))-(f5a_dark_a(2,:)));
f5b_light_a = (importdata('flame_5_2_light.txt')); f5b_dark_a = (importdata('flame_5_2_dark.txt'));
f5b_a = ((f5b_light_a(2,:))-(f5b_dark_a(2,:)));
f5c_light_a = (importdata('flame_5_3_light.txt')); f5c_dark_a = (importdata('flame_5_3_dark.txt'));
f5c_a = ((f5c_light_a(2,:))-(f5c_dark_a(2,:)));
f5d_light_a = (importdata('flame_5_4_light.txt')); f5d_dark_a = (importdata('flame_5_4_dark.txt'));
f5d_a = ((f5d_light_a(2,:))-(f5d_dark_a(2,:)));
f5e_light_a = (importdata('flame_5_5_light.txt')); f5e_dark_a = (importdata('flame_5_5_dark.txt'));
f5e_a = ((f5e_light_a(2,:))-(f5e_dark_a(2,:)));

%n2a_light_a = (importdata('nir_2_1_light.txt')); n2a_dark_a = (importdata('nir_2_1_dark.txt'));
%n2a_a = ((n2a_light_a(2,:))-(n2a_dark_a(2,:)));
n2b_light_a = (importdata('nir_2_2_light.txt')); n2b_dark_a = (importdata('nir_2_2_dark.txt'));
n2b_a = ((n2b_light_a(2,:))-(n2b_dark_a(2,:)));
n2c_light_a = (importdata('nir_2_3_light.txt')); n2c_dark_a = (importdata('nir_2_3_dark.txt'));
n2c_a = ((n2c_light_a(2,:))-(n2c_dark_a(2,:)));
n2d_light_a = (importdata('nir_2_4_light.txt')); n2d_dark_a = (importdata('nir_2_4_dark.txt'));
n2d_a = ((n2d_light_a(2,:))-(n2d_dark_a(2,:)));
n2e_light_a = (importdata('nir_2_5_light.txt')); n2e_dark_a = (importdata('nir_2_5_dark.txt'));
n2e_a = ((n2e_light_a(2,:))-(n2e_dark_a(2,:)));

n5a_light_a = (importdata('nir_5_1_light.txt')); n5a_dark_a = (importdata('nir_5_1_dark.txt'));
n5a_a = ((n5a_light_a(2,:))-(n5a_dark_a(2,:)));
n5b_light_a = (importdata('nir_5_2_light.txt')); n5b_dark_a = (importdata('nir_5_2_dark.txt'));
n5b_a = ((n5b_light_a(2,:))-(n5b_dark_a(2,:)));
n5c_light_a = (importdata('nir_5_3_light.txt')); n5c_dark_a = (importdata('nir_5_3_dark.txt'));
n5c_a = ((n5c_light_a(2,:))-(n5c_dark_a(2,:)));
n5d_light_a = (importdata('nir_5_4_light.txt')); n5d_dark_a = (importdata('nir_5_4_dark.txt'));
n5d_a = ((n5d_light_a(2,:))-(n5d_dark_a(2,:)));
n5e_light_a = (importdata('nir_5_5_light.txt')); n5e_dark_a = (importdata('nir_5_5_dark.txt'));
n5e_a = ((n5e_light_a(2,:))-(n5e_dark_a(2,:)));

% Percentages

%f2a_a_y = f2a_a./f_y_a;
f2b_a_y = f2b_a./f_y_a;
f2c_a_y = f2c_a./f_y_a;
f2d_a_y = f2d_a./f_y_a;
f2e_a_y = f2e_a./f_y_a;
f5a_a_y = f5a_a./f_y_a;
f5b_a_y = f5b_a./f_y_a;
f5c_a_y = f5c_a./f_y_a;
f5d_a_y = f5d_a./f_y_a;
f5e_a_y = f5e_a./f_y_a;
%n2a_a_y = n2a_a./n_y_a;
n2b_a_y = n2b_a./n_y_a;
n2c_a_y = n2c_a./n_y_a;
n2d_a_y = n2d_a./n_y_a;
n2e_a_y = n2e_a./n_y_a;
n5a_a_y = n5a_a./n_y_a;
n5b_a_y = n5b_a./n_y_a;
n5c_a_y = n5c_a./n_y_a;
n5d_a_y = n5d_a./n_y_a;
n5e_a_y = n5e_a./n_y_a;

fclose('all');

%% Epolin 2
fopen('/10_06_21_epolin_spectra/');

f_ref_light_b = importdata('flame_light.txt'); f_ref_dark_b = (importdata('flame_dark.txt'));
f_x_b = f_ref_light_b(1,:);
f_y_b = f_ref_light_b(2,:) - f_ref_dark_b(2,:);
nir_ref_light_b = importdata('nir_light.txt'); nir_ref_dark_b = importdata('nir_dark.txt');
n_x_b = nir_ref_light_b(1,:);
n_y_b = nir_ref_light_b(2,:) - nir_ref_dark_b(2,:);

%2mg
f2a_light = (importdata('flame_2_1_light.txt')); f2a_dark = (importdata('flame_2_1_dark.txt'));
f2a_b = ((f2a_light(2,:))-(f2a_dark(2,:)));
f2b_light = (importdata('flame_2_2_light.txt')); f2b_dark = (importdata('flame_2_2_dark.txt'));
f2b_b = ((f2b_light(2,:))-(f2b_dark(2,:)));
f2c_light = (importdata('flame_2_3_light.txt')); f2c_dark = (importdata('flame_2_3_dark.txt'));
f2c_b = ((f2c_light(2,:))-(f2c_dark(2,:)));
f2d_light = (importdata('flame_2_4_light.txt')); f2d_dark = (importdata('flame_2_4_dark.txt'));
f2d_b = ((f2d_light(2,:))-(f2d_dark(2,:)));
f2e_light = (importdata('flame_2_5_light.txt')); f2e_dark = (importdata('flame_2_5_dark.txt'));
f2e_b = ((f2e_light(2,:))-(f2e_dark(2,:)));

f5a_light = (importdata('flame_5_1_light.txt')); f5a_dark = (importdata('flame_5_1_dark.txt'));
f5a_b = ((f5a_light(2,:))-(f5a_dark(2,:)));
f5b_light = (importdata('flame_5_2_light.txt')); f5b_dark = (importdata('flame_5_2_dark.txt'));
f5b_b = ((f5b_light(2,:))-(f5b_dark(2,:)));
f5c_light = (importdata('flame_5_3_light.txt')); f5c_dark = (importdata('flame_5_3_dark.txt'));
f5c_b = ((f5c_light(2,:))-(f5c_dark(2,:)));
f5d_light = (importdata('flame_5_4_light.txt')); f5d_dark = (importdata('flame_5_4_dark.txt'));
f5d_b = ((f5d_light(2,:))-(f5d_dark(2,:)));
f5e_light = (importdata('flame_5_5_light.txt')); f5e_dark = (importdata('flame_5_5_dark.txt'));
f5e_b = ((f5e_light(2,:))-(f5e_dark(2,:)));

n2a_light = (importdata('nir_2_1_light.txt')); n2a_dark = (importdata('nir_2_1_dark.txt'));
n2a_b = ((n2a_light(2,:))-(n2a_dark(2,:)));
n2b_light = (importdata('nir_2_2_light.txt')); n2b_dark = (importdata('nir_2_2_dark.txt'));
n2b_b = ((n2b_light(2,:))-(n2b_dark(2,:)));
n2c_light = (importdata('nir_2_3_light.txt')); n2c_dark = (importdata('nir_2_3_dark.txt'));
n2c_b = ((n2c_light(2,:))-(n2c_dark(2,:)));
n2d_light = (importdata('nir_2_4_light.txt')); n2d_dark = (importdata('nir_2_4_dark.txt'));
n2d_b = ((n2d_light(2,:))-(n2d_dark(2,:)));
n2e_light = (importdata('nir_2_5_light.txt')); n2e_dark = (importdata('nir_2_5_dark.txt'));
n2e_b = ((n2e_light(2,:))-(n2e_dark(2,:)));

n5a_light = (importdata('nir_5_1_light.txt')); n5a_dark = (importdata('nir_5_1_dark.txt'));
n5a_b = ((n5a_light(2,:))-(n5a_dark(2,:)));
n5b_light = (importdata('nir_5_2_light.txt')); n5b_dark = (importdata('nir_5_2_dark.txt'));
n5b_b = ((n5b_light(2,:))-(n5b_dark(2,:)));
n5c_light = (importdata('nir_5_3_light.txt')); n5c_dark = (importdata('nir_5_3_dark.txt'));
n5c_b = ((n5c_light(2,:))-(n5c_dark(2,:)));
n5d_light = (importdata('nir_5_4_light.txt')); n5d_dark = (importdata('nir_5_4_dark.txt'));
n5d_b = ((n5d_light(2,:))-(n5d_dark(2,:)));
n5e_light = (importdata('nir_5_5_light.txt')); n5e_dark = (importdata('nir_5_5_dark.txt'));
n5e_b = ((n5e_light(2,:))-(n5e_dark(2,:)));

% Percentages

f2a_b_y = f2a_b./f_y_b;
f2b_b_y = f2b_b./f_y_b;
f2c_b_y = f2c_b./f_y_b;
f2d_b_y = f2d_b./f_y_b;
f2e_b_y = f2e_b./f_y_b;
f5a_b_y = f5a_b./f_y_b;
f5b_b_y = f5a_b./f_y_b;
f5c_b_y = f5a_b./f_y_b;
f5d_b_y = f5a_b./f_y_b;
f5e_b_y = f5a_b./f_y_b;
n2a_b_y = n2a_b./n_y_b;
n2b_b_y = n2b_b./n_y_b;
n2c_b_y = n2c_b./n_y_b;
n2d_b_y = n2d_b./n_y_b;
n2e_b_y = n2e_b./n_y_b;
n5a_b_y = n5a_b./n_y_b;
n5b_b_y = n5b_b./n_y_b;
n5c_b_y = n5c_b./n_y_b;
n5d_b_y = n5d_b./n_y_b;
n5e_b_y = n5e_b./n_y_b;

%% Variations

% flame
%f2a = abs(f2a_a_y - f2a_b_y);
f2b = abs(f2b_a_y - f2b_b_y);
f2c = abs(f2c_a_y - f2c_b_y);
f2d = abs(f2d_a_y - f2d_b_y);
f2e = abs(f2e_a_y - f2e_b_y);

f5a = abs(f5a_a_y - f5a_b_y);
f5b = abs(f5b_a_y - f5b_b_y);
f5c = abs(f5c_a_y - f5c_b_y);
f5d = abs(f5d_a_y - f5d_b_y);
f5e = abs(f5e_a_y - f5e_b_y);

% nir
%n2a = abs(n2a_a_y - n2a_b_y);
n2b = abs(n2b_a - n2b_b_y);
n2c = abs(n2c_a - n2c_b_y);
n2d = abs(n2d_a - n2d_b_y);
n2e = abs(n2e_a - n2e_b_y);

n5a = abs(n5a_a_y - n5a_b_y);
n5b = abs(n5b_a_y - n5b_b_y);
n5c = abs(n5c_a_y - n5c_b_y);
n5d = abs(n5d_a_y - n5d_b_y);
n5e = abs(n5e_a_y - n5e_b_y);

figure;
%subplot(2,2,1)
plot(f_x_a, f5b)
% hold on
% plot(f_x_a, f2c)
% hold on
% plot(f_x_a, f2d)
% hold on
% plot(f_x_a, f2e)

