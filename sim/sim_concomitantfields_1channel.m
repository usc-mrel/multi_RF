%% Simulation of concomitant fields using Bloch Siegert shift
%%% proposed correction method for 2D spiral with measured GIRFs
% single channel simulation at 0.55T

% Ziwei Zhao

%% Clean slate
% close all; 
% clear all;
% clc;

%% Setup paths
% setup_path;

%% field strengths and off-isocenter variables to change
% zoff = 00e-2;        % [m]
% concomitant_correct = 1;  % switch between original and proposed methods
% output_path = pwd; % sanity check

%% Constant definitions
gamma_uT = 267.5221;       % [rad/sec/uT]
gamma_mT = gamma_uT * 1e3; % [rad/sec/uT] * [1e3uT/mT] => *1e3 [rad/sec/mT]

flip = 90;                 % total flip angle [degrees]
flip = flip * pi / 180;

%% Design initial spiral trajectory
T          = 26e-3;         % pulse duration [sec]
dt         = 1e-6;          % RF/gradient raster time [sec] 3.2e-6 works
t          = [0:dt:T-dt]';  % seconds, time vector
dr         = 0.25;          % resolution of trajectory [cm]
kmax       = 1/2/dr;        % [cycles/cm], max radius in k-space
fov        = 12.8;          % XFOV of unaliased excitation [cm] 12.8 in the simulation

% spiral waveforms 
N = kmax/(1/fov);           % number of turns, no acceleration
k = 2*pi*kmax*(1-t/T).*exp(1i*2*pi*N*(1-t/T)); % [rad/cm] - kx = real(k); ky = imag(k);
k = [real(k) imag(k) zeros(length(k),1)];
k = k*100;                  % [rad/cm] -> [rad/m]

K = k;

%% Load in B0 & B1 field maps
% 15mm off isocenter
% load /Users/ziwei/Documents/matlab/multi_RF/b0b1_map_055T/F15_64/fieldmap_b0_fov220mm_matrix64_off_15cm_03312024.mat
% dicom_path = '/Users/ziwei/Documents/matlab/multi_RF/b0b1_map_055T/F15_64';

% 10mm off isocenter
% load /Users/ziwei/Documents/matlab/multi_RF/b0b1_map_055T/F10_64/fieldmap_b0_fov220mm_matrix64_off_10cm_03302024.mat
% dicom_path = '/Users/ziwei/Documents/matlab/multi_RF/b0b1_map_055T/F10_64';

% 5mm off isocenter
% load /Users/ziwei/Documents/matlab/multi_RF/b0b1_map_055T/F5_64/fieldmap_b0_fov220mm_matrix64_off_5cm_03312024.mat
% dicom_path = '/Users/ziwei/Documents/matlab/multi_RF/b0b1_map_055T/F5_64';

% 0mm isocenter
% load('/Users/ziwei/Documents/matlab/multi_RF/b0b1_map_055T/F0_64/fieldmap_b0_fov220mm_matrix64_02282024.mat');
% dicom_path = '/Users/ziwei/Documents/matlab/multi_RF/b0b1_map_055T/F0_64';

cd(dicom_path);
dicom_file = dir([dicom_path, '/*.IMA']);
b1_mag = double(dicomread(dicom_file.name))/10/79.99; % unitless scale 

% load('.../multi_RF/b0b1_map_055T/F15_64/mask_F15.mat'); % off isocenter 15cm
% load('.../multi_RF/b0b1_map_055T/F10_64/mask_F10.mat'); % off isocenter 10cm
% load('.../multi_RF/b0b1_map_055T/F5_64/mask_F5.mat');  % off isocenter  5cm
% load('/Users/ziwei/Documents/matlab/multi_RF/b0b1_map_055T/F0_64/mask_F0.mat');
% m = imresize(m,[64,64], Method='nearest'); % downsample the mask - only for F0 

% generate X Y and Z
pixel_spacing = 3.4375e-3; % [m]
FOV = 220e-3;              % [m]
base_resolution = FOV/pixel_spacing;

x  = linspace(-FOV/2, FOV/2, base_resolution); 
y  = x;

[X, Y] = meshgrid(x,y);
Y = flipud(Y);
Z = zeros(size(X));

%% field map
% 15mm 
% load .../multi_RF/b0b1_map_055T/F15_64/fieldmap_b0_fov220mm_matrix64_off_15cm_03312024.mat
% 10mm
% load .../multi_RF/b0b1_map_055T/F10_64/fieldmap_b0_fov220mm_matrix64_off_10cm_03302024.mat
% 5mm
% load .../multi_RF/b0b1_map_055T/F5_64/fieldmap_b0_fov220mm_matrix64_off_5cm_03312024.mat
% 0mm
% load('/Users/ziwei/Documents/matlab/multi_RF/b0b1_map_055T/F0_64/fieldmap_b0_fov220mm_matrix64_02282024.mat');

b0 = fieldmap .* m ./ (gamma_uT / (2*pi)); % [uT] 

%%
FOVx = max(X(:)) - min(X(:)); % [m]
FOVy = max(Y(:)) - min(Y(:)); % [m]
Z = Z * 0;

%% replace b1 mag
tx_phase = zeros(base_resolution);
tx_mag   = b1_mag; 
tx       = tx_mag .* exp(1i*tx_phase) .* m;

%% 
[N1, N2, Nc] = size(tx);
xlist = cat(2, X(:), Y(:), Z(:) + zoff); % N1 x N2 x 3
idx = find(m); % index of non-masked voxels

%% Define a square target
% Offsets [m]
xoff = +6e-3;  % [mm] -> [m] too minimum 
yoff = 0e-3;   % [mm] -> [m]

r0   = 15e-3;  % [mm] -> [m]
r0   = 1.74 * r0;

% Square beam
P  = double(((abs(X - xoff) < r0) & (abs(Y - yoff) < r0)));

% apodize
h  = fspecial('gaussian', 3, 0.5);
P_ = imfilter(P, h);

% Now scale to flip
P  = P_ * flip * 1j; % desired excitation pattern

%% apply GIRF at 0.55T
figure_out_transformation_matrix;
tRR  = 1;  
sR.R = R_gcs2dcs;
sR.T = 0.55;  % [T]
gcor = @(x)(apply_GIRF_tx(permute(x, [1 3 2]), dt, sR, tRR, root_path));

%% Example design: Include GIRF
% Set default options
if concomitant_correct 
    script_name = 'blochmex'; 
else   
    script_name = 'bloch'; 
end

opt        = reVERSE_init_test;
opt.dt     = dt;         % sampling dwell time [usec]
opt.lambda = 1;
opt.Nstop  = 20;
opt.show   = 1;

alpha = 0.5;
g     = 0;

T1 = 1e6; % T1 relaxation time [sec]
T2 = 1e6; % T2 relaxation time [sec]

B0 = 0.55;   % main magnetic strength [T]

mxyz_offcenter = zeros(N1, N2, 3, 'double');

% Set up function for system matrix
if concomitant_correct
    afun  = @(k, G)(STA_maxwell_system_matrix_con(xlist, k, G, B0, opt.dt * (1:size(k,1)), tx, b0, m, 'loopcalc'));
    tic;
    [bb, Gv] = reVERSE_GIRF_con(P(idx), K, afun, gcor, opt);
    Time_reVERSE_con = toc;
else
    afun2 = @(k)(STA_system_matrix(xlist, k, opt.dt * (1:size(k,1)), tx, b0, m, 'loopcalc'));
    tic;
    [bb,Gv] = reVERSE_GIRF(P(idx), K, afun2, gcor, opt);
    Time_reVERSE = toc;
end

% Outputs
rf = bb{end};   % [mT]
G  = Gv{end};   % [mT/m]

%% check slew rate and Gmax
% resample G based on sys.Gradrastertime
G_input = G(1:10:end,:)/10; 

% check Gmax and slewrate
G_check = (G_input(:,1) + 1i*G_input(:,2)) * 1e1; % [G/cm]
[ktmp,gtmp,s,m1,m2,t,v] = calcgradinfo(G_check, 1e-5);
max(abs(s))
max(G_input(:,1))
max(G_input(:,2))

duration = dt*1e3*length(G);

%% Perform Bloch simulation
N_lambda = length(idx); % number of voxels within a mask

mx0 = zeros([N1 N2], 'double');
my0 = zeros([N1 N2], 'double');
mz0 = zeros([N1 N2], 'double');
mz0(idx) = 1;

[Gedd, k] = gcor(G);  % [mT/m]

%%
% perform bloch equations
start_time = tic;
for lambda = 1 : N_lambda
    
    [idx1,idx2] = ind2sub([N1 N2], idx(lambda));
    fprintf('Performing Bloch simulation at (%3d/%3d)... ', lambda, N_lambda);
    
    % rf_out: [mT]   * [T/1e3mT] * [1e4G/T]             => *1e1 [G]
    % G_out : [mT/m] * [T/1e3mT] * [1e4G/T] * [m/1e2cm] => *1e-1[G/cm]

    rf_combined = sum(bsxfun(@times, rf, reshape(tx(idx1,idx2,:), [1 Nc])), 2); 
    df = b0(idx1, idx2) * gamma_uT / (2 * pi); % [uT] * [rad/sec/uT] * [cycle/2pi rad] => [Hz]
            
    % off iso-center
    dp = cat(3, X(idx1,idx2), Y(idx1,idx2), Z(idx1,idx2) + zoff) * 1e2; % [m] * [1e2cm/m] => [cm]
    
    [mx,my,mz] = bloch_maxwell(rf_combined*1e1, Gedd*1e-1, dt, T1, T2, df, dp, 0, B0, alpha, g, mx0(idx1,idx2), my0(idx1,idx2), mz0(idx1,idx2));
    
    mxyz_offcenter(idx1, idx2, :) = cat(3, mx, my, mz);
    fprintf('done! (%5.4f sec)\n', toc(start_time));
end

% save the results
cur_dir = pwd; 
cd(output_path);
save(sprintf('%s_055T_1tx_offc%.1fcm_iter%d_dur%.3f_dr%.2fcm_rf_g_dt1e-6_fixcoord_mxyz.mat', script_name, zoff*1e2, opt.Nstop, duration, dr), 'rf', 'G', 'dt', 'mxyz_offcenter');
cd(cur_dir);

% check NRMSE
mxy = squeeze(complex(mxyz_offcenter(:,:,1,:), mxyz_offcenter(:,:,2,:)));
mxy_ori = mxy;
NRMSE = sqrt(sum(sum((abs(mxy_ori) - abs(P_)).^2)))/ sqrt(sum(sum(abs(P_).^2)))

%% Display the excitation pattern
% load re-verse method results
if 0 
mxy_offcenter = complex(mxyz_offcenter(:,:,1,:), mxyz_offcenter(:,:,2,:));

block = 1.01 * complex(ones(N1,1, 'double'), ones(N1,1, 'double'));
im_montage = cat(1, imag(cat(2, 1j*P_.', block, mxy_offcenter(:,:,1,1).', block, mxy_offcenter(:,:,1,1).', block, mxy_offcenter(:,:,1,1).', block, mxy_offcenter(:,:,1,1).', ...
                    block, mxy_offcenter(:,:,1,1).')), ...
                    1.01 * ones(1, (N2+1)*6-1, 'double'), ...
                    real(cat(2, 1j*P_.', block, mxy_offcenter(:,:,1,1).', block, mxy_offcenter(:,:,1,1).', block, mxy_offcenter(:,:,1,1).', block, mxy_offcenter(:,:,1,1).', ...
                    block, mxy_offcenter(:,:,1,1).')), ...
                    1.01 * ones(1, (N2+1)*6-1, 'double'));                
                
FontSize = 14;
cmap = cat(1, jet(256), [1 1 1]);                

figure('Color', 'w', 'Position', [-1 2 1101 814]);
imagesc(abs(im_montage)*flip*180/pi); axis image off; colormap(gca, cmap); 
hc = colorbar;

% Draw arrows
cx = floor(N1/2) + 1;
cy = floor(N2/2) + 1;
arrow([cx cy], [cx+20 cy], 'Color', 'w', 'Length', 7, 'TipAngle', 25, 'Width', 1);
arrow([cx cy], [cx cy-20], 'Color', 'w', 'Length', 7, 'TipAngle', 25, 'Width', 1);
text(cx+20, cy+1 , '$\bf{x}$', 'Color', 'w', 'Interpreter', 'latex', 'FontSize', FontSize, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
text(cx+5 , cy-16, '$\bf{y}$', 'Color', 'w', 'Interpreter', 'latex', 'FontSize', FontSize, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

cx = floor(N1/2) + 1;
cy = floor(N2/2) + 1;
arrow([cx cy+(N1+1)*2], [cx+20 cy+(N1+1)*2], 'Color', 'w', 'Length', 7, 'TipAngle', 25, 'Width', 1);
arrow([cx cy+(N1+1)*2], [cx cy-20+(N1+1)*2], 'Color', 'w', 'Length', 7, 'TipAngle', 25, 'Width', 1);
text(cx+20, cy+1+(N1+1)*2 , '$\bf{x}$', 'Color', 'w', 'Interpreter', 'latex', 'FontSize', FontSize, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
text(cx+5 , cy-16+(N1+1)*2, '$\bf{y}$', 'Color', 'w', 'Interpreter', 'latex', 'FontSize', FontSize, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end