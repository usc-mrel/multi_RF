%% Simulation of concomitant fields using Bloch Siegert shift
%%% proposed correction method for 2D spiral with measured GIRFs
% Ziwei Zhao 

%% Clean slate
close all; 
clear all; 
clc;

%% Add paths
addpath(genpath('.../multi_RF/third_party/reVERSE-GIRF/original_source'));
addpath('.../multi_RF/third_party/lsqrSOL');
addpath(genpath('.../multi_RF/third_party/reVERSE-GIRF'));
addpath(genpath('.../multi_RF/thirdparty'));
addpath(genpath('.../multi_RF/third_party/Bloch_simulator'));

%% Constant definitions
gamma_uT = 267.5221;       % [rad/sec/uT]
gamma_mT = gamma_uT * 1e3; % [rad/sec/uT] * [1e3uT/mT] => *1e3 [rad/sec/mT]

flip = 90;                 % total flip angle [degrees]
flip = flip * pi / 180;

%% Design initial spiral trajectory
T          = 26e-3;         % pulse duration [sec]
dt         = 6.4e-6;        % RF/gradient raster time [sec]
t          = [0:dt:T-dt]';  % seconds, time vector
dr         = 0.25;          % resolution of trajectory [cm]
kmax       = 1/2/dr;        % [cycles/cm], max radius in k-space
fov        = 12.8;          % XFOV of unaliased excitation [cm] % 38.4 cm

beta       = 2;             % determines the spatial resolution of the selective volume
r0         = [0 0].';       % offset [m]

% spiral waveforms 
N = kmax/(1/fov);           % number of turns, no acceleration
k = 2*pi*kmax*(1-t/T).*exp(1i*2*pi*N*(1-t/T)); % [rad/cm] - kx = real(k); ky = imag(k);
k = [real(k) imag(k) zeros(length(k),1)];
k = k*100;                  % [rad/cm] -> [rad/m]

figure; plot(k(:,1), k(:,2), '-');
figure; plot(t, k(:,1)); hold on; plot(t, k(:,2));

K = k;

%% Load in B0 & B1 field maps + FOV info (example data from 7T 8ch head coil)
load('.../multi_RF/third_party/phase_relaxed_CPMG_excitation/b0_b1_maps_2d.mat');
% b0   : B0 map [uT], 64 x 64 (double)
% tx   : relative transmit sensitivity map, 64 x 64 x 8 (double, complex) 8 channels 
% X,Y,Z: spatial coordinates [m], 64 x 64 (double)
% m    : voxel mask, 64 x 64 (logical)

FOVx = max(X(:)) - min(X(:)); % [m]
FOVy = max(Y(:)) - min(Y(:)); % [m]
Z = Z * 0;
zoff = 10e-2;                 % [m]

[N1, N2, Nc] = size(tx);
xlist = cat(2, X(:), Y(:), Z(:) + zoff); % N1 x N2 x 3
idx = find(m); % index of non-masked voxels

%% Define a square target
% Offsets [m]
xoff = +6e-3;  % [mm] -> [m]
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

%% Define gradient correction model (GIRF)
girf = load('.../multi_RF/third_party/reVERSE-GIRF/GIRF_3T_London_20140729.mat');

% This correction gives G and k after convolution with GIRF
gcor = @(x)(gradient_distort_GIRF(x, girf.ff, girf.Hw, dt, 10));

%% Example design: Include GIRF
% Set default options
concomitant_correct = 0;    % switch between original and proposed methods
opt        = reVERSE_init;
dt         = opt.dt;        % sampling dwell time [usec]
opt.lambda = 1;
opt.Nstop  = 20;
opt.show   = 1;

B0_list   = [0.55];         % [tesla]
nr_B0     = length(B0_list);

alpha = 0.5;
g     = 0;

T1 = 1e6; % T1 relaxation time [sec]
T2 = 1e6; % T2 relaxation time [sec]

mxyz_offcenter = zeros(N1, N2, 3, nr_B0, 'double');

for ii = 1 : nr_B0
    
    B0 = B0_list(ii);   % main magnetic strength [T]

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
        
        mxyz_offcenter(idx1, idx2, :, ii) = cat(3, mx, my, mz);
        fprintf('done! (%5.4f sec)\n', toc(start_time));
    end
    
      % save the results
      cur_dir = pwd;    
      save(sprintf('bloch_B0%.2f_offc%.1fcm_iter%d_dur%.3f_dr%.2fcm_lambda%1.0f.mat', B0, zoff*1e2, opt.Nstop, duration, dr, opt.lambda), 'mxyz_offcenter');
      cd(cur_dir);
end

ind_t = 0:dt:(length(G)-1)*dt;
ind_t = ind_t * 1e3;

figure;  plot(ind_t, abs(rf(:,1)), 'LineWidth', 2);
hold on; plot(ind_t, abs(rf(:,2)), 'LineWidth', 2);
hold on; plot(ind_t, abs(rf(:,3)), 'LineWidth', 2);
hold on; plot(ind_t, abs(rf(:,4)), 'LineWidth', 2);
hold on; plot(ind_t, abs(rf(:,5)), 'LineWidth', 2);
hold on; plot(ind_t, abs(rf(:,6)), 'LineWidth', 2);
hold on; plot(ind_t, abs(rf(:,7)), 'LineWidth', 2);
hold on; plot(ind_t, abs(rf(:,8)), 'LineWidth', 2);

legend('Channel - 1', 'Channel - 2', 'Channel - 3', 'Channel - 4', 'Channel - 5',...
    'Channel - 6', 'Channel - 7', 'Channel - 8');
ylabel('Amplitude [mT]'); xlabel('time [ms]');
set(gca, 'FontSize', 16); title('Multi-channel RF @0.55T');

figure; subplot(1,2,1);
plot(ind_t, Gedd(:,1), 'LineWidth', 2);
hold on; plot(ind_t, Gedd(:,2), 'LineWidth', 2);
legend('Gx', 'Gy');
xlabel('time [ms]'); ylabel('[mT/m]');
xlim([0 18]); box off; grid on;
set(gca, 'FontSize', 18);
title('Gradient waveforms');
subplot(1,2,2);
plot(k(:,1), k(:,2), 'LineWidth', 2);
box off; grid on;
xlabel('[rads/m]'); ylabel('[rads/m]');
xlim([-1500 1500]); ylim([-1500 1500]);
set(gca, 'FontSize', 18); 
title('Excitation K-space');

% check NRMSE
mxy  = squeeze(complex(mxyz_offcenter(:,:,1,:), mxyz_offcenter(:,:,2,:)));
mxy_ori = mxy;
NRMSE = sqrt(sum(sum((abs(mxy_ori) - abs(P_)).^2)))/ sqrt(sum(sum(abs(P_).^2)));

%% Display the excitation pattern
% load re-verse method results
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
