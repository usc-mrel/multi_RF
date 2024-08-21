% demo_Pauly_1989_JMR_spiral2d_concomitant_fields.m
% Written by Namgyun Lee, Modified by Ziwei Zhao

%% Clean slate
close all; clear all; clc;

%% Add paths
setup_path;

%% B0 sweep 
B0_swp = [0.2 0.55 1.5 3 7];

% uncomment this line to generate dense samples to reproduce figure 2B
%B0_swp = [0.005 0.01 0.03 0.05 0.09 0.1 0.2 0.3 0.5 0.7 0.9 1.0 1.5 3.0 7.0];

for iB0 = 1:length(B0_swp)

%% Define parameters
XFOV       = 8e-2;       % XFOV of unaliased excitation [m]
dr         = 0.5e-2;     % resolution of a k-space trajectory [m]
T          = 8e-3;       % pulse duration [sec]
flip_angle = 90;         % flip angle [degrees]
beta       = 2;          % determines the spatial resolution of the selective volume
r0         = [1;1]*1e-2; % offset [m]

B0    = B0_swp(iB0);     % main magnetic stregnth [T]
alpha = 0.5;
g     = 0;

%% Define constants
% [Hz/mT] * [2pi rad/cycle] * [1e3mT/T] => [rad/sec/T]
gamma = 42577.46778 * (1e3 * 2 * pi); % gyromagnetic ratio for 1H [rad/sec/T]

%% Calculate a sampling interval in time [sec]
w0 = -gamma * B0;                    % [rad/sec/T] * [T] => [rad/sec]
dt_max = (2 * pi) /(2 * abs(w0));    % [sec]
dt = 7e-10;

%% Design a constant-angular-rate spiral excitation k-space trajectory
%--------------------------------------------------------------------------
% dr = 2*pi/(2*kmax), where dr in [m], kmax in [rad/m]
% => kmax = pi / dr 
% XFOV = 2*pi / dk, dk = 2*kmax / (2*n)
% XFOV = 2*pi / dk = 2*pi * (2*n) / (2*kmax) = 2*n*dr [m]
% n = XFOV / (2*dr)
%--------------------------------------------------------------------------
fprintf('Calculating a constant-angluar-rate spiral excitation k-space trajectory\n');
Nt = floor(T/dt) + 1;  % number of pulse samples
A  = pi / dr;          % maximum radius in k-space [rad/m]
n  = XFOV / (2 * dr);  % number of cycles
T  = dt * Nt;          % pulse duration [msec]

% Calculate a time axis [sec]
t  = (0:Nt-1).' * dt;

% Calculate the FWHM: FWHM = 2*sqrt(2*log(2))*sigma, sigma = sqrt(2)*beta/A
FWHM = 4 * sqrt(log(2)) * beta / A;

fprintf('Field of excitation (XFOV) = %3.2f [cm]\n', XFOV*1e2);
fprintf('spatial resolution dr      = %3.2f [cm]\n', dr*1e2);
fprintf('kmax = 1 / (2*dr)          = %3.2f [cycle/cm]\n', A/(2*pi*1e2));
fprintf('n cycles = XFOV / (2*dr)   = %3.2f\n', n);
fprintf('RF/gradient dwell time     = %4.3f [usec]\n', dt*1e6);
fprintf('Pulse duration             = %3.3f [msec]\n', T*1e3);
fprintf('FWHM                       = %3.3f [cm]\n', FWHM*1e2);

%% Calculate a k-space trajectory in the logical coordinate system [rad/m]
kx = A * (1 - t / T) .* cos(2 * pi * n * t / T);
ky = A * (1 - t / T) .* sin(2 * pi * n * t / T);
kz = zeros(Nt,1, 'double');

%% Calculate gradient waveforms in the logical coordinate system [G/cm]
% [rad/m] * [m/1e2cm] / ([rad/sec/T] * [T/1e4G] * [sec]) => [G/cm]
Gx = -(A * 1e-2)/ (gamma * T * 1e-4) * (2 * pi * n * (1 - t / T) .* sin(2 * pi * n * t / T) + cos(2 * pi * n * t / T));
Gy =  (A * 1e-2)/ (gamma * T * 1e-4) * (2 * pi * n * (1 - t / T) .* cos(2 * pi * n * t / T) - sin(2 * pi * n * t / T));
Gz = zeros(Nt,1, 'double');

%% Calculate an RF waveform that produces a cylindrical Gaussian weighting of k-space
B1 = (A * 1e-2) / (gamma * T * 1e-4) * exp(-beta^2 * (1 - t / T).^2) .* sqrt((2 * pi * n * (1 - t / T)).^2 + 1);

% Scale the RF to produce a given tip angle
%--------------------------------------------------------------------------
% flip angle for the kth sample = gamma * b1(k) * dt
% [rad] = [rad/sec/T] * [T/1e4G] * [G] * [sec]
% Thus, b1(k) = flip angle for the kth sample / (gamma * dt)
%         [G] = [rad] / ([rad/sec/T] * [T/1e4G] * [sec])
%--------------------------------------------------------------------------
scale_factor = (flip_angle * pi / 180) / (gamma * 1e-4 * sum(B1) * dt);
B1 = scale_factor * B1; % [G]

%% Shift the selective volume to another spatial position
B1 = B1 .* exp(-1j * (kx * r0(1) + ky * r0(2)));  % [G]

%% Calculate Cartesian spatial locations [m]
FOVx = 5e-2;      % FOV(row,x) of a desired pattern [m]
FOVy = 5e-2;      % FOV(col,y) of a desired pattern [m]
dx   = 0.2e-2;    % resolution(row,x) of a desired pattern [m]
dy   = 0.2e-2;    % resolution(col,y) of a desired pattern [m]
dz   = 0.2e-2;    % resolution(slice,z) of a desired pattern [m]
Nx   = FOVx / dx; % number of samples(row,x) in the desired pattern
Ny   = FOVy / dy; % number of samples(col,y) in the desired pattern
Nz   = 1;
Nv   = Nx * Ny;   % total number of voxels

z_offset = 20 * 1e-2;                         % [m] % please change it based on different settings
x_range = (-floor(Nx/2):ceil(Nx/2)-1).' * dx; % [m]
y_range = (-floor(Ny/2):ceil(Ny/2)-1).' * dy; % [m]
z_range = (-floor(Nz/2):ceil(Nz/2)-1).' * dz; % [m]
z_range = z_range + z_offset;
[x,y,z] = ndgrid(x_range, y_range, z_range);

%% Bloch simulation
df  = 0; % [Hz]
mx0 = zeros(Nx,Ny, 'double');
my0 = zeros(Nx,Ny, 'double');
mz0 =  ones(Nx,Ny, 'double');

mxyz_concomitant    = complex(zeros(Nx, Ny, 3, 'double'));
mxyz_rf             = complex(zeros(Nx, Ny, 3, 'double'));
mxyz_concomitant_rf = complex(zeros(Nx, Ny, 3, 'double'));
mxyz_bloch_siegert  = complex(zeros(Nx, Ny, 3, 'double'));

G = cat(2, Gx, Gy, Gz);

start_time = tic;
for idx2 = 1:Ny
    for idx1 = 1:Nx
        fprintf('Performing Bloch simulation (%3d/%3d,%3d/%3d)... ', idx1, Nx, idx2, Ny);
        dp = cat(2, x(idx1,idx2), y(idx1,idx2), z(idx1,idx2)) * 1e2; % [m] * [1e2cm/m] => [cm]

        %------------------------------------------------------------------
        % Calculate the concomitant fields in the rotating frame
        %------------------------------------------------------------------
        % [G/cm] * [1e2cm/m] * [m] => [G]
        Bcx_reference_frame = -0.5 * (Gz * 1e2) * x(idx1,idx2) + (Gx * 1e2) * z(idx1,idx2);
        Bcy_reference_frame = -0.5 * (Gz * 1e2) * y(idx1,idx2) + (Gy * 1e2) * z(idx1,idx2);
        Bcx_rotating_frame =  Bcx_reference_frame .* cos(w0 * t) + Bcy_reference_frame .* sin(w0 * t); 
        Bcy_rotating_frame = -Bcx_reference_frame .* sin(w0 * t) + Bcy_reference_frame .* cos(w0 * t); 
        Bc = -0.5 * (Gz * 1e2) * x(idx1,idx2) + (Gx * 1e2) * z(idx1,idx2) + 1j * (-0.5 * (Gz * 1e2) * y(idx1,idx2) + (Gy * 1e2) * z(idx1,idx2));
        Bc = Bc .* exp(-1j * w0 * t);

        %------------------------------------------------------------------
        % Perform Bloch simulation for concomitant fields
        %------------------------------------------------------------------
        [mx_concomitant,my_concomitant,mz_concomitant] = bloch(Bc, G, dt, 100, 100, df, dp, 0, mx0(idx1,idx2), my0(idx1,idx2), mz0(idx1,idx2));
        mxyz_concomitant(idx1,idx2,:) = cat(3, mx_concomitant, my_concomitant, mz_concomitant);

        %------------------------------------------------------------------
        % Perform Bloch simulation for RF
        %------------------------------------------------------------------
        [mx_rf,my_rf,mz_rf] = bloch(B1, G, dt, 100, 100, df, dp, 0, mx0(idx1,idx2), my0(idx1,idx2), mz0(idx1,idx2));
        mxyz_rf(idx1,idx2,:) = cat(3, mx_rf, my_rf, mz_rf);

        %------------------------------------------------------------------
        % Perform Bloch simulation for concomitant fields + RF
        %------------------------------------------------------------------
        [mx_concomitant_rf,my_concomitant_rf,mz_concomitant_rf] = bloch(B1+Bc, G, dt, 100, 100, df, dp, 0, mx0(idx1,idx2), my0(idx1,idx2), mz0(idx1,idx2));
        mxyz_concomitant_rf(idx1,idx2,:) = cat(3, mx_concomitant_rf, my_concomitant_rf, mz_concomitant_rf);

        %------------------------------------------------------------------
        % Perform Bloch simulation considering concomitant fields
        % as a Bloch-Siegert shift
        %------------------------------------------------------------------
        [mx_bloch_siegert,my_bloch_siegert,mz_bloch_siegert] = bloch_maxwell(B1, G, dt, 100, 100, df, dp, 0, B0, alpha, g, mx0(idx1,idx2), my0(idx1,idx2), mz0(idx1,idx2));
        mxyz_bloch_siegert(idx1,idx2,:) = cat(3, mx_bloch_siegert, my_bloch_siegert, mz_bloch_siegert);
        fprintf('done! (%5.4f sec)\n', toc(start_time));
    end
end

%% save magnetization obtained with concomitant fields, RF, and conconmitant fields and RF
mxy_concomitant    = complex(mxyz_concomitant(:,:,1), mxyz_concomitant(:,:,2));
mxy_rf             = complex(mxyz_rf(:,:,1), mxyz_rf(:,:,2));
mxy_concomitant_rf = complex(mxyz_concomitant_rf(:,:,1), mxyz_concomitant_rf(:,:,2));
mxy_bloch_siegert  = complex(mxyz_bloch_siegert(:,:,1), mxyz_bloch_siegert(:,:,2));

%% display NRMSE -- ziwei
diff_real = sum((real(mxy_bloch_siegert) - real(mxy_concomitant_rf)).^2) / sum((real(mxy_concomitant_rf).^2));
diff_imag = sum((imag(mxy_bloch_siegert) - imag(mxy_concomitant_rf)).^2) / sum((imag(mxy_concomitant_rf).^2));

diff_mag = sum((abs(mxy_bloch_siegert) - abs(mxy_concomitant_rf)).^2) / sum((abs(mxy_concomitant_rf).^2));
diff_phase = sum((angle(mxy_bloch_siegert)*180/pi - angle(mxy_concomitant_rf)*180/pi).^2) / sum(((angle(mxy_concomitant_rf)*180/pi).^2));

fprintf('main_field     = %4.3f [T]\n', B0);
fprintf('diff_real      = %4.3f \n', diff_real);
fprintf('diff_imag      = %3.3f \n', diff_imag);
fprintf('diff_mag       = %3.3f \n', diff_mag);
fprintf('diff_phase     = %3.3f \n', diff_phase);

diff.real  = diff_real;
diff.imag  = diff_imag;
diff.mag   = diff_mag;
diff.phase = diff_phase;

%% save results
cd('../figures/sim_results/Fig2_validation/');
% if sample densely from 0.005T to 7T
%cd('../figures/sim_results/Fig2_validation/densesamples/');

save(sprintf('results_%.3fT_7e_10dt_25Nxy.mat', B0), 'mxyz_concomitant', 'mxyz_rf', 'mxyz_concomitant_rf', 'mxyz_bloch_siegert', 'diff');

end
