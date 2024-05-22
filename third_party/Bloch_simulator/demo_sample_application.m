% demo_sample_application.m
% Written by Namgyun Lee
% Email: ggang56@gmail.com
% Started: 08/01/2018, Last modified: 08/01/2018

%% Clean slate
close all; clear all; clc; 

%% Add paths
addpath D:\mfiles_nam\rf_pulse_design_code\Bloch_simulator;
addpath D:\mfiles_nam\ee469b_rf_pulse_design_for_MRI_code\rf_tools_nam;
LineWidth = 1;

%% Define parameters
dt  = 4e-3;       % RF/gradient dwell time [msec] (e.g: 4 usec)
TBW = 12;         % time-bandwidth product [msec] * [kHz]
T   = 2.4;        % pulse duration [msec]
flip_angle = 60;  % flip angle [degrees]
slab_thickness = 1.2; % [cm]

nt = ceil(T / dt); % number of samples
if mod(nt,2) == 1, nt = nt + 1; end
T  = nt * dt;

%% Calculate the windowed sinc RF pulse
rf = (flip_angle*pi/180) * wsinc(TBW, nt);
rfs = [zeros(nt/2,1); rfscaleg(rf, T); zeros(nt/2,1)];

% [kHz/G] * [G/cm] * [cm] = [kHz]
% gamma_bar * G * slab_thickness = BW
BW = TBW / T; % bandwidth [kHz]
G = BW / (4.257 * slab_thickness); % [kHz] / ([kHz/G] * [cm]) => [G/cm]
g = [-G * ones(nt/2,1); G * ones(nt,1); -G * ones(nt/2,1)]; % [G/cm]

%% Simulate the Bloch equation over a frequency range and a spatial range
z = (-2:.1:2).';        % in cm
f = (-200:1:200).';     % in Hz
t = (1:2*nt).' * dt;    % in msec
nvz = length(z);
nf = length(f);

% Bloch Simulation
fprintf('Performing the Bloch simulation\n');
mx0 = zeros(nvz,nf, 'double'); my0 = zeros(nvz,nf, 'double'); mz0 = ones(nvz,nf, 'double');
[mx,my,mz] = bloch(rfs, g, t*1e-3, 100, 100, f, z, 0, mx0, my0, mz0); 

% Transverse Magnetization
mxy = mx + 1j*my;

%% Display RF and gradient waveforms
figure('Color', 'w');
subplot(2,1,1); plot(t, rfs, 'b', 'LineWidth', LineWidth); grid on;
xlabel('time (ms)'); ylabel('RF (G)'); xlim([0 t(end)]); ylim([-0.05 0.2]);
subplot(2,1,2); plot(t, g, 'b', 'LineWidth', LineWidth); grid on;
xlabel('time (ms)'); ylabel('Gradient (G/cm)'); xlim([0 t(end)]); ylim([-1.2 1.2]);

%% Display the simulation result
[Z,F] = ndgrid(z, f);
figure('Color', 'w'); surf(F, Z, abs(mxy), 'EdgeColor', 'none'); colormap(gray(256));
xlabel('Frequency (-200 to 200 Hz)'); ylabel('Position (-2 to 2 cm)'); view(0, 90); axis square;
