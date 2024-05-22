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
z = (0).';        % in cm
f = (0).';     % in Hz
t = (1:2*nt).' * dt;    % in msec
nvz = length(z);
nf = length(f);

% T1 and T2 relaxation times in seconds
T1 = 100e-3;
T2 = 50e-3;

T1 = 1e9;
T2 = 1e9;


% Bloch Simulation
fprintf('Performing the Bloch simulation\n');
mx0 = zeros(nvz,nf, 'double'); my0 = zeros(nvz,nf, 'double'); mz0 = -ones(nvz,nf, 'double');
tic; [mx,my,mz] = bloch(rfs, g, dt*1e-3, T1, T2, f, z, 2, mx0, my0, mz0); toc;

% Transverse Magnetization
mxy = mx + 1j*my;

%%
% plot
t = (1:length(mx)).' * dt;
figure;
plot(t*1e3, mx, t*1e3, my, t*1e3, mz);
legend('Mx', 'My', 'Mz');
xlabel('Time (ms)');

%% Display RF and gradient waveforms
figure('Color', 'w');
subplot(3,1,1); plot(t, rfs, 'b', 'LineWidth', LineWidth); grid on;
xlabel('time (ms)'); ylabel('RF (G)'); xlim([0 t(end)]); ylim([-0.05 0.2]);
subplot(3,1,2); plot(t, g, 'b', 'LineWidth', LineWidth); grid on;
xlabel('time (ms)'); ylabel('Gradient (G/cm)'); xlim([0 t(end)]); ylim([-1.2 1.2]);
subplot(3,1,3); plot(t*1e3, mx, t*1e3, my, t*1e3, mz, 'LineWidth', LineWidth); grid on;
legend('Mx', 'My', 'Mz'); xlabel('Time (ms)');




