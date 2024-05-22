% demo_test_orientation.m
% Written by Namgyun Lee
% Email: ggang56@gmail.com
% Started: 04/20/2017, Last modified: 04/20/2017
% From Miki Lustig's class homework

%% Clean slate
close all; clear all; clc;

dt = 4e-6;   % sampling rate [sec]
T  = 300e-3; % pulse duration [sec]
nt = T / dt; % number of pulse samples

%--------------------------------------------------------------------------
% flip angle = gamma * B1 * dt
%      [rad] = [Hz/mT] * [2*pi rad/cycle] * [mT/10G] * [G] * [sec]
% THus, B1 = flip angle / ((gamma*2*pi/10) * dt);
%      [G] = [rad] / ([Hz/mT] * [rad/cycle] * [mT/10G] * [sec])
%--------------------------------------------------------------------------
% impulse RF pulse (nt x 1) [G]
flip_angle = 90 * pi / 180; % [radians]
gamma = 42577.46778; % gyromagnetic ratio [Hz/mT]
b1 = zeros(nt,1);
%b1(1) = flip_angle / ((gamma*2*pi/10) * dt) * exp(1j * pi / 4);
% b1(1) = flip_angle / ((gamma*2*pi/10) * dt) *  exp(1j * 0 / 2);
% b1(1)
b1(1) = flip_angle / ((gamma*2*pi/10) * dt) *  exp(1j * pi / 2);
b1(1)
b1(1) = flip_angle / ((gamma*2*pi/10) * dt) *  exp(1j * pi);
b1(1)


% b1(1000) = 60 * pi / 180 / ((gamma*2*pi/10) * dt);
% b1(2000) = 60 * pi / 180 / ((gamma*2*pi/10) * dt);
% b1(3000) = 60 * pi / 180 / ((gamma*2*pi/10) * dt);
0
% no gradient (nt x 2) [G/cm]
g = zeros(nt,1);

% the spin is on-resonance (Nx1) [Hz]
df = 100;
% df = 100;

% the spin at iso-center
dp = 0;

% T1 and T2 relaxation times in seconds
t1 = 100e-3;
t2 = 50e-3;

% Start at Mz
mx0 = 0;
my0 = 0;
mz0 = 1;

% Simulate
tic; [mx,my,mz] = bloch(b1, g, dt, t1, t2, df, dp, 2, mx0, my0, mz0); toc;
%%
% plot
t = (1:nt).' * dt;
figure;
plot(t*1e3, mx, t*1e3, my, t*1e3, mz);
legend('Mx', 'My', 'Mz');
xlabel('Time (ms)');

return;

%% Visualization
visualizeMagn(zeros(size(mx)), zeros(size(mx)), [mx,my,mz], 100);
