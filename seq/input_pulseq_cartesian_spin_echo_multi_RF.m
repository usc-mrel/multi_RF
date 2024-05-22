% input_pulseq_cartesian_spin_echo_multislice_protocol1.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 02/05/2023, Last modified: 09/18/2023

%% Contrast
%--------------------------------------------------------------------------
% Common
%--------------------------------------------------------------------------
TR         = 1000e-3;  % Total acquisition time [sec]
TE         = 84e-3;    % Echo time [sec]
flip_angle = 90;       % Flip angle [deg]

%--------------------------------------------------------------------------
% Dynamic
%--------------------------------------------------------------------------
nr_averages = 1;       % Averages

%% Resolution
%--------------------------------------------------------------------------
% Common
%--------------------------------------------------------------------------
fov_read         = 256e-3;   % FoV read [m]
fov_phase        = 100;      % FoV phase [%]
slice_thickness  = 5e-3;     % Slice thickness [m]
base_resolution  = 256;      % Base resolution
phase_resolution = 100;      % Phase resolution [%]

%% Geometry 
%--------------------------------------------------------------------------
% Common
%--------------------------------------------------------------------------
nr_slices = 1;              % Number of slices in this slice group

%--------------------------------------------------------------------------
% Dist.Factor
% The gap between the slices of a slice group or the slabs of a slap group is
% expressed as a percentage of the slice thickness of slab thickness. At 100%,
% the gap between the slices/slabs is exactly one slice/slab thickness. Negative
% values lead to overlapping of slices.
%--------------------------------------------------------------------------
dist_factor = 700;            % Slice distance in this slice group [%]

%--------------------------------------------------------------------------
% Series
% The Series measurement parameters influence the sequence of slice processing.
% - Ascending: The slices are excited starting at the beginning of the slice
% or slab group (start -> end). e.g.) 1,2,3,4,5
% - Descending: The slices are excited starting at the end of the slice or slab
% group (end -> start). e.g.) 5,4,3,2,1
% - Interleaved: 1,3,5,2,4
%--------------------------------------------------------------------------
series = 'Interleaved'; % 'Ascending', 'Descending', 'Interleaved'

%% Sequence (Special)
%--------------------------------------------------------------------------
% Part 1
%--------------------------------------------------------------------------
bandwidth = 610;  % Bandwidth [Hz/Px]
bandwidth = 122;  % Bandwidth [Hz/Px]

%--------------------------------------------------------------------------
% Part 2
%--------------------------------------------------------------------------
grad_mode = 'Fast'; % Gradient mode: 'Whisper', 'Normal', 'Fast'

%% Set the readout oversampling factor
readout_os_factor = 1; % readout oversampling factor

%% Define imaging parameters
%--------------------------------------------------------------------------
% Define parameters for a hyperbolic secant pulse
%--------------------------------------------------------------------------
ref_duration  = 31e-3;       % RF duration [sec]
OMEGA0        = 520;         % [Hz]
mu            = 9.5;         % dimensionless parameter
ref_bandwidth = 810;         % bandwidth of an adiabatic pulse [Hz]

%% Define main field strength [T]
B0 = 0.55;

%% Define an output directory
output_directory = fullfile(pwd, mfilename);
