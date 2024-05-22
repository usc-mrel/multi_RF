% calculate_cartesian_spin_echo_imaging_parameters.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 07/06/2023, Last modified: 09/18/2023

%% Gradient mode
switch grad_mode
    case 'Fast'
        max_grad = 24;      % Max gradient strength [mT/m]
        max_slew = 180.18;  % Maximum slew rate [mT/m/ms]
    case 'Normal'
        max_grad = 22;      % Max gradient strength [mT/m]
        max_slew = 100;     % Maximum slew rate [mT/m/ms]
    case 'Whisper'
        max_grad = 22;      % Max gradient strength [mT/m]
        max_slew = 50;      % Maximum slew rate [mT/m/ms]
end

%--------------------------------------------------------------------------
% Derated
%--------------------------------------------------------------------------
max_grad_derated = max_grad / sqrt(3);
max_slew_derated = max_slew / sqrt(3);

%% Set system limits
sys = mr.opts('MaxGrad'       , max_grad, 'GradUnit', 'mT/m' , ...
              'MaxSlew'       , max_slew, 'SlewUnit', 'T/m/s', ...
              'rfRingdownTime', 20e-6 , ...
              'rfDeadTime'    , 100e-6, ...
              'adcDeadTime'   , 10e-6 , ...
              'B0'            , B0);

sys_derated = mr.opts('MaxGrad'       , max_grad_derated, 'GradUnit', 'mT/m' , ...
                      'MaxSlew'       , max_slew_derated, 'SlewUnit', 'T/m/s', ...
                      'rfRingdownTime', 20e-6 , ...
                      'rfDeadTime'    , 100e-6, ...
                      'adcDeadTime'   , 10e-6 , ...
                      'B0'            , B0);

%% Calculate the real dwell time [sec]
%--------------------------------------------------------------------------
% real dwell time [sec]
% IDEA p219: dRealDwellTime denotes the dwell time with oversampling
%--------------------------------------------------------------------------
% round-down dwell time to 100 ns (sys.adcRasterTime  = 100 ns)
real_dwell_time = round((1 / bandwidth) / (readout_os_factor * base_resolution) * 1e7) * 1e-7;

%% Calculate the readout bandwidth [Hz]
readout_bandwidth = 1 / real_dwell_time; % [Hz]

%% Calculate the number of phase-encoding steps
nr_phase_encoding_steps_1 = base_resolution * phase_resolution * 1e-2;

%% Calculate the duration of an ADC event
adc_samples  = base_resolution * readout_os_factor; % true number of samples
adc_duration = adc_samples * real_dwell_time;

%% Display input parameters
fprintf('----------------------- Resolution (Common) -----------------------\n');
fprintf(' FoV read             = %4.0f [mm]\n', fov_read * 1e3);
fprintf(' FoV phase            = %4.0f [%%]\n', fov_phase);
fprintf(' Slice thickness      = %4.0f [mm]\n', slice_thickness * 1e3);
fprintf(' Base resolution      = %4.0f\n'     , base_resolution);
fprintf(' Phase resolution     = %4.0f [%%]\n', phase_resolution);
fprintf('------------------------------------------------------------------\n');

fprintf('----------------------- Geometry (Common) -----------------------\n');
fprintf(' Slices               = %d\n', nr_slices);
fprintf(' Series               = %s\n', series);
fprintf('------------------------------------------------------------------\n');

fprintf('----------------------- Sequence (Special) -----------------------\n');
fprintf(' bandwidth             = %6.0f [Hz/Px]\n', bandwidth);
fprintf(' Acquisition duration  = %6.4f [ms]\n', adc_duration * 1e3);
fprintf(' Real dwell time       = %6.0f [ns]\n', real_dwell_time * 1e9);
fprintf(' Max gradient strength = %6.0f [mT/m]\n', max_grad);
fprintf(' Max slew rate         = %6.0f [mT/m/ms]\n', max_slew);
fprintf('------------------------------------------------------------------\n');
