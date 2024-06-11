% Demo_pulseq_Pauly_1989_JMR_modified.m
% Written by Nam Gyun Lee, modified by Ziwei Zhao

%% Clean slate
close all; clear all; clc;

%% Set source directories
addpath(genpath('.../multi_RF/third_party/pulseq'));
addpath(genpath('.../multi_RF/third_party/mintgrad'));

%% Load an input file
input_pulseq_cartesian_spin_echo_multi_RF;

%% Calculate cartesian spin echo imaging parameters
calculate_cartesian_spin_echo_imaging_parameters;
% pause;

%% Make an output directory
mkdir(output_directory);

%% Adiabatic (hyperbolic secant) SE pulse sequence
%-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%                                                                                                                  TR 
%    |<------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------>|
%    |                                                                                             TE                                                                                                  |                              |
%    |         |<------------------------------------------------------------------------------------------------------------------------------------------------------------------------->|----------------------------->|
%    |         |                   tau                                       tau                                                                 TE - (2 * tau)                            |                              |
%    |         |<--------------------------------------->|<--------------------------------------->|<------------------------------------------------------------------------------------->|                              |
%    |         |         |     |           |             |             |           |               |           |           |             |             |           |             |         |         |                    |
%    |         |         |     |           |           rf_ref          |           |               |           |           |           rf_ref          |           |             |         |         |                    |
%    |         |         |     |           |                           |           |               |           |           |                           |           |             |         |         |                    |
%    |         |         |     |           |            ___            |           |               |           |           |            ___            |           |             |         |         |                    |
%    | |       |rf_ex  | |     |           |           / | \           |           |               |           |           |           / | \           |           |             |         |         |                    |
%    | |      _|_      | |     |           | |        /  |  \        | |           |               |           |           | |        /  |  \        | |           |             |         |         |                    |
% RF | |     / | \     | |     |           | |       /   |   \       | |           |               |           |           | |       /   |   \       | |           |             |         |         |                    |
% ___|_|    /  |  \    |_|_____|___________|_|______/    |    \______|_|___________|_______________|___________|___________|_|______/    |    \______|_|___________|_____________|_________|_________|____________________|
%    | |\__/   |   \__/| |     |           | |<----------|---------->| |           |               |           |           | |<----------|---------->| |           |             |         |         |                    |
%    | |       |       | |     |         ->| |<-         |           | |           |               |           |         ->| |<-         |           | |           |             |         |         |                    |
%    | |       |       | |     |           | |rfDeadTime |           | |           |               |           |           | |rfDeadTime |           | |           |             |         |         |                    |
%    | |       |       | |     left_crusher| |           |           | |right_crusher              |           left_crusher| |           |           | |right_crusher            |         |         | gz_spoiler         |
%    | |               | |     |    ___    | |                       | |    ___    |               |           |    ___    | |                       | |    ___    |             |         |         |    ___    |        |
%    | |            gz_rephaser|   /   \   | |         gz_ref        | |   /   \   |               |           |   /   \   | |         gz_ref        | |   /   \   |             |         |         |   /   \   |        |
%    | |_______|_______| |     |  /     \  | |___________|___________| |  /     \  |               |<--------->|  /     \  | |___________|___________| |  /     \  |             |         |         |  /     \  |        |
% Gz | /       | gz_ex \ |     | /       \ | /           |           \ | /       \ |               |           | /       \ | /           |           \ | /       \ |             |         |         | /       \ |        |
% ___|/        |        \|     |/         \|/            |            \|/         \|_______________|___________|/         \|/            |            \|/         \|_____________|_________|_________|/         \|________|
%    |         |         \    /|           |                           |           |               |           |           |                           |           |             |         |         |           |        |
%    |         |         |\__/ |           |<------------------------->|           |               |           |           |<------------------------->|           |             |         |         |           |        |
%    |         |         |     |           |                           |           |    delay1     |   delay2  |           |                           |           |             |         |         |           |        |
%    |<----------------->|<--->|<--------->|<------------------------->|<--------->|<------------->|<--------->|<--------->|<------------------------->|<--------->|<--------------------->|         |<--------->|        |
%    |  mr.calcDuration  |     |           |      mr.calcDuation       |           |               |           |           |      mr.calcDuration      |           |             |         |         |           |        |
%    |      (gz_ex)      |     |           |         (gz_ref)          |           |               |           |           |         (gz_ref)          |           |         gx.riseTime   |         |           |        |
%    |                   |     |           |                           |           |               |           |           |                           |           |           ->| |<-     |       | |                    |
%    |                   |     |           |                           |           |               |           |           |                           |           |             | |_______|_______| |                    |
% Gx |                   |     |           |                           |           |               |           |           |                           |           |             | / |     |gx   | \ |                    |
% ___|___________________|_____|___________|___________________________|___________|_______________|___________|___________|___________________________|___________|____         |/| |     |     | |\|____________________|
%    |                   |     |           |                           |           |               |           |           |                           |           |    \        / | |     |     | | |                    |
%    |                   |     |           |                           |           |               |           |           |                           |           |     \______/| | |     |     | | |                    |
%    |                   |     |           |                           |           |               |           |           |                           |           | gx_prephaser| |<------------->| |<------------------>|
%    |                   |     |           |                           |           |               |           |           |                           |           |             | |   flat_time   | |      delayTR       |
%    |                   |     |           |                           |           |               |           |           |                           |           |             | | |xxxxxxxxxxx| | |                    |
%    |                   |     |           |                           |           |               |           |           |                           |   delay3  |<----------->| | |    ADC      | |                    |
%    |                   |     |           |                           |           |               |           |           |                           |           |   |shift_adc->| |<-   |       | |                    |
%    |                   |     |           |                           |           |               |           |           |                           |           |  _______    |                   |                    |
% Gy |                   |     |           |                           |           |               |           |           |                     gy_phase_encoding | /_______\   |                   |                    |
% ___|___________________|_____|___________|___________________________|___________|_______________|___________|___________|___________________________|___________|/_________\__|___________________|____________________|
%    |                   |     |           |                           |           |               |           |           |                           |           |\ _______ /  |                   |                    |
%    |                   |     |           |                           |           |               |           |           |                           |           | \_______/   |                   |                    |
%    |<----------------->|<--->|<--------->|<------------------------->|<--------->|<------------->|<--------->|<--------->|-------------------------->|<--------->|<----------->|<----------------->|<------------------>|
%    |      block 1      | b2  |  block 3  |          block 4          |  block 5  |     block 6   |   block 7 |  block 8  |          block 9          |  block 10 |   block 11  |      block 12     |      block 13      |
%----o-------------------+-----+-----------+---------------------------+-----------+---------------+-----------+-----------+---------------------------+-----------+-------------+-------------------+-------------------> t
%--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% tau    : tau = gz_ex.flatTime / 2 + gz_ex.fallTime + mr.calcDuration(gz_rephaser) + mr.calcDuration(left_crusher) + mr.calcDuration(rf_ref) / 2
% delay1 : tau = mr.calcDuration(rf_ref) / 2 + mr.calcDuration(right_crusher) + delay1
% delay2 : (TE - (2 * tau)) / 2 = delay2 + mr.calcDuration(left_crusher) + mr.calcDuration(rf_ref) / 2
% delay3 : (TE - (2 * tau)) / 2 =  mr.calcDuration(rf_ref) / 2 + mr.calcDuration(right_crusher) + delay3 + mr.calcDuration(gx)/ 2
% minTE  : minTE = tau * 2 + delay2 + mr.calcDuration(left_crusher) + mr.calcDuration(rf_ref) + mr.calcDuration(right_crusher) + min(mr.calcDuration(gy_phase_encoding), mr.calcDuration(gx_prephaser)) + mr.calcDuration(gx)/ 2
% delayTR: TR = (gz_ex.delay + gz_ex.riseTime + gz_ex.flatTime / 2 + TE + mr.calcDuration(gx) / 2 + delayTR) * nr_slices

%% Define parameters - original multi RF settings
load('.../multi_RF/figures/sim_results/phantom_055T/03312024_15cm/bloch_055T_1tx_offc15.0cm_iter20_dur18.840_dr0.25cm_rf_g_dt1e-6_fixcoord_mxyz.mat');
% rf : Nt by 1 in [mT]
% G  : Nt by 3 in [mT/m]
% dt : sampling rate [s]

% resample G based on sys.Gradrastertime
G_input = G(1:10:end,:); 

% check Gmax and slewrate
G_check = (G_input(:,1) + 1i*G_input(:,2)) * 1e1; % [G/cm]
[ktmp,gtmp,s,m1,m2,t,v] = calcgradinfo(G_check, sys.gradRasterTime);

%% Convert rf shape to [Hz] 
% [mT] * [Hz/T] * [T/1e3mT] -> [Hz]
multirf_ex = rf * sys.gamma * 1e-3;

%% Create an excitation RF event - multi RF
delay_rf = size(G_check,1) * sys.gradRasterTime - length(multirf_ex) * sys.rfRasterTime;
delay_rf = delay_rf + sys.gradRasterTime - 0e-6; % 
rf_sse2d = mr.makeArbitraryRf(multirf_ex, flip_angle * pi / 180, 'system', sys, 'delay', sys.rfDeadTime + delay_rf);
rf_sse2d.ringdownTime = rf_sse2d.ringdownTime + 0e-6;

%% Create spiral (readout + rewinder) gradient events ([PE,RO,SL] = [y,x,z] in Pulseq)
%--------------------------------------------------------------------------
% Create a spiral gradient event on the PE direction (PRS)
% [PE,RO,SL] = [y,x,z] in Pulseq
%--------------------------------------------------------------------------
% [mT/m] * [T/1e3mT] * [Hz/T] => [Hz/m]
% gy_delay  = 10e-6; % [s]
gy_spiral = mr.makeArbitraryGrad('y', G_input(:,1) * sys.gamma * 1e-3, 'system', sys, 'delay', sys.rfDeadTime + sys.gradRasterTime); % PE => 'y'

%--------------------------------------------------------------------------
% Create a spiral gradient event on the RO direction (PRS)
% [PE,RO,SL] = [y,x,z] in Pulseq
%--------------------------------------------------------------------------
% [mT/m] * [T/1e3mT] * [Hz/T] => [Hz/m]
gx_spiral = mr.makeArbitraryGrad('x', G_input(:,2) * sys.gamma * 1e-3, 'system', sys, 'delay', sys.rfDeadTime + sys.gradRasterTime); % RO => 'x'

%% 180 - Calculate a hyperbolic secant pulse
%--------------------------------------------------------------------------
% Calculate the modulation anglular frequency [rad/sec]
%--------------------------------------------------------------------------
beta = ref_bandwidth * pi / mu; % modulation anglular frequency [Hz] * [rad/cycle] => [rad/sec]

%--------------------------------------------------------------------------
% Calculate the maximum RF amplitude [uT]
%--------------------------------------------------------------------------
A0 = OMEGA0 / sys.gamma * 1e6; % [Hz] / [Hz/T] * [1e6uT/T] => * 1e6 [uT]

%--------------------------------------------------------------------------
% Calculate the number of RF samples
%--------------------------------------------------------------------------
rf_samples = ceil(ref_duration / sys.rfRasterTime);

%--------------------------------------------------------------------------
% Calculate a time axis for an RF pulse [sec]
%--------------------------------------------------------------------------
t = (-floor(rf_samples/2):ceil(rf_samples/2)-1).' * sys.rfRasterTime; % [sec]

%--------------------------------------------------------------------------
% Calculate the amplitude modulation A(t)
%--------------------------------------------------------------------------
am_shape = A0 * sech(beta * t); % [uT]

%--------------------------------------------------------------------------
% Calculate the phase modulation phi(t)
%--------------------------------------------------------------------------
pm_shape = mu * log(sech(beta * t)) + mu * log(A0);

%--------------------------------------------------------------------------
% Calculate a complex-valued hyperbolic secant pulse [uT]
%--------------------------------------------------------------------------
rf_shape = am_shape .* exp(1j * pm_shape); % [uT]

%% Create a hyperbolic secant pulse event [Hz] and the corresponding slice-selection gradient event [Hz/m]
grad_amplitude = ref_duration / slice_thickness; % [Hz] / [m] => [Hz/m]
% [Hz/T] * [T/1e6uT] * [uT] * [sec] * [2*pi rad/cycle] => [rad]
flip_angle_ref = (sys.gamma * 1e-6 * 2 * pi) * abs(sum(rf_shape)) * sys.rfRasterTime; % [rad]
[rf_ref, gz_ref, delay] = mr.makeArbitraryRf(rf_shape, flip_angle_ref, 'bandwidth', ref_bandwidth, 'sliceThickness', slice_thickness, 'system', sys);

%--------------------------------------------------------------------------
% Set the phase offset of an RF pulse
%--------------------------------------------------------------------------
rf_ref.phaseOffset = 0; % pi / 2; => no difference?

%% Create crusher gradient events for the refocusing pulse
%--------------------------------------------------------------------------
% phase = gamma * area * slice_thickness
% [cycle] = [Hz/T] * [T/1e3mT] * [mT/m] * [sec] * [m]
% (gamma * area) = phase / (slice_thickness)
% 'Area' = gamma * area: [Hz/T] * [T/1e3mT] * [mT/m] * [sec] => [cycle/m] 
%--------------------------------------------------------------------------
left_crusher  = mr.makeTrapezoid('z', 'Area', 4 / slice_thickness, 'system', sys); % 4 cycles over slice thickness
right_crusher = left_crusher;

%% Create a spoiler event on each axis
%--------------------------------------------------------------------------
% gradient amplitude in [Hz/m] = [Hz/T] * [mT/m]
% gradient area in [Hz/m] * [sec] => [Hz/m*sec] = [cycle/m]
% area = 4 / slice_thickness: 4 cycle dephasing across the slice thickness
%--------------------------------------------------------------------------
gx_spoiler = mr.makeTrapezoid('x', 'Area', 4 / slice_thickness, 'system', sys_derated); % 4 cycles over slice thickness
gy_spoiler = mr.makeTrapezoid('y', 'Area', 4 / slice_thickness, 'system', sys_derated); % 4 cycles over slice thickness
gz_spoiler = mr.makeTrapezoid('z', 'Area', 4 / slice_thickness, 'system', sys_derated); % 4 cycles over slice thickness

%% Calculate the duration of a trapezoid (TotalTime) [sec]
flat_time = ceil(adc_duration / sys.gradRasterTime) * sys.gradRasterTime; % [sec]

%% Create a readout gradient event ([PE,RO,SL] = [y,x,z] in Pulseq)
% Siemens GCS (PRS) 'RO' => Pulseq 'x'
deltak_read = 1 / (fov_read * readout_os_factor); % [cycle/m]
gx = mr.makeTrapezoid('x', 'FlatArea', adc_samples * deltak_read, 'FlatTime', flat_time, 'system', sys);
gx_prephaser = mr.makeTrapezoid('x', 'Area', -gx.area / 2, 'system', sys_derated);

%--------------------------------------------------------------------------
% Decrease the slew rate of a descending ramp
%--------------------------------------------------------------------------
gx_rise_time_derated = ceil(gx.amplitude / sys.gamma / max_slew_derated / sys.gradRasterTime) * sys.gradRasterTime;
gx_fall_time_derated = ceil(gx.amplitude / sys.gamma / max_slew_derated / sys.gradRasterTime) * sys.gradRasterTime;

gx.riseTime = gx_rise_time_derated;
gx.fallTime = gx_fall_time_derated;

%% Create a list of phase-encoding gradient areas ([PE,RO,SL] = [y,x,z] in Pulseq)
% Siemens GCS (PRS) 'PE' => Pulseq 'y'
deltak_phase = 1 / (fov_read * fov_phase * 1e-2); % [cycle/m]
phase_areas  = (-floor(nr_phase_encoding_steps_1/2):ceil(nr_phase_encoding_steps_1/2)-1).' * deltak_phase; % [cycle/m]

%% Flip the sign of gradients along the PE direction [PE,RO,SL]
% Note: for 3D encoding, SL needs to be flipped.
% This step is necessary to make use of coordinate transformations defined by
% Siemens on the datasets acquired with Pulseq.
phase_areas = flip(phase_areas,1); % PE

%% Create a phase-encoding gradient event ([PE,RO,SL] = [y,x,z] in Pulseq)
gy_phase_encoding = mr.makeTrapezoid('y', 'Area', max(abs(phase_areas)), 'system', sys_derated);

%% Create an ADC readout event
% NOT WORKING?? A BUG?? Here, adc_delay must be a multiple of 1 us (=1000 ns) instead of 100 ns.
% shift_adc = sys.adcDeadTime + round(((TE2 - TE1) - adc_duration) / 2 / sys.adcRasterTime) * sys.adcRasterTime; % [sec]
% shift_adc = round((flat_time - adc_duration) / 2 / (sys.adcRasterTime * 10)) * (sys.adcRasterTime * 10); % [sec]
% adc_delay = gx.riseTime + shift_adc;
adc_delay  = round((gx.riseTime + (flat_time - adc_duration) / 2) / sys.rfRasterTime) * sys.rfRasterTime;
adc        = mr.makeAdc(adc_samples, 'Dwell', real_dwell_time, 'delay', adc_delay, 'system', sys);

%% Calculate timing (need to decide on the block structure already)
%--------------------------------------------------------------------------
% tau    : tau = gz_ex.flatTime / 2 + gz_ex.fallTime + mr.calcDuration(gz_rephaser) + mr.calcDuration(left_crusher) + mr.calcDuration(rf_ref) / 2
% delay1 : tau = mr.calcDuration(rf_ref) / 2 + mr.calcDuration(right_crusher) + delay1
% delay2 : (TE - (2 * tau)) / 2 = delay2 + mr.calcDuration(left_crusher) + mr.calcDuration(rf_ref) / 2
% delay3 : (TE - (2 * tau)) / 2 = mr.calcDuration(rf_ref) / 2 + mr.calcDuration(right_crusher) + delay3 + mr.calcDuration(gx)/ 2
% minTE  : minTE = tau * 2 + delay2 + mr.calcDuration(left_crusher) + mr.calcDuration(rf_ref) + mr.calcDuration(right_crusher) + min(mr.calcDuration(gy_phase_encoding), mr.calcDuration(gx_prephaser)) + mr.calcDuration(gx)/ 2
% delayTR: TR = (gz_ex.delay + gz_ex.riseTime + gz_ex.flatTime / 2 + TE + mr.calcDuration(gx) / 2 + delayTR) * nr_slices
%--------------------------------------------------------------------------
tau = mr.calcDuration(left_crusher) + mr.calcDuration(rf_ref) / 2;
delay1 = tau - (mr.calcDuration(rf_ref) / 2 + mr.calcDuration(right_crusher));
delay2 = (TE - (2 * tau)) / 2 - (mr.calcDuration(left_crusher) + mr.calcDuration(rf_ref) / 2);
delay3 = (TE - (2 * tau)) / 2 - (mr.calcDuration(rf_ref) / 2 + mr.calcDuration(right_crusher) + mr.calcDuration(gx)/ 2);
minTE = tau * 2 + delay2 + mr.calcDuration(left_crusher) + mr.calcDuration(rf_ref) + mr.calcDuration(right_crusher) + min(mr.calcDuration(gy_phase_encoding), mr.calcDuration(gx_prephaser)) + mr.calcDuration(gx)/ 2;
delayTR = TR / nr_slices - (mr.calcDuration(rf_sse2d) + TE + mr.calcDuration(gx) / 2);

delay1  = round(delay1 / sys.gradRasterTime) * sys.gradRasterTime;
delay2  = round(delay2 / sys.gradRasterTime) * sys.gradRasterTime;
delay3  = round(delay3 / sys.gradRasterTime) * sys.gradRasterTime;
delayTR = round(delayTR / sys.gradRasterTime) * sys.gradRasterTime;

assert(delayTR > 0);

%% Calculate the slice ordering
if strcmp(series, 'Ascending')
    slice_ordering = (1:nr_slices).';
elseif strcmp(series, 'Descending')
    slice_ordering = (nr_slices:-1:1).';
elseif strcmp(series, 'Interleaved')
    slice_ordering = cat(1, (1:2:nr_slices).', (2:2:nr_slices).');
end

%% Calculate a list of frequency offsets for multi-slice imaging
% if mod(nr_slices,2) == 0 % even
%     freq_offset_list = ((-floor(nr_slices/2):ceil(nr_slices/2)-1).' + 0.5) * gz_ex.amplitude * slice_thickness * (1 + dist_factor * 1e-2); % [Hz/m] * [m] => [Hz]
% else % odd
%     freq_offset_list = ((-floor(nr_slices/2):ceil(nr_slices/2)-1).') * gz_ex.amplitude * slice_thickness * (1 + dist_factor * 1e-2); % [Hz/m] * [m] => [Hz]
% end
% freq_offset_list = flip(freq_offset_list,1);

%% Create a sequence object
seq = mr.Sequence(sys);
start_time = tic;

%% Define sequence blocks
% all LABELS / counters and flags are automatically initialized to 0 in the beginning, no need to define initial 0's  
% so we will just increment LIN after the ADC event (e.g. during the spoiler)
%--------------------------------------------------------------------------
% ISMRMRD header
% uint16_t kspace_encode_step_1;    /**< e.g. phase encoding line number */
% uint16_t kspace_encode_step_2;    /**< e.g. partition encoding number */
% uint16_t average;                 /**< e.g. signal average number */
% uint16_t slice;                   /**< e.g. imaging slice number */
% uint16_t contrast;                /**< e.g. echo number in multi-echo */
% uint16_t phase;                   /**< e.g. cardiac phase number */
% uint16_t repetition;              /**< e.g. dynamic number for dynamic scanning */
% uint16_t set;                     /**< e.g. flow encoding set */
% uint16_t segment;                 /**< e.g. segment number for segmented acquisition */
% uint16_t user[ISMRMRD_USER_INTS]; /**< Free user parameters */
%--------------------------------------------------------------------------
lbl_inc_lin   = mr.makeLabel('INC', 'LIN', 1); % lin == line
lbl_inc_par   = mr.makeLabel('INC', 'PAR', 1); % par == partition
lbl_inc_avg   = mr.makeLabel('INC', 'AVG', 1); % avg == average
lbl_inc_slc   = mr.makeLabel('INC', 'SLC', 1); % slc == slice

lbl_reset_lin = mr.makeLabel('SET', 'LIN', 0);
lbl_reset_par = mr.makeLabel('SET', 'PAR', 0);
lbl_reset_avg = mr.makeLabel('SET', 'AVG', 0);
lbl_reset_slc = mr.makeLabel('SET', 'SLC', 0);

%--------------------------------------------------------------------------
% Average (AVG)
%--------------------------------------------------------------------------
for avg = 1:nr_averages

    %----------------------------------------------------------------------
    % Update AVG counter
    %----------------------------------------------------------------------
    lbl_set_avg = mr.makeLabel('SET', 'AVG', avg - 1);

    %----------------------------------------------------------------------
    % Phase encoding line number (LIN)
    %----------------------------------------------------------------------
    for lin = 1:nr_phase_encoding_steps_1

        %------------------------------------------------------------------
        % Update LIN counter
        %------------------------------------------------------------------
        lbl_set_lin = mr.makeLabel('SET', 'LIN', lin - 1);

        %------------------------------------------------------------------
        % Define a phase encoding gradient event
        %------------------------------------------------------------------
        gy_phase_encoding = mr.makeTrapezoid('y', 'Area', phase_areas(lin), 'Duration', mr.calcDuration(gx_prephaser), 'system', sys_derated);

        %------------------------------------------------------------------
        % Slice (SLC)
        %------------------------------------------------------------------
        for slc = 1:nr_slices
            tstart = tic;
            fprintf('%s:(AVG=%2d/%2d)(LIN=%3d/%3d)(SLC=%2d/%2d): Defining blocks for SE-CARTESIAN data acquisitions... ', datetime, avg, nr_averages, lin, nr_phase_encoding_steps_1, slc, nr_slices);

            %--------------------------------------------------------------
            % Set the frequency offset of an RF pulse event [Hz]
            %--------------------------------------------------------------
            % rf_ex.freqOffset = freq_offset_list(slice_ordering(slc));

            %--------------------------------------------------------------
            % Update SLC counter
            %--------------------------------------------------------------
            lbl_set_slc = mr.makeLabel('SET', 'SLC', slc - 1);

            %--------------------------------------------------------------
            % Add a new block to the sequence (Block 1)
            %--------------------------------------------------------------
            rf_sse2d.freqOffset = 0;

            seq.addBlock(rf_sse2d, gx_spiral, gy_spiral, lbl_set_avg, lbl_set_lin, lbl_set_slc);

            %--------------------------------------------------------------
            % Add a new block to the sequence (Block 2)
            %--------------------------------------------------------------
%           % seq.addBlock(mr.scaleGrad(gz_rephaser,0));

            %--------------------------------------------------------------
            % Add a new block to the sequence (Block 3)
            %--------------------------------------------------------------
            seq.addBlock(left_crusher);

            %--------------------------------------------------------------
            % Add a new block to the sequence (Block 4)
            %--------------------------------------------------------------
            seq.addBlock(rf_ref, gz_ref);

            %--------------------------------------------------------------
            % Add a new block to the sequence (Block 5)
            %--------------------------------------------------------------
            seq.addBlock(right_crusher);

            %--------------------------------------------------------------
            % Add a new block to the sequence (Block 6)
            %--------------------------------------------------------------
            % seq.addBlock(mr.makeDelay(delay1));

            %--------------------------------------------------------------
            % Add a new block to the sequence (Block 7)
            %--------------------------------------------------------------
            seq.addBlock(mr.makeDelay(delay2));

            %--------------------------------------------------------------
            % Add a new block to the sequence (Block 8)
            %--------------------------------------------------------------
            seq.addBlock(left_crusher);

            %--------------------------------------------------------------
            % Add a new block to the sequence (Block 9)
            %--------------------------------------------------------------
            seq.addBlock(rf_ref, gz_ref);

            %--------------------------------------------------------------
            % Add a new block to the sequence (Block 10)
            %--------------------------------------------------------------
            seq.addBlock(right_crusher);

            %--------------------------------------------------------------
            % Add a new block to the sequence (Block 11)
            %--------------------------------------------------------------
            seq.addBlock(mr.align('left', mr.makeDelay(delay3), gy_phase_encoding, 'right', gx_prephaser));

            %--------------------------------------------------------------
            % Add a new block to the sequence (Block 12)
            %--------------------------------------------------------------
            seq.addBlock(gx, adc);

            %--------------------------------------------------------------
            % Add a new block to the sequence (Block 13)
            %--------------------------------------------------------------
            seq.addBlock(mr.makeDelay(delayTR), gx_spoiler, gz_spoiler, mr.scaleGrad(gy_phase_encoding,-1));
            fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
        end % SLC
    end % LIN
end % AVG

%% check whether the timing of the sequence is correct
[ok, error_report] = seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

%% prepare sequence export
seq.setDefinition('Averages', nr_averages);
seq.setDefinition('BaseResolution', base_resolution);
seq.setDefinition('FOV', [fov_read fov_read * fov_phase * 1e-2 slice_thickness]);
seq.setDefinition('MaxGrad', max_grad);
seq.setDefinition('MaxSlew', max_slew);
seq.setDefinition('Name', 'Multislice Cartesian SE');
seq.setDefinition('PhaseEncodingSteps1', nr_phase_encoding_steps_1);
seq.setDefinition('ReadoutOSFactor', readout_os_factor);
seq.setDefinition('RealDwellTime', real_dwell_time);
seq.setDefinition('Slices', nr_slices);
seq.setDefinition('TE', TE);
seq.setDefinition('TR', TR);
seq.setDefinition('DistFactor', dist_factor);

seq_filename = sprintf('03312024_sse2d_se_multirf_ori_off15_%3.0fx%3.0fmm_%dx%d_%gmm_DF%d_TE%d_TR%2.0f_%3.0fHz_slc%d_avg%d_osf%d_fov12_8_y1x2_nodelayrf.seq', fov_read * 1e3, fov_read * fov_phase * 1e-2 * 1e3, base_resolution, nr_phase_encoding_steps_1, slice_thickness * 1e3, dist_factor, TE * 1e3, TR * 1e3, bandwidth, nr_slices, nr_averages, readout_os_factor);
seq_file = fullfile(output_directory, seq_filename);
seq.write(seq_file);   % Output sequence for scanner

%% plot sequence and k-space diagrams
seq.plot('timeRange', [0 1] * TR, 'label', 'LIN');

%% very optional slow step, but useful for testing during development e.g. for the real TE, TR or for staying within slewrate limits  
if 0
    rep = seq.testReport;
    fprintf([rep{:}]);
end