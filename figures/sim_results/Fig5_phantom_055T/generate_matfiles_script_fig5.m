%% Figure 5A script to generate .mat files 
% sweep through different off-isocenters

clear all; close all; clc;

cur_dir = pwd;
cd('../../../sim');
setup_path;

file_path = fullfile(root_path, '/multi_RF/figures/sim_results/Fig5_phantom_055T/');
B0_map_path = fullfile(root_path, '/multi_RF/b0b1_map_055T/');

%% 0cm isocenter
zoff = 0;
% load related measurements
load(fullfile(B0_map_path, '/F0_64/fieldmap_b0_fov220mm_matrix64_02282024.mat')); % B0 field map
load(fullfile(B0_map_path, '/F0_64/mask_F0.mat')); % mask
m = imresize(m,[64,64], Method='nearest'); % downsample mask
dicom_path = fullfile(B0_map_path, '/F0_64'); 
output_path = [file_path, '03312024_00cm'];

% generate blochmex .mat files
concomitant_correct = 1; 
sim_concomitantfields_1channel;
concomitant_correct = 0; 
sim_concomitantfields_1channel;

%% 5cm isocenter
zoff = 5e-2;        % [m]
% load related measurements
load(fullfile(B0_map_path, '/F5_64/fieldmap_b0_fov220mm_matrix64_off_5cm_03312024.mat')); % B0 field map
load(fullfile(B0_map_path, '/F5_64/mask_F5.mat')); % mask
dicom_path = fullfile(B0_map_path, '/F5_64'); 
output_path = [file_path, '03312024_05cm'];

% generate blochmex .mat files
concomitant_correct = 1; 
sim_concomitantfields_1channel;
concomitant_correct = 0; 
sim_concomitantfields_1channel;

%% 10cm isocenter
zoff = 10e-2;        % [m]
% load related measurements
load(fullfile(B0_map_path, '/F10_64/fieldmap_b0_fov220mm_matrix64_off_10cm_03302024.mat')); % B0 field map
load(fullfile(B0_map_path, '/F10_64/mask_F10.mat')); % mask
dicom_path = fullfile(B0_map_path, '/F10_64'); 
output_path = [file_path, '03312024_10cm'];

% generate blochmex .mat files
concomitant_correct = 1; 
sim_concomitantfields_1channel;
concomitant_correct = 0; 
sim_concomitantfields_1channel;

%% 15cm isocenter
zoff = 15e-2;        % [m]
% load related measurements
load(fullfile(B0_map_path, '/F15_64/fieldmap_b0_fov220mm_matrix64_off_15cm_03312024.mat')); % B0 field map
load(fullfile(B0_map_path, '/F15_64/mask_F15.mat')); % mask
dicom_path = fullfile(B0_map_path, '/F15_64'); 
output_path = [file_path, '03312024_15cm'];

% generate blochmex .mat files
concomitant_correct = 1; 
sim_concomitantfields_1channel;
concomitant_correct = 0; 
sim_concomitantfields_1channel;
