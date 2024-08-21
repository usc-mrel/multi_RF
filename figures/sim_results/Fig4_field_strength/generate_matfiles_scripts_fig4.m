%% Figure 4 script to generate .mat files 
% sweep through different off-isocenters

clear all; close all; clc;

cur_dir = pwd;
cd('../../../sim');
setup_path;

output_path = fullfile(root_path, '/multi_RF/figures/sim_results/Fig4_field_strength/');

B0_list = [0.2 0.55 1.5 3 7];  % [tesla]
zoff = [10]*1e-2;    % [m]
T2  = 1e6;
fov = 12.8;

for i = 1:length(B0_list)
    
    B0 = B0_list(i);

    % generate blochmex .mat files
    concomitant_correct = 1; 
    sim_concomitantfields_8channel;
    
    % generate bloch .mat files
    concomitant_correct = 0; 
    sim_concomitantfields_8channel;

end

