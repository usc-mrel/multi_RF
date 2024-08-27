%% Figure 3 script to generate .mat files 
% sweep through different off-isocenters

clear all; close all; clc;

cur_dir = pwd;
cd('../../../sim');
setup_path;

output_path = fullfile(root_path, '/multi_RF/figures/sim_results/Fig3_distance/');
paths.output_path = output_path;
paths.root_path = root_path;

display_figure = 0; % uncomment to generate Supporting Figure S1

B0_list = [0.55];         % [tesla]
zoff_list = [0, 05, 10, 15, 20, 25, 30]*1e-2;    % [m]
T2 = 1e6;
fov = 12.8;


for i = 1 : length(zoff_list)
    
    zoff = zoff_list(i);
    B0 = B0_list(1);

    % generate blochmex .mat files
    concomitant_correct = 1; 
    sim_concomitantfields_8channels(B0, zoff, T2, fov, concomitant_correct, display_figure, paths)

    % generate bloch .mat files
    concomitant_correct = 0; 
    sim_concomitantfields_8channels(B0, zoff, T2, fov, concomitant_correct, display_figure, paths)

end
