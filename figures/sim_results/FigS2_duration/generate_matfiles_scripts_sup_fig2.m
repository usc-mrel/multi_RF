%% Supporting Figure S2 script to generate .mat files
% sweep through different off-isocenters

clear all; close all; clc;

cur_dir = pwd;
cd('../../../sim');
setup_path;
output_path = fullfile(root_path, '/multi_RF/figures/sim_results/FigS2_duration/');

paths.output_path = output_path;
paths.root_path = root_path;
display_figure  = 0;

B0_list = [0.55];         % [tesla]
zoff_list = [10]*1e-2;    % [m]
T2_list = [1e6 40e-3];    % [s]
fov_list = [12.8 12.8*3]; % [cm]

zoff = zoff_list(1);
B0 = B0_list(1);

for i = 1 : length(T2_list)
    for j = 1 : length(fov_list)
    
        T2 = T2_list(i);
        fov = fov_list(j);
    
        % generate blochmex .mat files
        concomitant_correct = 1; 
        sim_concomitantfields_8channels(B0, zoff, T2, fov, concomitant_correct, display_figure, paths);
    
        % generate bloch .mat files
        concomitant_correct = 0; 
        sim_concomitantfields_8channels(B0, zoff, T2, fov, concomitant_correct, display_figure, paths);
    
    end
end
