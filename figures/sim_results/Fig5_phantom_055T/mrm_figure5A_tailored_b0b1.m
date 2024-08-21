% Figure 5A generating figures from simulation results at 0.55T  
% 
% Ziwei Zhao 
% 05122024
clear all;
close all;
clc;

cur_dir = pwd;
cd('../../../sim');
setup_path;

cd(fullfile(root_path, '/multi_RF/figures/sim_results/Fig5_phantom_055T/'));

folder = dir('./*cm');
cur_dir = pwd;
ori_mxyz = zeros(64, 64, 3, length(folder));
pro_mxyz = zeros(64, 64, 3, length(folder));

for i_folder = 1:length(folder)
    
    cd(folder(i_folder).name);
    
    % load results
    ori_results = dir('bloch_*_mxyz.mat');
    pro_results = dir('blochmex_*_mxyz.mat');

    load(ori_results(1).name);
    ori_mxyz(:,:,:,i_folder) = mxyz_offcenter;
    load(pro_results(1).name);
    pro_mxyz(:,:,:,i_folder) = mxyz_offcenter;

    cd(cur_dir);
end

%% generate
% load re-verse method results
flip = pi/2;
mxy_ori = abs(complex(ori_mxyz(:,:,1,:), ori_mxyz(:,:,2,:)));
mxy_pro = abs(complex(pro_mxyz(:,:,1,:), pro_mxyz(:,:,2,:)));

mxy_ori = mxy_ori(9:64-8, 9:64-8, :, :);
mxy_pro = mxy_pro(9:64-8, 9:64-8, :, :);

N1 = 48; N2 = 48;
block = 1.01 * ones(N1,1, 'double');

tmp_ori = cat(2, mxy_ori(:,:,1,1).', block, mxy_ori(:,:,1,2).', block, mxy_ori(:,:,1,3).', block, mxy_ori(:,:,1,4).');
tmp_pro = cat(2, mxy_pro(:,:,1,1).', block, mxy_pro(:,:,1,2).', block, mxy_pro(:,:,1,3).', block, mxy_pro(:,:,1,4).');

im_montage = cat(1, tmp_ori, 1.01 * ones(1, (N1+1)*4-1, 'double'), tmp_pro);                
                
FontSize = 26;
cmap = cat(1, jet(256), [1 1 1]);                
figure('Color', 'w', 'Position', [-1 2 1101 814]);
imagesc(abs(im_montage)); axis image off; colormap(gca, cmap); 
hc = colorbar; %clim([0 1]);

set(hc,'FontSize', FontSize);
ylabel(hc, 'Scaled Magnetization [%]'); % ^\circ

text(N2/2+(N2+1)*1.5  , -10, 'B_0 = 0.55T - 1 channel simulation', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', FontSize+2, 'FontWeight','bold');
text(N2/2         , 0, {sprintf('z = %dcm', 0)},  'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', FontSize);
text(N2/2+(N2+1)  , 0, {sprintf('z = %dcm', 5)}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', FontSize);
text(N2/2+(N2+1)*2, 0, {sprintf('z = %dcm', 10)}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', FontSize);
text(N2/2+(N2+1)*3, 0, {sprintf('z = %dcm', 15)}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', FontSize);

%% NRMSE 
NRMSE_ori = [0.0581, 0.1908, 0.2995, 0.3821];
NRMSE_pro = [0.0581, 0.0591, 0.0767, 0.1218];

figure;  plot(NRMSE_ori, 'LineWidth', 3, 'Color', [0.4940 0.1840 0.5560]);
hold on; plot(NRMSE_pro, 'LineWidth', 3, 'Color', [0.8500 0.3250 0.0980]);
set(gca, 'FontSize', 26); grid on;
xticklabels({'0cm', '5cm', '10cm', '15cm', '20cm', '25cm', '30cm'});
xlim([0.8 4.2]);
ylim([0 0.5]);
box off;
ax = gca;
ax.LineWidth = 2;
title('NRMSE'); 
legend('Original', 'Proposed');

