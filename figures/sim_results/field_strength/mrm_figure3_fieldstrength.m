% supporting figure 1 - isocenter distance
% re-verse method vs. proposed max-bloch
% ziwei 02192024

% load re-verse method results
% size: [Nx Ny Nz dur]
% B0:  0.55T
% off: 10cm
clear all;
close all;
clc;

addpath(genpath('.../multi_RF/third_party/'));

%% load exciatation patterns for reVERSE method
reverse_girf = zeros(80,80,3,5);
tmp = dir('.../multi_RF/figures/sim_results/field_strength/bloch_B0*.mat');
for i = 1:length(tmp)
    load(tmp(i).name);
    reverse_girf(:,:,:,i) = mxyz_offcenter;
end

% load proposed method results
reverse_max = zeros(80,80,3,5);
tmp = dir('.../multi_RF/figures/sim_results/field_strength/blochmex_B0*.mat');
for i = 1:length(tmp)
    load(tmp(i).name);
    reverse_max(:,:,:,i) = mxyz_offcenter;
end

%% parameters 
load('.../multi_RF/third_party/phase_relaxed_CPMG_excitation/b0_b1_maps_2d.mat');

[N1,N2,Nc] = size(tx);

xoff = +6e-3;   % [mm] -> [m]
yoff = 0e-3;    % [mm] -> [m]
r0   = 15e-3;   % [mm] -> [m]
r0   = 1.74 * r0;

% Square beam
P = double(((abs(X - xoff) < r0) & (abs(Y - yoff) < r0)));

% apodize 
h = fspecial('gaussian', 3, 0.5);
P_ = imfilter(P,h);

flip = 90;                 % total flip angle [degrees]
flip = flip * pi / 180;

%% B0 dependence figure
% zoff = 10cm
mxy_reverse  = squeeze(complex(reverse_girf(:,:,1,:), reverse_girf(:,:,2,:)));
mxy_maxverse = squeeze(complex(reverse_max(:,:,1,:), reverse_max(:,:,2,:)));

block = 1.01 * complex(ones(N1,1, 'double'), ones(N1,1, 'double'));
im_montage = cat(1, imag(cat(2, 1j*P_.', block, mxy_reverse(:,:,1).',  block, mxy_reverse(:,:,2).',  block, mxy_reverse(:,:,3).',  block, mxy_reverse(:,:,4).', block, mxy_reverse(:,:,5).')), ...
                    1.01 * ones(1, (N2+1)*6-1, 'double'), ...
                    real(cat(2, 1j*P_.', block, mxy_reverse(:,:,1).',  block, mxy_reverse(:,:,2).',  block, mxy_reverse(:,:,3).',  block, mxy_reverse(:,:,4).', block, mxy_reverse(:,:,5).')), ...
                    1.01 * ones(1, (N2+1)*6-1, 'double'), ...
                    imag(cat(2, 1j*P_.', block, mxy_maxverse(:,:,1).', block, mxy_maxverse(:,:,2).', block, mxy_maxverse(:,:,3).', block, mxy_maxverse(:,:,4).', block, mxy_maxverse(:,:,5).')), ...
                    1.01 * ones(1, (N2+1)*6-1, 'double'), ...
                    real(cat(2, 1j*P_.', block, mxy_maxverse(:,:,1).', block, mxy_maxverse(:,:,2).', block, mxy_maxverse(:,:,3).', block, mxy_maxverse(:,:,4).', block, mxy_maxverse(:,:,5).')));

FontSize = 18;
cmap = cat(1, jet(256), [1 1 1]);                
figure('Color', 'w', 'Position', [-1 2 1101 814]);
imagesc(abs(im_montage)); axis image off; colormap(gca, cmap); 
hc = colorbar;

set(hc,'FontSize', FontSize);
ylabel(hc, 'Scaled Magnetization [%]');

text((N2+1)*3 , -15, 'z = 10mm', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', FontSize+2, 'FontWeight','bold');

text(N2/2         , 0, 'Target', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', FontSize);
text(N2/2+(N2+1)  , 0, {sprintf('%gT', 0.2)},  'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', FontSize);
text(N2/2+(N2+1)*2, 0, {sprintf('%gT', 0.55)}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', FontSize);
text(N2/2+(N2+1)*3, 0, {sprintf('%gT', 1.5)},  'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', FontSize);
text(N2/2+(N2+1)*4, 0, {sprintf('%gT', 3)},    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', FontSize);
text(N2/2+(N2+1)*5, 0, {sprintf('%gT', 7)},    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', FontSize);

% Draw arrows
cx = floor(N1/2) + 1;
cy = floor(N2/2) + 1;
arrow([cx cy], [cx+20 cy], 'Color', 'w', 'Length', 7, 'TipAngle', 25, 'Width', 1);
arrow([cx cy], [cx cy-20], 'Color', 'w', 'Length', 7, 'TipAngle', 25, 'Width', 1);
text(cx+20, cy+1 , '$\bf{x}$', 'Color', 'w', 'Interpreter', 'latex', 'FontSize', FontSize, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
text(cx+5 , cy-16, '$\bf{y}$', 'Color', 'w', 'Interpreter', 'latex', 'FontSize', FontSize, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

cx = floor(N1/2) + 1;
cy = floor(N2/2) + 1;
arrow([cx cy+(N1+1)*2], [cx+20 cy+(N1+1)*2], 'Color', 'w', 'Length', 7, 'TipAngle', 25, 'Width', 1);
arrow([cx cy+(N1+1)*2], [cx cy-20+(N1+1)*2], 'Color', 'w', 'Length', 7, 'TipAngle', 25, 'Width', 1);
text(cx+20, cy+1+(N1+1)*2 , '$\bf{x}$', 'Color', 'w', 'Interpreter', 'latex', 'FontSize', FontSize, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
text(cx+5 , cy-16+(N1+1)*2, '$\bf{y}$', 'Color', 'w', 'Interpreter', 'latex', 'FontSize', FontSize, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

%% NRMSE plot 
% re-verse method 
for i = 1:5 % field strengths 
        mxy_ori = mxy_reverse(:,:,i);
        mxy_pro = mxy_maxverse(:,:,i);
        NRMSE_ori(i) = sqrt(sum(sum((abs(mxy_ori) - abs(P_)).^2))) / sqrt(sum(sum(abs(P_).^2)));
        NRMSE_pro(i) = sqrt(sum(sum((abs(mxy_pro) - abs(P_)).^2))) / sqrt(sum(sum(abs(P_).^2)));
end

figure;  plot(NRMSE_ori, 'LineWidth', 3, 'Color', [0.4940 0.1840 0.5560]);
hold on; plot(NRMSE_pro, 'LineWidth', 3, 'Color', [0.8500 0.3250 0.0980]);
set(gca, 'FontSize', 20); grid on;
xticklabels({'0.2T', '0.55T', '1.5T', '3.0T', '7.0T'});
xlim([0.8 5]);
ylim([0 0.5]); box off; 
ax = gca;
ax.LineWidth = 2;
title('NRMSE'); legend('Original', 'Proposed');
    