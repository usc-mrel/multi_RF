% figure 5 duration impacts
% re-verse method vs. proposed max-bloch
% ziwei 02192024

% load re-verse method results
% size: [Nx Ny Nz dur]
% B0:  0.55T
% off: 10cm
clear all;
close all;
clc;

%% load exciatation patterns for reVERSE method
reverse_girf = zeros(80,80,3,4);
index = [1 3 2 4];
tmp = dir('/Users/ziwei/Documents/matlab/STA_maxwell/sim_code_MRM/sim_results/duration/bloch_*dr0.25cm*.mat');
for i = 1:length(tmp)
    load(tmp(index(i)).name);
    reverse_girf(:,:,:,i) = mxyz_offcenter;
end

% load proposed method results
reverse_max = zeros(80,80,3,4);
tmp = dir('/Users/ziwei/Documents/matlab/STA_maxwell/sim_code_MRM/sim_results/duration/blochmex_*dr0.25cm*.mat');
for i = 1:length(tmp)
    load(tmp(index(i)).name);
    reverse_max(:,:,:,i) = mxyz_offcenter;
end

%% parameters 
load /Users/ziwei/Documents/matlab/STA_maxwell/spiral2d_rf_pulse_design/phase_relaxed_CPMG_excitation/b0_b1_maps_2d.mat;

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

%% duration depenence figure
mxy_reverse  = complex(reverse_girf(:,:,1,:), reverse_girf(:,:,2,:));
mxy_maxverse = complex(reverse_max(:,:,1,:), reverse_max(:,:,2,:));

block = 1.01 * complex(ones(N1,1, 'double'), ones(N1,1, 'double'));
im_montage = cat(1, imag(cat(2, mxy_reverse(:,:,1,1).', block, mxy_reverse(:,:,1,2).', block, mxy_reverse(:,:,1,3).', block, mxy_reverse(:,:,1,4).')), ...
                    1.01 * ones(1, (N2+1)*4-1, 'double'), ...
                    real(cat(2, mxy_reverse(:,:,1,1).', block, mxy_reverse(:,:,1,2).', block, mxy_reverse(:,:,1,3).', block, mxy_reverse(:,:,1,4).')), ...
                    1.01 * ones(1, (N2+1)*4-1, 'double'), ...
                    imag(cat(2, mxy_maxverse(:,:,1,1).', block, mxy_maxverse(:,:,1,2).', block, mxy_maxverse(:,:,1,3).', block, mxy_maxverse(:,:,1,4).')), ...
                    1.01 * ones(1, (N2+1)*4-1, 'double'), ...
                    real(cat(2, mxy_maxverse(:,:,1,1).', block, mxy_maxverse(:,:,1,2).', block, mxy_maxverse(:,:,1,3).', block, mxy_maxverse(:,:,1,4).')));

FontSize = 18;
cmap = cat(1, jet(256), [1 1 1]);                
figure('Color', 'w', 'Position', [-1 2 1101 814]);
imagesc(abs(im_montage)); axis image off; colormap(gca, cmap); 
hc = colorbar;

set(hc,'FontSize', FontSize);
ylabel(hc, 'Scaled Magnetization [%]');

text((N2+1)*2 , -15, 'B_0 = 0.55T, z = 10cm - 8 channel simulation', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', FontSize+2, 'FontWeight','bold');

text(N2/2         , 0, {'17.66ms'}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', FontSize);
text(N2/2+(N2+1)  , 0, {'52.07ms'}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', FontSize);
text(N2/2+(N2+1)*2, 0, {'17.66ms'},  'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', FontSize);
text(N2/2+(N2+1)*3, 0, {'52.07ms'},  'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', FontSize);

%% NRMSE 
for i = 1:4 % field strengths 
        mxy_ori = mxy_reverse(:,:,i);
        mxy_pro = mxy_maxverse(:,:,i);
%         NRMSE_ori_real(i) = sqrt(sum(sum((real(mxy_ori) - real(P_)).^2))) / sqrt(sum(sum(abs(P_).^2)));
%         NRMSE_ori_imag(i) = sqrt(sum(sum((imag(mxy_ori) - imag(P_)).^2))) / sqrt(sum(sum(abs(P_).^2)));
%         NRMSE_pro_real(i) = sqrt(sum(sum((real(mxy_pro) - real(P_)).^2))) / sqrt(sum(sum(abs(P_).^2)));
%         NRMSE_pro_imag(i) = sqrt(sum(sum((imag(mxy_pro) - imag(P_)).^2))) / sqrt(sum(sum(abs(P_).^2)));
        NRMSE_ori(i) = sqrt(sum(sum((abs(mxy_ori) - abs(P_)).^2))) / sqrt(sum(sum(abs(P_).^2)));
        NRMSE_pro(i) = sqrt(sum(sum((abs(mxy_pro) - abs(P_)).^2))) / sqrt(sum(sum(abs(P_).^2)));
end

figure;  plot([1, 2], NRMSE_ori(1:2), 'LineWidth', 3, 'Color', [0.4940 0.1840 0.5560]);
hold on; plot([1, 2], NRMSE_pro(1:2), 'LineWidth', 3, 'Color', [0.8500 0.3250 0.0980]);
hold on; plot([1, 2], NRMSE_ori(3:4), 'LineWidth', 3, 'Color', [0.6350 0.0780 0.1840]);
hold on; plot([1, 2], NRMSE_pro(3:4), 'LineWidth', 3, 'Color', [0.9290 0.6940 0.1250]);

set(gca, 'FontSize', 20); grid on;
xticks([1 2]);
xticklabels({'17.66ms', '52.07ms'});
xlim([0.8 2.2]);
ylim([0.1 0.5]); box off; 
ax = gca;
ax.LineWidth = 2;
title('NRMSE'); legend('Original    - T2=inf', 'Proposed - T2=inf', 'Original    - T2=40ms', 'Proposed - T2=40ms');



    