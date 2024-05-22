% Supporting figure 1 - isocenter distance
% re-verse method vs. proposed max-bloch
% Ziwei 02192024

% load re-verse method results
% size: [Nx Ny Nz dur]
% B0:  0.55T
% off: 10cm
clear all;
close all;
clc;

%% load exciatation patterns for reVERSE method
reverse_girf = zeros(80,80,3,7);
tmp = dir('/Users/ziwei/Documents/matlab/STA_maxwell/sim_code_MRM/sim_results/distance/bloch_B00.55_offc*.mat');
for i = 1:length(tmp)
    load(tmp(i).name);
    reverse_girf(:,:,:,i) = mxyz_offcenter;
end

% load proposed method results
reverse_max = zeros(80,80,3,7);
tmp = dir('/Users/ziwei/Documents/matlab/STA_maxwell/sim_code_MRM/sim_results/distance/blochmex_B00.55_offc*.mat');
for i = 1:length(tmp)
    load(tmp(i).name);
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
P    = double(((abs(X - xoff) < r0) & (abs(Y - yoff) < r0)));

% apodize 
h    = fspecial('gaussian', 3, 0.5);
P_   = imfilter(P,h);

flip = 90;                 % total flip angle [degrees]
flip = flip * pi / 180;

%% B0 dependence figure
% zoff = 10cm
% mxy_reverse  = squeeze(complex(reverse_girf(:,:,1,:), reverse_girf(:,:,2,:)));
% mxy_maxverse = squeeze(complex(reverse_max(:,:,1,:), reverse_max(:,:,2,:)));
% 
% block = 1.01 * complex(ones(N1,1, 'double'), ones(N1,1, 'double'));
% im_montage = cat(1, imag(cat(2, 1j*P_.', block, mxy_reverse(:,:,1).', block, mxy_reverse(:,:,2).', block, mxy_reverse(:,:,3).')), ...
%                     1.01 * ones(1, (N2+1)*4-1, 'double'), ...
%                     real(cat(2, 1j*P_.', block, mxy_reverse(:,:,1).', block, mxy_reverse(:,:,2).', block, mxy_reverse(:,:,3).')), ...
%                     1.01 * ones(1, (N2+1)*4-1, 'double'), ...
%                     imag(cat(2, 1j*P_.', block, mxy_maxverse(:,:,1).', block, mxy_maxverse(:,:,2).', block, mxy_maxverse(:,:,3).')), ...
%                     1.01 * ones(1, (N2+1)*4-1, 'double'), ...
%                     real(cat(2, 1j*P_.', block, mxy_maxverse(:,:,1).', block, mxy_maxverse(:,:,2).', block, mxy_maxverse(:,:,3).')));
% 
% FontSize = 18;
% cmap = cat(1, jet(256), [1 1 1]);                
% figure('Color', 'w', 'Position', [-1 2 1101 814]);
% imagesc(abs(im_montage)*flip*180/pi); axis image off; colormap(gca, cmap); 
% hc = colorbar;
% 
% set(hc,'FontSize', FontSize);
% ylabel(hc, 'Flip angle [degree]');
% 
% text((N2+1)*2 , -15, 'B0 = 0.55T z = 10mm', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', FontSize+2, 'FontWeight','bold');
% 
% text(N2/2         , 0, 'Target', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', FontSize);
% text(N2/2         , 0, sprintf('re-VERSE'), 'Color', 'yellow', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', 'FontSize', FontSize);
% text(N2/2+(N2+1)  , 0, {sprintf('X%g', 4)}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', FontSize);
% text(N2/2+(N2+1)*2, 0, {sprintf('X%g', 2)},  'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', FontSize);
% text(N2/2+(N2+1)*3, 0, {sprintf('X%g', 1)},    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', FontSize);
% 
% text(N2/2, (N1+1)*2, sprintf('Proposed'), 'Color', 'yellow', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', 'FontSize', FontSize);
% text(-3, N2/2          , 'Imaginary', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', FontSize, 'rotation', 90);
% text(-3, N2/2+N2+1    ,  'Real'     , 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', FontSize, 'rotation', 90);
% text(-3, N2/2+(N2+1)*2,  'Imaginary', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', FontSize, 'rotation', 90);
% text(-3, N2/2+(N2+1)*3,  'Real'     , 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', FontSize, 'rotation', 90);
% 
% % Draw arrows
% cx = floor(N1/2) + 1;
% cy = floor(N2/2) + 1;
% %arrow([cx cy], [cx+20 cy], 'Color', 'w', 'Length', 7, 'TipAngle', 25, 'Width', 1);
% %arrow([cx cy], [cx cy-20], 'Color', 'w', 'Length', 7, 'TipAngle', 25, 'Width', 1);
% text(cx+20, cy+1 , '$\bf{x}$', 'Color', 'w', 'Interpreter', 'latex', 'FontSize', FontSize, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
% text(cx+5 , cy-16, '$\bf{y}$', 'Color', 'w', 'Interpreter', 'latex', 'FontSize', FontSize, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
% 
% cx = floor(N1/2) + 1;
% cy = floor(N2/2) + 1;
% %arrow([cx cy+(N1+1)*2], [cx+20 cy+(N1+1)*2], 'Color', 'w', 'Length', 7, 'TipAngle', 25, 'Width', 1);
% %arrow([cx cy+(N1+1)*2], [cx cy-20+(N1+1)*2], 'Color', 'w', 'Length', 7, 'TipAngle', 25, 'Width', 1);
% text(cx+20, cy+1+(N1+1)*2 , '$\bf{x}$', 'Color', 'w', 'Interpreter', 'latex', 'FontSize', FontSize, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
% text(cx+5 , cy-16+(N1+1)*2, '$\bf{y}$', 'Color', 'w', 'Interpreter', 'latex', 'FontSize', FontSize, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

%% off-isocenter figure
mxy_reverse  = complex(reverse_girf(:,:,1,:), reverse_girf(:,:,2,:));
mxy_maxverse = complex(reverse_max(:,:,1,:), reverse_max(:,:,2,:));

block = 1.01 * complex(ones(N1,1, 'double'), ones(N1,1, 'double'));
im_montage = cat(1, imag(cat(2, mxy_reverse(:,:,1,1).', block, mxy_reverse(:,:,1,2).', block, mxy_reverse(:,:,1,3).', block, mxy_reverse(:,:,1,4).', block, mxy_reverse(:,:,1,5).', block, mxy_reverse(:,:,1,6).', block, mxy_reverse(:,:,1,7).')), ...
                    1.01 * ones(1, (N2+1)*7-1, 'double'), ...
                    real(cat(2, mxy_reverse(:,:,1,1).', block, mxy_reverse(:,:,1,2).', block, mxy_reverse(:,:,1,3).', block, mxy_reverse(:,:,1,4).', block, mxy_reverse(:,:,1,5).', block, mxy_reverse(:,:,1,6).', block, mxy_reverse(:,:,1,7).')), ...
                    1.01 * ones(1, (N2+1)*7-1, 'double'), ...
                    imag(cat(2, mxy_maxverse(:,:,1,1).', block, mxy_maxverse(:,:,1,2).', block, mxy_maxverse(:,:,1,3).', block, mxy_maxverse(:,:,1,4).', block, mxy_reverse(:,:,1,5).', block, mxy_reverse(:,:,1,6).', block, mxy_reverse(:,:,1,7).')), ...
                    1.01 * ones(1, (N2+1)*7-1, 'double'), ...
                    real(cat(2, mxy_maxverse(:,:,1,1).', block, mxy_maxverse(:,:,1,2).', block, mxy_maxverse(:,:,1,3).', block, mxy_maxverse(:,:,1,4).', block, mxy_reverse(:,:,1,5).', block, mxy_reverse(:,:,1,6).', block, mxy_reverse(:,:,1,7).')));

FontSize = 18;
cmap = cat(1, jet(256), [1 1 1]);                
figure('Color', 'w', 'Position', [-1 2 1101 814]);
imagesc(abs(im_montage)); axis image off; colormap(gca, cmap); 
hc = colorbar;

set(hc,'FontSize', FontSize);
ylabel(hc, 'Scaled Magnetization [%]'); % ^\circ

text(N2/2+(N2+1)*3  , -10, 'B_0 = 0.55T', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', FontSize+2, 'FontWeight','bold');
% text(N2/2         , 0, sprintf('re-VERSE'), 'Color', 'yellow', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', 'FontSize', FontSize);
text(N2/2         , 0, {sprintf('z = %dcm', 0)},  'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', FontSize);
text(N2/2+(N2+1)  , 0, {sprintf('z = %dcm', 5)}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', FontSize);
text(N2/2+(N2+1)*2, 0, {sprintf('z = %dcm', 10)}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', FontSize);
text(N2/2+(N2+1)*3, 0, {sprintf('z = %dcm', 15)}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', FontSize);
text(N2/2+(N2+1)*4, 0, {sprintf('z = %dcm', 20)}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', FontSize);
text(N2/2+(N2+1)*5, 0, {sprintf('z = %dcm', 25)}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', FontSize);
text(N2/2+(N2+1)*6, 0, {sprintf('z = %dcm', 30)}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', FontSize);

% text(N2/2, (N1+1)*2, sprintf('Proposed'), 'Color', 'yellow', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', 'FontSize', FontSize);
% text(-3, N2/2          , 'Imaginary', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', FontSize, 'rotation', 90);
% text(-3, N2/2+N2+1    ,  'Real'     , 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', FontSize, 'rotation', 90);
% text(-3, N2/2+(N2+1)*2,  'Imaginary', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', FontSize, 'rotation', 90);
% text(-3, N2/2+(N2+1)*3,  'Real'     , 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', FontSize, 'rotation', 90);

%% NRMSE 
% re-verse method 
for i = 1:7 % field strengths
        mxy_ori = mxy_reverse(:,:,i);
        mxy_pro = mxy_maxverse(:,:,i);
        NRMSE_ori(i) = sqrt(sum(sum((abs(mxy_ori) - abs(P_)).^2)))/ sqrt(sum(sum(abs(P_).^2)));
        NRMSE_pro(i) = sqrt(sum(sum((abs(mxy_pro) - abs(P_)).^2)))/ sqrt(sum(sum(abs(P_).^2)));
end

figure;  plot(NRMSE_ori, 'LineWidth', 3, 'Color', [0.4940 0.1840 0.5560]);
hold on; plot(NRMSE_pro, 'LineWidth', 3, 'Color', [0.8500 0.3250 0.0980]);
set(gca, 'FontSize', 20); grid on;
xticklabels({'0cm', '5cm', '10cm', '15cm', '20cm', '25cm', '30cm'});
xlim([0.8 7]);
ylim([0 0.62]);
box off;
ax = gca;
ax.LineWidth = 2;
title('NRMSE'); 
legend('Original', 'Proposed');
