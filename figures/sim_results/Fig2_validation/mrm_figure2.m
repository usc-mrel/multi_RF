% Figure 2: validate the approximation of BS shift 
% display: 
%        ideal |  0.05T  |  0.5T  |   3T   |
%   mag          ori/err   ori/err  ori/err
% phase

%%%%%%%% Ziwei: 11242020
files = dir('results_*');

for i = 1 : length(files)
    load(files(i).name);
    mxy_concomitant_rf  = complex(mxyz_concomitant_rf(:,:,1), mxyz_concomitant_rf(:,:,2));
    mxy_bloch_siegert   = complex(mxyz_bloch_siegert(:,:,1), mxyz_bloch_siegert(:,:,2));
    
    mxy_con(:,:,i) = mxy_concomitant_rf;
    mxy_bs(:,:,i)  = mxy_bloch_siegert;
    
end

mxy_rf = complex(mxyz_rf(:,:,1), mxyz_rf(:,:,2));
% errors amp & phase
err_mag   = abs(mxy_con) - abs(mxy_bs);
err_phase = angle(mxy_con) - angle(mxy_bs);

FOVx = 5e-2;      
dx   = 0.2e-2;    
Nx   = FOVx / dx;
Ny   = Nx; 

%% Showing direct calculation of con. field vs. BS approximation
block = 1.5 * complex(ones(Nx,1, 'double'), ones(Nx,1, 'double'));
mxy_montage1 = cat(2, mxy_rf, block, mxy_con(:,:,1), block,mxy_con(:,:,2), block, mxy_con(:,:,3), block, mxy_con(:,:,4), block, mxy_con(:,:,5));

block = NaN * complex(ones(Nx,1, 'double'), ones(Nx,1, 'double'));
mxy_montage_angle1 = cat(2, mxy_rf, block, mxy_con(:,:,1), block, mxy_con(:,:,2), block, mxy_con(:,:,3), block, mxy_con(:,:,4), block, mxy_con(:,:,5));

block2 =  5 * ones(Nx,Nx, 'double');
block  =  5 * ones(Nx, 1, 'double');
mxy_montage2 = cat(2, block2, block, mxy_bs(:,:,1), block, mxy_bs(:,:,2), block, mxy_bs(:,:,3), block, mxy_bs(:,:,4), block, mxy_bs(:,:,5));

block2 = NaN * ones(Nx,Nx, 'double');
block  = NaN * ones(Nx, 1, 'double');
mxy_montage_angle2 = cat(2, block2, block, mxy_bs(:,:,1), block, mxy_bs(:,:,2), block, mxy_bs(:,:,3), block, mxy_bs(:,:,4), block, mxy_bs(:,:,5));

cmap1 = cat(1, jet(256), [1 1 1]);
cmap2 = cat(1, [1 1 1], hsv(256));
FontSize = 18;

%% Figure 2A: excitation profile comparison
figure('Color', 'w', 'Position', [0 2 1009 814]);
h1 = subplot(4,1,1); 
imagesc(abs(mxy_montage1)); axis image off; colormap(gca, cmap1); 
hc1 = colorbar;
set(hc1, 'FontSize', FontSize);
caxis([min(abs(mxy_con(:))) 1.05]);

text(Ny/2         , 0, '  ', 'Color', 'k', 'FontSize', FontSize, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
text(Ny/2+(Ny+1)  , 0, '0.2T' ,      'Color', 'k', 'FontSize', FontSize, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
text(Ny/2+(Ny+1)*2, 0, '0.55T'  ,    'Color', 'k', 'FontSize', FontSize, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
text(Ny/2+(Ny+1)*3, 0, '1.5T'   ,    'Color', 'k', 'FontSize', FontSize, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
text(Ny/2+(Ny+1)*4, 0, '3T'    ,     'Color', 'k', 'FontSize', FontSize, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
text(Ny/2+(Ny+1)*5, 0, '7T'    ,     'Color', 'k', 'FontSize', FontSize, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
text(0, Nx/2, 'Magnitude', 'Color', 'k', 'FontSize', FontSize, 'Rotation', 90,  'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');

h2 = subplot(4,1,2);
imagesc(abs(mxy_montage2)); axis image off; colormap(gca, cmap1); 
caxis([min(abs(mxy_con(:))) 1.05]);
hc2 = colorbar;
set(hc2, 'FontSize', FontSize);
%text(h2, 0, Nx/2, 'Errors', 'Color', 'k', 'FontSize', FontSize, 'Rotation', 90, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');

h3 = subplot(4,1,3); 
imagesc(angle(mxy_montage_angle1)*180/pi); axis image off; colormap(gca, cmap2);
caxis([-180 180]);
hc3 = colorbar;
set(hc3, 'FontSize', FontSize);
text(h3, 0, Nx/2, 'Phase', 'Color', 'k', 'FontSize', FontSize, 'Rotation', 90, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');

h4 = subplot(4,1,4);
imagesc(angle(mxy_montage_angle2)*180/pi); axis image off; colormap(gca, cmap2); 
caxis([-180 180]);
hc4 = colorbar;
set(hc4, 'FontSize', FontSize);

%% Figure 2B: NRMSE plot using denser samples
cd('./densesamples');
files = dir('results_*');

for i = 1:length(files)

    load(files(i).name);
    mxy_concomitant_rf  = complex(mxyz_concomitant_rf(:,:,1), mxyz_concomitant_rf(:,:,2));
    mxy_bloch_siegert   = complex(mxyz_bloch_siegert(:,:,1), mxyz_bloch_siegert(:,:,2));
    
    mxy_con(:,:,i) = mxy_concomitant_rf;
    mxy_bs(:,:,i)  = mxy_bloch_siegert;

    NRMSE(i) = sqrt(sum(sum((abs(mxy_con(:,:,i)) - abs(mxy_bs(:,:,i))).^2)))/ sqrt(sum(sum(abs(mxy_bs(:,:,i)).^2)));
end
t = [0.005 0.01 0.03 0.05 0.09 0.1 0.2 0.3 0.5 0.7 0.9 1.0 1.5 3 7];


NRMSE = [NRMSE(1:2) NRMSE(4) NRMSE(6:7) NRMSE(9) NRMSE(12:15)];
t     = [t(1:2) t(4) t(6:7) t(9) t(12:15)];

% log scale of NRMSE 
figure; semilogx(t,log2(NRMSE*100),'--','LineWidth',2, 'color', 'k');
hold on; 
semilogx(t,log2(NRMSE*100),'*','LineWidth',2,'MarkerSize',16,'color','blue');
xlim([0.004 8]); ylim([-4 9]);
xlabel('field strengths [T]');
set(gca,'FontSize',20);
grid on;
yticks([log10(1) log10(2) log10(5) log10(10) log10(20) log10(80) log10(250)]);
yticks([log2(1) log2(2) log2(5) log2(10) log2(20) log2(80) log2(250)]);
yticklabels({'1%','2%','5%', '10%', '20%', '80%', '250%'});
xticks((t));
ax = gca;
ax.LineWidth = 1.5;

%% linear scale of NRMSE
% figure; semilogx(t, NRMSE*100, '--','LineWidth',2, 'color', 'k');
% hold on; 
% semilogx(t, NRMSE*100, '*','LineWidth',2,'MarkerSize',16,'color','blue');
% xlim([0.004 8]);
% xlabel('field strengths [T]');
% set(gca,'FontSize',20);
% grid on;
% yticks([(1) (2) (5) (10) (20) (80) (250)]);
% yticklabels({'1%','2%','5%', '10%', '20%', '80%', '250%'});
% xticks((t));
% ax = gca;
% ax.LineWidth = 1.5;
