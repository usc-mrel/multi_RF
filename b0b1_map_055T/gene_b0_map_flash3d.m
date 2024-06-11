% demo B0 maps from FLASH sequence
%
% Ziwei Zhao 01032023
% recon from rawdata

clear all;
close all;

addpath(genpath('.../multi_RF/third_party/mapVBVD'));
addpath(genpath('.../multi_RF/sim/funcs'));

%% Setup BART library 
cur_dir  = pwd;
base_dir = ['/Users/' getenv('USER')];
bart_dir = [base_dir '/Documents/matlab/bart-0.8.00']; % change this to your installation
cd(bart_dir);
startup;
cd(cur_dir);

%% 3D FLASH Recon
root_path = '...';
vol_num   = '03312024_off15cm';
TEs = [37 47 57 67 77] * 1e-3; % [s]
Necho = 5;

rawdata_path = fullfile(root_path, vol_num, '/raw_data');
cd(rawdata_path);

% read in siemens .dat data
rawdata = dir(fullfile(rawdata_path, './*flash_2*.dat')); 

for i = 1 : length(rawdata)
    
    tmp_raw=rawdata(i).name;
    obj=mapVBVD(tmp_raw,'ignoreSeg');
    raw2=obj{1,2}.image{''};
    Dim=obj{1,2}.image.sqzDims;

    [nread, nc, npheco, nslices] = size(raw2);
    recon_imgc = zeros(nread, npheco, nc, nslices);

    %% generate per coil images using fft
    recon_nslices = 1 / sqrt(size(raw2, 4)) * fftshift(fft(ifftshift(raw2, 4), [], 4), 4);

    for inslice = 1:nslices
        for icoil = 1:nc
            recon_imgc(:,:,icoil,inslice) = fft2c(squeeze(recon_nslices(:,icoil,:,inslice)));
        end
    end

    %% remove readout oversampling
    recon_img_ori(:,:,:,:,i) = recon_imgc(size(recon_imgc,1)/4+1:size(recon_imgc,1)/4*3,:,:,:);

end

%% coil maps using ESPIRIT - only first kspace data 
% recalculate kspace
ksp_tmp = zeros(size(recon_img_ori(:,:,:,:,1)));

for inslice = 1:nslices
    for icoil = 1:nc
        recon_img2d(:,:,icoil,inslice) = fft2c(squeeze(recon_img_ori(:,:,icoil,inslice, 1)));
    end
end

ksp_ori = 1/sqrt(size(recon_img2d, 4)) * fftshift(fft(ifftshift(recon_img2d, 4), [], 4), 4);

calib_size = [11 11 0];
% [22 22 18] for 128 by 128 by 104
cmaps = bart(sprintf('ecalib -t 0.02 -c0 -r%d:%d -m1 -I -d2', calib_size(1), calib_size(2)), permute(ksp_ori, [1 2 4 3]));

for i = 1:length(rawdata)
    recon_img(:,:,:,i) = squeeze(sum(conj(cmaps).* permute(recon_img_ori(:,:,:,:,i), [1 2 4 3]), 4)) ./ squeeze(sqrt(sum(abs(cmaps).^2, 4)));
end

for iecho = 1:5
    recon_img_up(:,:,1,iecho) = imresize(recon_img(:,:,1,iecho), 1.5);
end

as(recon_img_up);

%% phase unwrapping
img_phase = angle(recon_img); % [rad]
img_phase = unwrap(img_phase, [], 4); % N1 x N2 x N3 x Ne 
N1 = size(img_phase, 1);
N2 = size(img_phase, 2);
N3 = size(img_phase, 3);

%% normalize the phase
img_phase = img_phase - repmat(img_phase(:,:,:,1), [1 1 1 Necho]);

%% check phase images
img_phase_deg = img_phase * 180 / pi;
as(img_phase_deg)

%% Perform a linear least-squares fit of the phase-time curve
y = reshape(permute(img_phase, [4 1 2 3]), [Necho N1*N2*N3]); % [rad]
A = cat(2, TEs', ones(Necho,1));
x_ls = inv(A.' * A) * (A.' * y); % 2 x N1*N2*N3

fieldmap = reshape(x_ls(1,:) / (2 * pi), [N1 N2 N3]); % [Hz]

as(fieldmap);

figure; imshow(fieldmap, [-5 5]);
set(gca, 'FontSize', 16, ...
    'XColor', 'w', 'YColor', 'w', 'GridColor', 'w', ...
    'Color', 'k'); % has some problems..
axis image;
colormap(hot(256));
colorbar;

% fieldmap_tmp = imrotate(fieldmap, 270); % should match with the B1 map orientation
save('fieldmap_b0_fov220mm_matrix64_off_15cm_03312024.mat', 'fieldmap', 'recon_img');
