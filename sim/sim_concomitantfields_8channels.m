%% Simulation of concomitant fields using Bloch Siegert shift
%%% proposed correction method for 2D spiral with measured GIRFs
% Ziwei Zhao 

function sim_concomitantfields_8channels(B0, zoff, T2, fov, concomitant_correct, display_figure, paths)
    
    %% Constant definitions
    gamma_uT = 267.5221;       % [rad/sec/uT]

    flip = 90;                 % total flip angle [degrees]
    flip = flip * pi / 180;

    %% Design initial spiral trajectory
    T          = 26e-3;         % pulse duration [sec]
    dt         = 6.4e-6;        % RF/gradient raster time [sec]
    t          = [0:dt:T-dt]';  % seconds, time vector
    dr         = 0.25;          % resolution of trajectory [cm]
    kmax       = 1/2/dr;        % [cycles/cm], max radius in k-space

    % spiral waveforms 
    N = kmax/(1/fov);           % number of turns, no acceleration
    k = 2*pi*kmax*(1-t/T).*exp(1i*2*pi*N*(1-t/T)); % [rad/cm] - kx = real(k); ky = imag(k);
    k = [real(k) imag(k) zeros(length(k),1)];
    k = k*100;                  % [rad/cm] -> [rad/m]

    K = k;

    %% Load in B0 & B1 field maps + FOV info (example data from 3T 8ch head coil)
    load(fullfile(paths.root_path, '/multi_RF/third_party/phase_relaxed_CPMG_excitation/b0_b1_maps_2d.mat'));

    Z = Z * 0;
    [N1, N2, Nc] = size(tx);
    xlist = cat(2, X(:), Y(:), Z(:) + zoff); % N1 x N2 x 3
    idx = find(m); % index of non-masked voxels

    %% Define a square target
    xoff = +6e-3;  % [mm] -> [m]
    yoff = 0e-3;   % [mm] -> [m]

    r0   = 15e-3;  % [mm] -> [m]
    r0   = 1.74 * r0;

    % Square beam
    P  = double(((abs(X - xoff) < r0) & (abs(Y - yoff) < r0)));

    % apodize
    h  = fspecial('gaussian', 3, 0.5);
    P_ = imfilter(P, h);

    % Now scale to flip
    P  = P_ * flip * 1j; % desired excitation pattern

    %% Define gradient correction model (GIRF)
    girf = load(fullfile(paths.root_path, '/multi_RF/third_party/reVERSE-GIRF/GIRF_3T_London_20140729.mat'));

    % This correction gives G and k after convolution with GIRF
    gcor = @(x)(gradient_distort_GIRF(x, girf.ff, girf.Hw, dt, 10));

    %% Example design: Include GIRF
    if concomitant_correct 
        script_name = 'blochmex'; 
    else   
        script_name = 'bloch'; 
    end

    opt        = reVERSE_init;
    dt         = opt.dt;        % sampling dwell time [usec]
    opt.lambda = 1;
    opt.Nstop  = 20;
    opt.show   = 1;

    alpha = 0.5;
    g     = 0;

    T1 = 1e6; % T1 relaxation time [sec]

    mxyz_offcenter = zeros(N1, N2, 3, 'double');

    % scale the b0 map 
    % b0 = b0 / 3.0 * B0; % please uncomment to generate Figure 4

    % Set up function for system matrix
    if concomitant_correct
        afun  = @(k, G)(STA_maxwell_system_matrix_con(xlist, k, G, B0, opt.dt * (1:size(k,1)), tx, b0, m, 'loopcalc'));
        tic;
        [bb, Gv] = reVERSE_GIRF_con(P(idx), K, afun, gcor, opt);
        Time_reVERSE_con = toc;
    else
        afun2 = @(k)(STA_system_matrix(xlist, k, opt.dt * (1:size(k,1)), tx, b0, m, 'loopcalc'));
        tic;
        [bb,Gv] = reVERSE_GIRF(P(idx), K, afun2, gcor, opt);
        Time_reVERSE = toc;
    end

    % Outputs
    rf = bb{end};   % [mT]
    G  = Gv{end};   % [mT/m]

    duration = dt*1e3*length(G);

    %% Perform Bloch simulation
    N_lambda = length(idx); % number of voxels within a mask

    mx0 = zeros([N1 N2], 'double');
    my0 = zeros([N1 N2], 'double');
    mz0 = zeros([N1 N2], 'double');
    mz0(idx) = 1;

    [Gedd, k] = gcor(G);  % [mT/m]

    %%
    % perform bloch equations
    start_time = tic;
    for lambda = 1 : N_lambda

        [idx1,idx2] = ind2sub([N1 N2], idx(lambda));
        fprintf('Performing Bloch simulation at (%3d/%3d)... ', lambda, N_lambda);

        % rf_out: [mT]   * [T/1e3mT] * [1e4G/T]             => *1e1 [G]
        % G_out : [mT/m] * [T/1e3mT] * [1e4G/T] * [m/1e2cm] => *1e-1[G/cm]

        rf_combined = sum(bsxfun(@times, rf, reshape(tx(idx1,idx2,:), [1 Nc])), 2); 
        df = b0(idx1, idx2) * gamma_uT / (2 * pi); % [uT] * [rad/sec/uT] * [cycle/2pi rad] => [Hz]
                
        % off iso-center
        dp = cat(3, X(idx1,idx2), Y(idx1,idx2), Z(idx1,idx2) + zoff) * 1e2; % [m] * [1e2cm/m] => [cm]
        
        [mx,my,mz] = bloch_maxwell(rf_combined*1e1, Gedd*1e-1, dt, T1, T2, df, dp, 0, B0, alpha, g, mx0(idx1,idx2), my0(idx1,idx2), mz0(idx1,idx2));
        
        mxyz_offcenter(idx1, idx2, :) = cat(3, mx, my, mz);
        fprintf('done! (%5.4f sec)\n', toc(start_time));
    end

    % save the results
    cur_dir = pwd;  
    cd(paths.output_path); % this is set up from the script
    save(sprintf('%s_B0%04.2f_offc%04.1fcm_iter%d_dur%.3f_dr%.2fcm_lambda%1.0f_T2%02.0fms.mat', script_name, B0, zoff*1e2, opt.Nstop, duration, dr, opt.lambda, T2*1e3), 'mxyz_offcenter');
    cd(cur_dir);

    %%
    % Check NRMSE
    mxy  = squeeze(complex(mxyz_offcenter(:,:,1,:), mxyz_offcenter(:,:,2,:)));
    mxy_ori = mxy;
    NRMSE = sqrt(sum(sum((abs(mxy_ori) - abs(P_)).^2)))/ sqrt(sum(sum(abs(P_).^2)));
    fprintf('NRMSE: %f\n', NRMSE);

    if display_figure
        ind_t = 0:dt:(length(G)-1)*dt;
        ind_t = ind_t * 1e3;
        
        figure;  plot(ind_t, abs(rf(:,1)), 'LineWidth', 2);
        hold on; plot(ind_t, abs(rf(:,2)), 'LineWidth', 2);
        hold on; plot(ind_t, abs(rf(:,3)), 'LineWidth', 2);
        hold on; plot(ind_t, abs(rf(:,4)), 'LineWidth', 2);
        hold on; plot(ind_t, abs(rf(:,5)), 'LineWidth', 2);
        hold on; plot(ind_t, abs(rf(:,6)), 'LineWidth', 2);
        hold on; plot(ind_t, abs(rf(:,7)), 'LineWidth', 2);
        hold on; plot(ind_t, abs(rf(:,8)), 'LineWidth', 2);
        
        legend('Channel - 1', 'Channel - 2', 'Channel - 3', 'Channel - 4', 'Channel - 5',...
            'Channel - 6', 'Channel - 7', 'Channel - 8');
        ylabel('Amplitude [mT]'); xlabel('time [ms]');
        set(gca, 'FontSize', 16); title('Multi-channel RF @0.55T');
        
        figure; subplot(1,2,1);
        plot(ind_t, Gedd(:,1), 'LineWidth', 2);
        hold on; plot(ind_t, Gedd(:,2), 'LineWidth', 2);
        legend('Gx', 'Gy');
        xlabel('time [ms]'); ylabel('[mT/m]');
        xlim([0 18]); box off; grid on;
        set(gca, 'FontSize', 18);
        title('Gradient waveforms');
        subplot(1,2,2);
        plot(k(:,1), k(:,2), 'LineWidth', 2);
        box off; grid on;
        xlabel('[rads/m]'); ylabel('[rads/m]');
        xlim([-1500 1500]); ylim([-1500 1500]);
        set(gca, 'FontSize', 18); 
        title('Excitation K-space');
    end

end
