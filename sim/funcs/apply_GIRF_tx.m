function [GPred, kPred_tx] = apply_GIRF_tx(gradients_nominal, dt, R, tRR)
% function [kPred, GPred] = apply_GIRF_tx(gradients_nominal, dt, R, tRR)
%
% % UNITS
% gradients_nominal [ G/cm ]
% dt                [ s ]
%
% Hack to R to handle field strength (backwards compatibility)
% R.R = rotation matrix;
% R.T = field strength {}
%
% tRR is sub-dwell-time offset (-1 to 1) [0]?


% handle "nasty" co-opting of R-variable to include field info.
if isstruct(R)
    field_T = R.T;
    R = R.R;
else
    field_T = 0.55;
end

%% LOAD GIRF (Field/scanner dependent)
% Field selection needs to be handled upstream, or include the headers
% in to this file

if field_T == 1.4940
    % 1.5T Aera (NHLBI 2016)
    girf_file = 'GIRF_20160501.mat';
elseif field_T == 0.55
    % 0.55T Aera (NHLBI 2018)
    girf_file = '.../multi_RF/b0b1_map_055T/GIRF_20200221_Duyn_method_coil2.mat';
    % girf_file = ones(3800,3);
end

% === Load file ===
try
    load(girf_file);
    disp(['Using ' girf_file]);
    %GIRF(:,2) = GIRF(:,1); % Hack by NGL
catch
    disp('Couldn''t find the GIRF file, using ones..');
    GIRF = ones(3800,3);
end

%%
dtGIRF = 10e-6;
dtSim  = dt;
l_GIRF = length(GIRF);
[samples, interleaves, ~] = size(gradients_nominal);

% if readout is real long, need to pad the GIRF measurement
if samples*dt > dtGIRF*l_GIRF
    disp('readout length > 38ms, zero-padding GIRF');
    pad_factor = 1.5 * (samples * dt) / (dtGIRF * l_GIRF); % 1.5 factor to ease calculations below
    new_GIRF = zeros(round(l_GIRF * pad_factor), 3);

    for i = 1:3
        fft_GIRF = fftshift(ifft(ifftshift(GIRF(:,i))));
        zeropad = round( abs( (l_GIRF-length(new_GIRF) ) /2 ));
        temp = zeros(length(new_GIRF),1);
        % smoothing of padding:
        H_size = 200;
        H = hanningt(H_size);
        fft_GIRF(1:(H_size/2)) = fft_GIRF(1:(H_size/2)).*reshape(H(1:(H_size/2)),size(fft_GIRF(1:(H_size/2))));
        fft_GIRF(end-(H_size/2 - 1):end) = fft_GIRF(end-(H_size/2 - 1):end).*reshape(H((H_size/2 + 1):H_size),size(fft_GIRF(end-(H_size/2 - 1):end)) );

        temp((1+zeropad):(zeropad + l_GIRF))= fft_GIRF;
        % figure, plot(real(temp))
        new_GIRF(:,i) = fftshift(fft(fftshift(temp)));
    end

    old_GIRF = GIRF;
    GIRF = new_GIRF;
    l_GIRF = length(GIRF);

end

%%   GIRF prediction  
%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 4
    tRR = 0;
end

ADCshift = 0.85e-6 + 0.5 * dt + tRR * dt; % NCO Clock shift % 

%%
clear G0 GNom GPred kNom kPred

% allocation
Nominal   = zeros(samples, 3, 'double');
Predicted = zeros(samples, 3, 'double');
GNom      = zeros(samples, 3, interleaves, 'double');
GPred     = zeros(samples, 3, interleaves, 'double');
kNom      = zeros(samples, 3, interleaves, 'double');
kPred     = zeros(samples, 3, interleaves, 'double');

% GIRF process
for l = 1:interleaves
    % Select current spiral arm
    G0(:,1) = gradients_nominal(:,l,1);
    G0(:,2) = gradients_nominal(:,l,2);
    G0(:,3) = zeros(size(gradients_nominal,1),1);
    % Rotate into physical coordinates
    G0 = (R * G0.').';

    %--Loop through x,y,z gradient trajectories--%
    for ax = 1:3

        % Zeropad in time domain to match frequency resolution of GIRF (match readout length)
        L = round(dtGIRF * l_GIRF / dtSim); % when waveform not at GRT
        G = zeros(L,1, 'double');

        % REQUIRES SOME INTELLIGENCE TO DETERMINE WHICH SPIRAL FORM

        N = length(G0);
        index_range = (-floor(N/2):ceil(N/2)-1).' + floor(L/2) + 1;
        G(index_range) = G0(:,ax);

        % Make a waveform periodic by returning to zero
        H = G(index_range(N)) * hanningt(400);
        G(index_range(N)+(1:length(H)*0.5)) = H(length(H)*0.5+1:end);
        %FFT nominal gradient
        dw = 1 / (L * dtSim); % frequency resolution [Hz]
        w = (-floor(L/2):ceil(L/2)-1).' * dw; % [Hz]
        I = fftshift(fft(ifftshift(G)));

        %Zeropad GIRF and I to bandwidth of sampling (when waveform not at GRT)
        GIRF1 = zeros(L,1, 'double');

        if dt > dtGIRF
            % RR crop
            GIRF1 = GIRF(round(l_GIRF/2 - L/2 + 1):round(l_GIRF/2 + L/2),ax);
            % RR .. padding
            temp = hanningt(10);
            GIRF1(1) = 0; GIRF1(end) = 0;
            GIRF1(2:round(length(temp)/2) + 1) = GIRF1(2:round(length(temp)/2)+1).*reshape(temp(1:round(length(temp)/2)),size(GIRF1(2:round(length(temp)/2)+1)));
            GIRF1(end-round(length(temp)/2):end-1) = GIRF1(end-round(length(temp)/2):end-1).*reshape(temp((round(length(temp)/2) + 1):length(temp)),size(GIRF1(end-round(length(temp)/2):end-1)));
        else
            index_range = (-floor(l_GIRF/2):ceil(l_GIRF/2)-1).' + floor(L/2) + 1;
            GIRF1(index_range) = GIRF(:,ax);
        end

        % Predict Gradient and apply clock shift
        P = I .* GIRF1 .* exp(1j * ADCshift * 2 * pi * w);

        % zeropad to required bandwidth (when waveform is at GRT)
        BW = 1 / dt;
        L  = round(BW/dw);    % number of points required

        PredGrad = zeros(L,1, 'double');
        NomGrad  = zeros(L,1, 'double');
        zeropad  = round(abs((length(G)-L)/2)); % amount of zeropadding

        PredGrad((1+zeropad):(zeropad+length(P)))= P;
        NomGrad((1+zeropad):(zeropad+length(I)))= I;

        % FFT back to time domain
        PredGrad = fftshift(ifft(ifftshift(PredGrad))) * L / length(G);
        NomGrad  = fftshift(ifft(ifftshift(NomGrad)))  * L / length(G);

        %Correct polarity of gradients
        multiplier = zeros(length(PredGrad),1);
        for i = 1:length(PredGrad)
            if real(PredGrad(i))>0; multiplier(i) = 1;
            else; multiplier(i) = -1;
            end
        end
        PredGrad = abs(PredGrad).*multiplier;

        multiplier = zeros(length(NomGrad),1);
        for i = 1:length(NomGrad)
            if real(NomGrad(i))>0; multiplier(i) = 1;
            else; multiplier(i) = -1;
            end
        end
        NomGrad = abs(NomGrad).*multiplier;

        index_range = (-floor(samples/2):ceil(samples/2)-1).' + floor(L/2) + 1;
        Nominal(:,ax) = NomGrad(index_range);
        Predicted(:,ax) = PredGrad(index_range);
    end

    %rotate back to logical coordinates
    GNom(:,:,l)= (R.'*Nominal.').';
    GPred(:,:,l) = (R.'*Predicted.').';

    %Integrate to get k-space trajectory from gradient
    kNom(:,:,l)  = cumsum(GNom(:,:,l));
    kPred(:,:,l) = cumsum(GPred(:,:,l));
end

% add ktx_pred
gamma_mT = 2*pi*42577.46778; % rad s^-1 mT^-1 [rad/mT]
M = length(GPred);
t = (0:M-1)*dt;
kPred_tx = -gamma_mT*flipud(cumtrapz(t,flipud(GPred))); % rad/mT * mT/m = rad/m


end

function wind = hanningt(windowLength)
% If license toolbox being an idiot..
% ripped from: https://www.mathworks.com/matlabcentral/fileexchange/48925-hann-window

if license('checkout','Signal_Toolbox')
    wind = hanning(windowLength);
else
    N = windowLength - 1;
    num = linspace(0,N,windowLength);
    wind =  0.5*(1 - cos(2*pi*num/N));
end

% matlab:
% function w = sym_hanning(n)
% %   SYM_HANNING Symmetric Hanning window.
% %   SYM_HANNING Returns an exactly symmetric N point window by evaluating
% %   the first half and then flipping the same samples over the other half.
%
% if ~rem(n,2)
%    % Even length window
%    half = n/2;
%    w = calc_hanning(half,n);
%    w = [w; w(end:-1:1)];
% else
%    % Odd length window
%    half = (n+1)/2;
%    w = calc_hanning(half,n);
%    w = [w; w(end-1:-1:1)];
% end
%
% %---------------------------------------------------------------------
% function w = calc_hanning(m,n)
% %   CALC_HANNING Calculates Hanning window samples.
% %   CALC_HANNING Calculates and returns the first M points of an N point
% %   Hanning window.
%
% w = .5*(1 - cos(2*pi*(1:m)'/(n+1)));

end