%% 22-4-15: Initialize options
function opt = reVERSE_init_055T

opt = struct;

%%% Lambda value (regularization)
opt.lambda = 1;

%%% Number of channels
opt.Nc = 1;

%%% B1 limit for VERSE
% [uT] * [mT/1e3uT] * [T/1e3mT] * [1e4G/T] => *1e1[G]
opt.b1_limit = 13.3e-3 * 10;  % <--- limit in Gauss
opt.b1_alpha = 0.9;         % <--- At each iteration of VERSE we limit to b1_limit*b1_alpha
opt.bmax     = opt.b1_limit * opt.b1_alpha;  %<--- overall limit is limit*alpha

%%% Max iterations for reVERSE
opt.Nstop = 10;

%%% Sampling dwell, sec
opt.dt = 6.4e-6;
%opt.dt = 8e-9;

%%% Max and slew rate limits
opt.Gmax = 24;          % mT/m
opt.Smax = 180.18;      % T/m/s

%%% Display information
opt.show = false;
end