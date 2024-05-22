function Afull = STA_maxwell_system_matrix_con(x,k,G,B0,t,Tx,b0,mask,varargin)

% function Afull = STA_maxwell_system_matrix(x,k,t,Tx,b0,mask,varargin)
% Implementation of Grissom's Spatial Domain STA method (doi:10.1002/mrm.20978)
%
% Note: this code is old and not optimized. It is provided only for
% demonstration purposes but should not be viewed as 'state of the art'
%
% Function returns a system matrix Afull
%
%   Inputs:     x - list of voxel coordinates [x(:) y(:) z(:)] can be 2D,
%               3D, or 4D (4th dimension is angular frequency)
%               k - k-space trajectory with same dimensions as x (i.e. if
%               x,y,z specified then k needs kx,ky,kz
%               G - Gradient waveforms (Gx, Gy, Gz) [mT/m]
%               t - time (used for B0 effects) [sec]
%               Tx - transmit sensitivity (no units)
%               b0 - B0 map, uT
%               mask - mask for Tx and b0
%               extra arguments allow setting quadrature (i.e. single
%               channel TX) and to loop the matrix concatenation (slower
%               but less memory reqd)
%
%
% Shaihan Malik 2009 - 2015
% Modified by Namgyun Lee 05/18/2020
% Modified by Ziwei Zhao  12/02/2020
%
% 3-2-09: SJM. Function to define system matrix for STA
% 3-12-09:  Add ability to correct for staircasing of grads
% 30-1-13:  Change default frequency convention to -iwt.

con      = consts;
gamma_uT = con.gamma_uT; % units rad s^-1 uT^-1
gamma_mT = gamma_uT * 1e3;
quad     = false;
two_d    = false;
loopcalc = false; % 9-2-11: allow S*A operation to be done in loop
quad3d   = false;

for ii = 1:length(varargin)
    % 7-11-09, quad mode => Nc=1
    if strcmpi(varargin{ii},'quad')
        quad = true;
    end
    % 28-1-11: force 2d
    if strcmpi(varargin{ii},'2d')
        two_d = true;
    end
    % 9-2-11: allow S*A operation to be done in loop
    if strcmpi(varargin{ii},'loop')||strcmpi(varargin{ii},'loopcalc')
        loopcalc = true;
    end
    % 10-1-13: 3D quad mode
    if strcmpi(varargin{ii},'quad3d')
        quad3d = true;
    end
end
No=1;
switch ndims(Tx)
    case 2
        [Nx Ny] = size(Tx);Nc=1;
        Nz=1;No=1;
    case 3
        if quad
            [Nx Ny No] = size(Tx);
            Nc=1;
        else
            if quad3d
                [Nx Ny Nz] = size(Tx);
                Nc=1;
            else
                [Nx Ny Nc] = size(Tx);
                Nz=1;No=1;
            end
        end
    case 4
        if ~two_d %edit 28-1-11
        [Nx Ny Nz Nc] = size(Tx);
        if (Nc>10)||quad
            No=Nc;Nc=1;
        end

        % 24-11-10: edit this, z dim removed for 3d shim
        if (Ny/Nx)>1.5
            No=Nz;Nz=1;
        end
        else
            [Nx Ny No Nc] = size(Tx);
            Nz=1;
        end
    case 5
        [Nx Ny Nz No Nc] = size(Tx);
end

if isstruct(x)
    xlist=x.list;
else
    xlist=x; % can specify list diectly
end
if isstruct(k)
    klist=k.list;
else
    klist=k;
end

%%% 30-1-13: Add minus sign for frequency axis
if (No>1)
    % assume frequency is last dimension in xlist
    xlist(:,end) = -xlist(:,end);
end

%% set up mask
if ~exist('mask','var')
    mask = ones([Nx Ny Nz No]);
end
idx    = find(mask(:));
Nidx   = length(idx);
% generate indices for multicoil variables
N      = Nx*Ny*Nz*No;
idx_mc = logical(repmat(mask,[1 1 1 1 Nc]));

%% look in args...
% assume dt = t(2)-t(1)
if numel(t)~=1
    dt = t(2)-t(1);
    T = t(end);
else
    dt=1;T=0; % arbitrary
end

%% set up concomitant field terms 
xx    = xlist(idx,1);      % [m]
yy    = xlist(idx,2);      % [m]
zz    = xlist(idx,3);      % [m]
gx    = G(:,1).';       % [mT/m]
gy    = G(:,2).';       % [mT/m]
gz    = G(:,3).';       % [mT/m]
B0    = B0*1e3;         % [mT]!!

% % second order: G*xlist [Nx * Nt]
% Con_term1 = 1/(2*B0) * (1/4*xx.^2*gz.^2 + 1/4*yy.^2*gz.^2 + ...
%                 zz.^2*(gx.^2+gy.^2)) + ...
%             1/(B0)   * (-1/2*(yy.*zz)*(gy.*gz) - 1/2*(xx.*zz)*(gx.*gz));
% 
% Con_term2 = -1/(2*B0.^2) * (1/4*xx.^3*(gz.^2.*gx) + 1/4*yy.^3*(gz.^2.*gy) + ...
%                 zz.^3*(gz.*(gx.^2+gy.^2))) + ...
%             -1/(2*B0.^2) * (1/4*(xx.^2.*yy)*(gz.^2.*gy) + ...
%                 (xx.^2.*zz)*(1/4*gz.^3-gx.^2.*gz));
% 
% Con_term3 = -1/(2*B0.^2) * (yy.^2.*xx) * (1/4*gz.^2.*gx);
% Con_term4 = -1/(2*B0.^2) * (yy.^2.*zz) * (1/4*gz.^3-gy.^2.*gz);
% 
% Con_term5 = -1/(2*B0.^2) * (zz.^2.*xx) * (gx.*(gx.^2+gy.^2)-gz.^2.*gx);
% Con_term6 = -1/(2*B0.^2) * (zz.^2.*yy) * (gy.*(gx.^2+gy.^2)-gz.^2.*gy);
% Con_term7 =  1/(B0.^2)   * (xx.*yy.*zz) * (gx.*gy.*gz);
% 
% Con_term = Con_term1 + Con_term2 + Con_term3 + Con_term4 + ...
%            Con_term5 + Con_term6 + Con_term7; % mT
       
%% spatial basis
h = zeros(length(xx),15);
h(:,1)  = xx.^2;
h(:,2)  = yy.^2;
h(:,3)  = zz.^2;
h(:,4)  = yy .* zz;
h(:,5)  = xx .* zz;
h(:,6)  = xx.^3;
h(:,7)  = yy.^3;
h(:,8)  = zz.^3;
h(:,9)  = xx.^2 .* yy;
h(:,10) = xx.^2 .* zz;
h(:,11) = xx.* yy.^2;
h(:,12) = yy.^2 .* zz;
h(:,13) = xx.* zz.^2;
h(:,14) = yy.* zz.^2;
h(:,15) = xx .* yy .* zz;
      
%%  higher order gradient terms
c = zeros(length(gx), 15);
 
c(:,1)  = ( 1 / (8 * B0) * gz.^2);                                     % x2
c(:,2)  = ( 1 / (8 * B0) * gz.^2);                                     % y2
c(:,3)  = ( 1 / (2 * B0) * (gx.^2 + gy.^2));                           % z2
c(:,4)  = (-1 / (2 * B0) * gy .* gz);                                  % yz
c(:,5)  = (-1 / (2 * B0) * gx .* gz);                                  % xz
c(:,6)  = (-1 / (8 * B0^2) * gx .* gz.^2);                             % x3
c(:,7)  = (-1 / (8 * B0^2) * gy .* gz.^2);                             % y3
c(:,8)  = (-1 / (2 * B0^2) * gz .* (gx.^2 + gy.^2));                   % z3
c(:,9)  = (-1 / (8 * B0^2) * gy .* gz.^2);                            % x2y
c(:,10) = (-1 / (2 * B0^2) * (1 / 4 * gz.^3 - gx.^2 .* gz));          % x2z
c(:,11) = (-1 / (8 * B0^2) * gx .* gz.^2);                            % xy2
c(:,12) = (-1 / (2 * B0^2) * (1 / 4 * gz.^3 - gy.^2 .* gz));          % y2z
c(:,13) = (-1 / (2 * B0^2) * (gx .* (gx.^2 + gy.^2) - gx .* gz.^2));  % xz2
c(:,14) = (-1 / (2 * B0^2) * (gy .* (gx.^2 + gy.^2) - gy .* gz.^2));  % yz2
c(:,15) = (1  / (B0^2)     * gx .* gy .* gz);                         % xyz

% kk = cumsum(gamma_mT * c * dt); % check gamma [rad/s/mT]; dt [s]
for ic = 1:15
    gg = c(:,ic);
    kk(:,ic) = -gamma_mT*flipud(cumtrapz(t, flipud(gg)));
end

phi = kk * h.';
E   = exp(1i * phi);  % should be minus?

%% set up b0 - assume in uT
dB0mat = 1e-3 * squeeze(b0(idx))*(t-T); % convert from uT to mT
M0     = 1;
% by default use minus sign for frequency
dB0mat = 1i * gamma_mT * M0 * dt * exp(-1i * gamma_mT * dB0mat);

%% ---- Fourier Matrix -----
F = xlist(idx,:) * klist';        % this is a matrix of dot products
F = exp(1i*F);
A = dB0mat .* F .* E.';           % element wise multiplication

%% do block diagonal repeat for coils
%%% 28-10-13: no need to do block repeat if using loopcalc
if ~loopcalc
a = cell(1,Nc);
if ~strcmp(class(A),'single')
    [a{:}] = deal(sparse(A)); % local machine
else
    [a{:}] = deal(A);         % this is an efficient block diag trick
end
%if ~loopcalc
    Ab = (blkdiag(a{:}));
%end
end
%% ---- sensitivity information -------
% use sparse matrices
row_idx = repmat((1:Nidx).',[Nc 1]);
col_idx = (1:(Nidx*Nc)).';
% 28-10-13: only do this if not loopcalc
if ~loopcalc
    S = sparse(row_idx,col_idx,squeeze(double(Tx(idx_mc))),Nidx,Nidx*Nc,Nidx*Nc);
end

%% system matrix A
if ~loopcalc
    Afull = S*Ab;
else
    M = size(A,2);
    Afull=zeros([Nidx Nc*M]);
    for ii=1:Nc
        idx2 = (1:M) + (ii-1)*M;
        tmp=squeeze(Tx(idx+(ii-1)*N));
        S = repmat(tmp(:),[1 M]);
        Afull(:,idx2) = S.*A;
    end
end
%%
end
