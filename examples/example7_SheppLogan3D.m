%%% Spectral Tomography
%%% 3D Shepp Logan phantom of size 128 x 128 x 128
%%% 3 materials 
%%% 100 channels 
%%% Full-view setting: parallel beam geometry, 60 angles from 0 to 180 deg
%%% Non-inverse crime configuration: measurements are simulated on 2x
%%% high-resolution grid and corrupted with mild Poisson noise.
%%% The phantom by Matthias Schabel can be found at: 
%%% https://nl.mathworks.com/matlabcentral/fileexchange/9416-3d-shepp-logan-phantom
%%% Please note:
%%% The spectral reconstruction algorithms requires around 32-48 GB memory
%%% and hence, we advice to run this script only if you have PC with 64 GB
%%% RAM. Otherwise, please reduce the size of phantom to either 32 or 64
%%% px.
%
% Authors:
%   Ajinkya Kadu,
%       Centrum Wiskunde & Informatica, Amsterdam (aak@cwi.nl)
%   Mathé Zeegers, 
%       Centrum Wiskunde & Informatica, Amsterdam (M.T.Zeegers@cwi.nl)

clc; clearvars; close all;

% setting random stream
myStream = RandStream('mt19937ar', 'Seed', 10);
RandStream.setGlobalStream(myStream)

% directory for saving results
saveDir = [pwd '\results\'];
if ~exist(saveDir, 'dir'), mkdir(saveDir); end

%% settings

n       = 128;            % spatial resolution of phantom
phantom = 'SheppLogan3D'; % phantom name

%% load phantom (for measurement and inversion)

% load matrices for measurements (phantom discretized on 2x resolution grid)
[Am, Fm, Tm, Qm] = loadSpectralPhantom(phantom, 2*n);

% load matrices for inversion
[A, F, T, Q] = loadSpectralPhantom(phantom, n);

% number of materials (in phantom)
k = size(A, 2);

% problem structure
strProb.A = A;
strProb.F = F;
strProb.k = k;
strProb.T = T;
strProb.n = [n n n];

%% generate spectral measurements (on 2x resolution grid)

% spectral volume
Um     = Am * Fm;

% tomography operator
n_ang  = 60;
theta  = linspace(0, pi, n_ang);  

volGeo = astra_create_vol_geom(2*n, 2*n, 2*n);
projGeo= astra_create_proj_geom('parallel3d', 2, 2, n, n, theta);
Wm     = opTomo('cuda', projGeo, volGeo);

% spectral tomographic measurements with noise
Y0     = (Wm * Um)/2;      % the factor of 2 division required for 2x grid
Y      = 0*Y0;
for i = 1:size(Y0, 2)
   Y(:, i) = astra_add_noise_to_sino(Y0(:, i), Qm(i));
end

%% (non-inverse-crime) forward operator to be used for inversion

volGeo = astra_create_vol_geom(n, n, n);
projGeo= astra_create_proj_geom('parallel3d', 1, 1, n, n, theta);
W      = opTomo('cuda', projGeo, volGeo);

% (only for checking) measurments with forward operator
U  = A * F;
Y1 = W * U;

strProb.angles = theta;
strProb.Y      = Y;
strProb.Yclean = Y1;

clear Am Fm Wm Um projGeo volGeo Y0 Y1

%% Spectral Tomography algorithms

%%% Reconstruction and then Unmixing
[Aru, Fru] = RU(Y, W, k);

strRU.A = Aru;
strRU.F = Fru;

%%% Unmixing then Reconstruction
[Aur, Fur] = UR(Y, W, k);

strUR.A = Aur;
strUR.F = Fur;

%%% cJoint
Jopt.rho      = 1e-2;
Jopt.iterMax  = 500;
[Aj, Fj, histJ] = cJoint(Y, W, k, Jopt);

strJ.A = Aj;
strJ.F = Fj;

%%% ADJUST
Dopt.rho     = 1e-2;
Dopt.iterMax = 500;
[Ad, Fd, Rd, hist] = ADJUST(Y, W, k, T, Dopt);

strD.A = Ad;
strD.F = Fd;

%% plots and tables

strOut = computeResults(strProb, strRU, strUR, strJ, strD);

%% save results

save([saveDir 'example_SheppLogan3D.mat'], 'strOut');