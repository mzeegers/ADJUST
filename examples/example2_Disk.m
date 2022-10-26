%%% Spectral Tomography
%%% Disk phantom of size 512 x 512
%%% 8 materials 
%%% 100 channels 
%%% Full-view setting: 180 angles from 0 to 180 degrees
%%% Non-inverse crime configuration: measurements are simulated on 2x
%%% high-resolution grid and corrupted with mild Poisson noise.
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

n       = 512;            % spatial resolution of phantom
phantom = 'Disks';        % phantom name

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
strProb.n = [n n];

%% generate spectral measurements (on 2x resolution grid)

% spectral volume
Um     = Am * Fm;

% tomography operator
n_ang  = 180;
theta  = linspace(0, pi, n_ang); 

volGeo = astra_create_vol_geom(2*n, 2*n);
projGeo= astra_create_proj_geom('parallel', 2, n, theta);
Wm     = opTomo('cuda', projGeo, volGeo);

% spectral tomographic measurements with noise
Y0     = (Wm * Um)/2;      % the factor of 2 division required for 2x grid
Y      = 0*Y0;
for i = 1:size(Y0,2)
   Y(:,i) = astra_add_noise_to_sino(Y0(:,i), Qm(i));
end

%% (non-inverse-crime) forward operator to be used for inversion

volGeo = astra_create_vol_geom(n, n);
projGeo= astra_create_proj_geom('parallel', 1, n, theta);
W      = opTomo('cuda', projGeo, volGeo);

% (only for checking) measurments with forward operator
U  = A * F;
Y1 = W * U;

strProb.angles = theta;
strProb.Y      = Y;
strProb.Yclean = Y1;

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
Jopt.iterMax  = 1000;
Jopt.showIter = 1;
[Aj, Fj, histJ] = cJoint(Y, W, k, Jopt);

strJ.A = Aj;
strJ.F = Fj;

%%% ADJUST
Dopt.rho      = 1e-2;
Dopt.iterMax  = 1000;
Dopt.showIter = 1;
[Ad, Fd, Rd, hist] = ADJUST(Y, W, k, T, Dopt);

strD.A = Ad;
strD.F = Fd;

%% plots and tables

strOut = computeResults(strProb, strRU, strUR, strJ, strD);

%% save results

save([saveDir 'example_Disk.mat'], 'strOut');