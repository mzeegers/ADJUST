%%% Spectral Tomography
%%% Thorax phantom of size 512 x 512
%%% 5 materials (3 are identified)
%%% 100 channels 
%%% Full-view setting: 180 angles from 0 to 180
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
phantom = 'Thorax';       % phantom name

%% load phantom

% load matrices for inversion 
[A, F, T, Q] = loadSpectralPhantom(phantom, n);

% number of materials (to be reconstructed)
k = 3;

% problem structure
strProb.A = A;
strProb.F = F;
strProb.k = k;
strProb.T = T;
strProb.n = [n n];

%% generate spectral measurements (on 2x resolution grid)

% create matrices for measurements (phantom discretized on 2x resolution grid, 
% by upsampling the original phantom)
Am = zeros(size(A,1) * 4, size(A,2));
for i = 1:size(A,2) 
    Am(:,i) = reshape(imresize(reshape(A(:,i), n, n), 2, 'nearest'),...
                    n*n*4, 1);
end

% spectral volume
Um     = Am * F;

% tomography operator
n_ang  = 180;
theta  = linspace(0, pi, n_ang); 

volGeo = astra_create_vol_geom(2*n, 2*n);
projGeo= astra_create_proj_geom('parallel', 2, n, theta);
Wm     = opTomo('cuda', projGeo, volGeo);

% spectral tomographic measurements with noise
Y0     = (Wm * Um)/2;         % the factor of 2 division required for 2x grid
Y      = 0*Y0;
for i = 1:size(Y0,2)
   Y(:,i) = astra_add_noise_to_sino(Y0(:,i), Q(i));
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
Jopt.iterMax  = 2000;
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

A1 = [A(:,1), A(:,2), A(:,3) + A(:,4) + A(:,5)];
strProb.A = A1;
strOut = computeResults(strProb, strRU, strUR, strJ, strD);

%% save results

save([saveDir 'example_Thorax.mat'], 'strOut');