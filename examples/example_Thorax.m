%%% Spectral Tomography
%%% Thorax phantom of size 512 x 512
%%% 5 materials (3 soft + 2 hard)
%%% 100 channels 
%%% Full-view setting: 60 angles from 0 to 180 degrees
%
% Authors:
%   Ajinkya Kadu,
%       Centrum Wiskunde & Informatica, Amsterdam (aak@cwi.nl)
%   Mathé Zeegers, 
%       Centrum Wiskunde & Informatica, Amsterdam (M.T.Zeegers@cwi.nl)

clc; clearvars; close all;

% setting random stream
myStream = RandStream('mt19937ar','Seed',10);
RandStream.setGlobalStream(myStream)

%% settings

%Thorax phantom size is fixed to 512
n         = 512; 
phantom   = 'Thorax';

%% load phantom

% load matrices
[A,F,T,Q] = loadSpectralPhantom(phantom);


%% generate spectral measurements

U0    = A*F;

% tomography
n_ang  = 60;
theta  = linspace(0,180,n_ang); 

volGeo = astra_create_vol_geom(n, n);
projGeo= astra_create_proj_geom('parallel', 1, n,theta);
W      = opTomo('cuda', projGeo, volGeo);

% measurements with noise
Y0    = W*U0;
Y     = 0*Y0;
for i=1:size(Y0,2)
   Y(:,i) = astra_add_noise_to_sino(Y0(:,i),Q(i));
end

%% problem setting
% Please note that we only distinguish hard materials from soft materials.
% There are 3 soft materials and 2 hard materials. Hence, we use k = 3
% (2 hard and 1 combined soft) here.

% number of materials
k = 3;

% convert 3 softer ones to 1 soft
[A_Thorax,F_Thorax]  = modifyThorax(A,F);

strProb.A = A_Thorax;
strProb.F = F_Thorax;
strProb.k = k;
strProb.T = T;
strProb.n = [n n];

%% Spectral Tomography algorithms

%%% Reconstruction and then Unmixing
[Aru,Fru] = RU(Y,W,k);

strRU.A   = Aru;
strRU.F   = Fru;

%%% Unmixing then Reconstruction
[Aur,Fur] = UR(Y,W,k);

strUR.A   = Aur;
strUR.F   = Fur;

%%% cJoint
Jopt.rho      = 1e-2;
Jopt.iterMax  = 500;
[Aj,Fj,histJ] = cJoint(Y,W,k,Jopt);

strJ.A = Aj;
strJ.F = Fj;

%%% ADJUST
Dopt.rho        = 1e-2;
Dopt.iterMax    = 10;
[Ad,Fd,Rd,hist] = ADJUST(Y,W,k,T,Dopt);

strD.A = Ad;
strD.F = Fd;

%% plots and tables

strOut = computeResults(strProb,strRU,strUR,strJ,strD);

