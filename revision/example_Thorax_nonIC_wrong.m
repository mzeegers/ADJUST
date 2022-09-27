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

% load matrices for inversion
[A_orig,F_orig,T,Q] = loadSpectralPhantom(phantom, n); % n does not do anything here
[A,F]  = modifyThorax(A_orig,F_orig);

% create matrices for measurements
Am = zeros(size(A,1)*4, size(A,2));
for i=1:size(A,2) 
    Am(:,i) = reshape(imresize(reshape(A(:,i),n,n),2), n*n*4, 1);
end

Fm = F;
Tm = T;
Qm = Q;

% Please note that we only distinguish hard materials from soft materials.
% There are 3 soft materials and 2 hard materials. Hence, we use k = 3
% (2 hard and 1 combined soft) here.
k = 3;

strProb.A = A;
strProb.F = F;
strProb.k = k;
strProb.T = T;
strProb.n = [n n];
        
%% generate spectral measurements

U0    = Am*Fm/2;

% tomography
n_ang  = 180;
theta  = linspace(0,180,n_ang); 

volGeo = astra_create_vol_geom(2*n, 2*n);
projGeo= astra_create_proj_geom('parallel', 2, n,theta);
Wm     = opTomo('cuda', projGeo, volGeo);

% measurements with noise
Y0    = Wm*U0;
Y     = 0*Y0;
for i=1:size(Y0,2)
   Y(:,i) = astra_add_noise_to_sino(Y0(:,i),Qm(i));
end

%% (non-inverse-crime) forward operator to be used for inversion

volGeo = astra_create_vol_geom(n, n);
projGeo= astra_create_proj_geom('parallel', 1, n,theta);
W      = opTomo('cuda', projGeo, volGeo);

% (only for checking) measurments with forward operator
U1    = A*F;
Y1    = W*U1;

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
Jopt.iterMax  = 2000;
[Aj,Fj,histJ] = cJoint(Y,W,k,Jopt);

strJ.A = Aj;
strJ.F = Fj;
%%
%%% ADJUST
Dopt.rho        = 1e-2;
Dopt.iterMax    = 500;
[Ad,Fd,Rd,hist] = ADJUST(Y,W,k,T,Dopt);

strD.A = Ad;
strD.F = Fd;

%% plots and tables

strOut = computeResults(strProb,strRU,strUR,strJ,strD);
save('results/example_Thorax_nonInversCrime.mat', 'strOut');
