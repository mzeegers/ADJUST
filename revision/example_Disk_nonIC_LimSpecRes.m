%%% Spectral Tomography
%%% Disk phantom of size 256 x 256
%%% 8 materials 
%%% 100 channels 
%%% Full-view setting: 60 angles from 0 to 180
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

n         = 512;
phantom   = 'Disks';

%% load phantom

% load matrices for measurements
[Am,Fm,Tm,Qm] = loadSpectralPhantom(phantom,2*n);

% load matrices for inversion
[A,F,T,Q] = loadSpectralPhantom(phantom,n);

k = size(A,2);

%% generate spectral measurements

U0    = Am*Fm/2;

% tomography
n_ang  = 180;
theta  = linspace(0,pi,n_ang); 

volGeo = astra_create_vol_geom(2*n, 2*n);
projGeo= astra_create_proj_geom('parallel', 2, n,theta);
Wm     = opTomo('cuda', projGeo, volGeo);

% measurements with noise
Y0    = Wm*U0;

%Resize the data and the dictionary
%Resize energy bins to 10 in total
Yfullres = Y0;
Tfullres = T;
Ffullres = F;
Qfullres = Qm;

%Y = imresize(Y,[size(Y,1) size(Y,2)/10], 'nearest');
%T = imresize(T,[size(T,1) size(T,2)/10], 'nearest');
%F = imresize(F,[size(F,1) size(F,2)/10], 'nearest');

AVG_COLS = 10;
% Dimension over which to average
DIM = 2; % Columns
%F2 = Ffullres;
% Use filter to calculate the moving average across EVERY combination of columns
F_moving_avg = filter(ones(1,AVG_COLS)/AVG_COLS,1,Ffullres,[],DIM);
Y0_moving_avg = filter(ones(1,AVG_COLS)/AVG_COLS,1,Y0,[],DIM);
T_moving_avg = filter(ones(1,AVG_COLS)/AVG_COLS,1,T,[],DIM);
Q_moving_avg = filter(ones(1,AVG_COLS)/AVG_COLS,1,Qm,[],1);
% Grab only the column averages that were actually wanted
F = F_moving_avg(:,AVG_COLS:AVG_COLS:end);
Y0 = Y0_moving_avg(:,AVG_COLS:AVG_COLS:end);
T = T_moving_avg(:,AVG_COLS:AVG_COLS:end);
Qm = Q_moving_avg(AVG_COLS:AVG_COLS:end);

Y     = 0*Y0;
for i=1:size(Y0,2)
   Y(:,i) = astra_add_noise_to_sino(Y0(:,i),Qm(i));
end

strProb.A = A;
strProb.F = F;
strProb.Ffullres = Ffullres;
strProb.k = k;
strProb.T = T;
strProb.Qfullres = Qfullres;
strProb.Tfullres = Tfullres;
strProb.Yfullres = Yfullres;
strProb.n = [n n];

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

%% ADJUST
Dopt.rho        = 1e-2;
Dopt.iterMax    = 500;
[Ad,Fd,Rd,hist] = ADJUST(Y,W,k,T,Dopt);

strD.A = Ad;
strD.F = Fd;

%% plots and tables

strOut = computeResults(strProb,strRU,strUR,strJ,strD);
save('results/example_Disk_nonInversCrime_LimSpecRes.mat', 'strOut');
