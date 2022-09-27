%%% X-ray test script (with dataset for F and T)
%%% Dictionary approach
%%% Shepp Logan phantom of size 128 x 128
%%% 6 materials 
%%% 100 channels 
%%% Full-view setting: 60 angles from 0 to 180


%clc;
clearvars; close all;

% setting random stream
myStream = RandStream('mt19937ar','Seed',10);
RandStream.setGlobalStream(myStream)

%% generate U

% load matrices
phantomId = 5;
n = 512;
%Thard = readmatrix('MatrixTXrayHard.csv');
%Tsoft = readmatrix('MatrixTXraySoft.csv');
F     = readmatrix('Thorax_MatrixFXray100chan_17matBoneBloodIodineSofttissueBloodLungsIodine_20to80kV_CORRECTED.csv');
T     = readmatrix('Thorax_MatrixTXray100chan_60mat_20to80kV_CORRECTED.csv');
Q     = readmatrix('VectorIXrayMatWEn100Min20Max80.csv');

% adding Iodine
T     = [T;F(2,:)];
% T     = transpose(licols(T'));

[A,F] = loadPhantomrev(phantomId,n,T,F);

% scale matrix A
A = A(:,1:size(F,1));
A = A/max(A(:));

% create matrices for measurements
Am = zeros(size(A,1)*4, size(A,2));
for i=1:size(A,2) 
    Am(:,i) = reshape(imresize(reshape(A(:,i),n,n),2, 'nearest'), n*n*4, 1);
end

% starts from 8 kV
kC = 1;
F = F(:,kC:end);
T = T(:,kC:end);

% generate U
U0    = Am*F/2;

% tomography
%N     = sqrt(size(A,1));
n_ang = 180;
theta = linspace(0,180,n_ang); 

volGeo= astra_create_vol_geom(n, n);
projGeo= astra_create_proj_geom('parallel', 1, n,theta);
W     = opTomo('cuda', projGeo, volGeo);

volGeo = astra_create_vol_geom(2*n, 2*n);
projGeo= astra_create_proj_geom('parallel', 2, n,theta);
Wm     = opTomo('cuda', projGeo, volGeo);

% measurements with noise
Y0    = Wm*U0;
Y     = 0*Y0;
for i=1:size(Y0,2)
   Y(:,i) = astra_add_noise_to_sino(Y0(:,i),Q(i));
end

% sizes
[m,c] = size(U0);          % number of images x number of channels
k = 3;

strProb.A = Am;
strProb.F = F;
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

%%
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

A1 = [A(:,1) A(:,2) A(:,3)+A(:,4)+A(:,5)];
strProb.A = A1;

strOut = computeResults(strProb,strRU,strUR,strJ,strD);
save('results/example_Thorax_nonInversCrime.mat', 'strOut');



