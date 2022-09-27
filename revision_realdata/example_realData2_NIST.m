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

globalDataDir = [pwd '/data/'];

%% settings

n         = 255;
phantom   = 'Real';

%% reconstruct particles

% data for gold particle
datapath = [globalDataDir 'Spectral_CT_W_Au_Pb/'];

darkfield = double(imread(strcat(datapath, 'di000000.tif')));
flatfield = double(imread(strcat(datapath, 'io000000.tif')));
data = zeros(128,257,600);
datacor = zeros(128,257,600);
for i=1:size(data,3)
    data(:,:,i) = imread(strcat(datapath, 'scan_', sprintf('%06d',i-1), '.tif'));
    datacor(:,:,i) = log((flatfield-darkfield)./(data(:,:,i) - darkfield));
end

data = flipdim(data,3);
datacor = flipdim(datacor,3);

%% make reconstruction of every channel
n_ang = 600;
theta = linspace(0,2*pi,n_ang); 

DETSIZE = 0.80;
SOD = 57.999001;
SDD = 883.081116;
MAGN = SDD/SOD;
CONV = MAGN/DETSIZE;

volGeo = astra_create_vol_geom(255, 255);
projGeo = astra_create_proj_geom('fanflat', DETSIZE*CONV, 257, theta, SOD*CONV, (SDD-SOD)*CONV);
%W     = opTomo('cuda', projGeo, volGeo);

cor_shift = 257/2-127.379570;
projGeo_cor = astra_geom_postalignment(projGeo, cor_shift);

recsParticles = zeros(128,255,255);

for i=1:128   
    fprintf('running spectral idx=%d \n',i);
    ALG = 'SIRT_CUDA';
    ITER = 100;

    sinogram_id = astra_mex_data2d('create', '-sino', projGeo_cor, permute(datacor(i,:,:), [3 2 1]));
    rec_id = astra_mex_data2d('create', '-vol', volGeo);
    cfg = astra_struct(ALG);
    cfg.ReconstructionDataId = rec_id;
    cfg.ProjectionDataId = sinogram_id;

    alg_id = astra_mex_algorithm('create', cfg);

    astra_mex_algorithm('iterate', alg_id, ITER);

    recsParticles(i,:,:) = astra_mex_data2d('get', rec_id);
    
end


%% forward operator

n_ang  = 600;
theta  = linspace(0,2*pi,n_ang); 

DETSIZE = 0.80;
SOD = 57.999001;
SDD = 883.081116;
MAGN = SDD/SOD;
CONV = MAGN/DETSIZE;

volGeo = astra_create_vol_geom(255, 255);
projGeo = astra_create_proj_geom('fanflat', DETSIZE*CONV, 257, theta, SOD*CONV, (SDD-SOD)*CONV);

cor_shift = 257/2-127.379570;
projGeo_cor = astra_geom_postalignment(projGeo, cor_shift);
W     = opTomo('cuda', projGeo_cor, volGeo);

Y = permute(datacor, [3 2 1]);
Y = reshape(Y, [] , size(Y,3));

n = 255;

%% dictionary

T1 = readmatrix([globalDataDir 'NISTequivalentofRealData_MatrixTXray128chan_3_MinEn20.35_MaxEn161.15_mat_CORRECTED.csv']);
T1 = 0.01*T1;
T1(4,:) = recsParticles(:,95,200); %quartz

T0 = zeros(3,128);
T0(1,:) = recsParticles(:,24,155); %gold
T0(2,:) = recsParticles(:,109,122); %thungsten
T0(3,:) = recsParticles(:,207,106); %lead
T0(4,:) = recsParticles(:,95,200); %quartz


Y0 = Y;
%%

k = 4;

Y = Y0(:,30:94);
T = T1(:,30:94);
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
Jopt.iterMax  = 50;
Jopt.showIter = 2;
[Aj,Fj,histJ] = cJoint(Y,W,k,Jopt);

strJ.A = Aj;
strJ.F = Fj;

%% ADJUST

Dopt.rho        = 1e-2;
Dopt.iterMax    = 1000;
Dopt.innerIter  = 2;
Dopt.showIter   = 2;
[Ad,Fd,Rd,hist] = ADJUST(Y,W,k,T,Dopt);

strD.A = Ad;
strD.F = Fd;

%% plots and tables

strProb.A = Ad;
strProb.F = Fd;

strOut = computeResults(strProb,strRU,strUR,strJ,strD);
save([pwd '/results/example_realData2.mat'], 'strOut');
