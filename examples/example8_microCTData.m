%%% Spectral Tomography
%%% MicroCT data of size 1 x 257
%%% 4 materials 
%%% 64 channels (reduced from 128)
%%% Full-view setting: 600 angles from 0 to 360
%%% Data acquired by Jonathan Sittner et al. (2022)
%%% Data taken from: https://rodare.hzdr.de/record/1627
%%% More information: https://doi.org/10.1002/xrs.3200
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

% directory for saving results
saveDir = [pwd '\results\'];
if ~exist(saveDir,'dir'), mkdir(saveDir); end

%% settings

n         = 255;            % spatial resolution the data
phantom   = 'Real';

%% data preprocessing

% data with particles - change the path to location of the data
datapath = '/Spectral_CT_W_Au_Pb/';

% dark and flatfield correction
darkfield = double(imread(strcat(datapath, 'di000000.tif')));
flatfield = double(imread(strcat(datapath, 'io000000.tif')));
data = zeros(128,257,600);
datacor = zeros(128,257,600);
for i=1:size(data,3)
    data(:,:,i) = imread(strcat(datapath, 'scan_', sprintf('%06d',i-1), '.tif'));
    datacor(:,:,i) = log((flatfield-darkfield)./(data(:,:,i) - darkfield));
end

% flip the sinograms for correct alignment with angles
data = flipdim(data,3);
datacor = flipdim(datacor,3);

%% make reconstruction of every spectral channel
n_ang = 600;
theta = linspace(0,2*pi,n_ang); 

% values taken from the settings file
DETSIZE = 0.80;
SOD = 57.999001;
SDD = 883.081116;
MAGN = SDD/SOD;
CONV = MAGN/DETSIZE;

volGeo = astra_create_vol_geom(255, 255);
projGeo = astra_create_proj_geom('fanflat', DETSIZE*CONV, 257, theta, SOD*CONV, (SDD-SOD)*CONV);

cor_shift = 257/2-127.379570;
projGeo_cor = astra_geom_postalignment(projGeo, cor_shift);

recsParticles = zeros(128,255,255);

ALG = 'SIRT_CUDA';
ITER = 100;

for i=1:128   
    fprintf('running spectral idx=%d \n',i);

    sinogram_id = astra_mex_data2d('create', '-sino', projGeo_cor, permute(datacor(i,:,:), [3 2 1]));
    rec_id = astra_mex_data2d('create', '-vol', volGeo);
    cfg = astra_struct(ALG);
    cfg.ReconstructionDataId = rec_id;
    cfg.ProjectionDataId = sinogram_id;

    alg_id = astra_mex_algorithm('create', cfg);

    astra_mex_algorithm('iterate', alg_id, ITER);

    recsParticles(i,:,:) = astra_mex_data2d('get', rec_id);
    
end


%% forward operator to be used for inversion

W     = opTomo('cuda', projGeo_cor, volGeo);

Y0 = permute(datacor, [3 2 1]);
Y0 = reshape(Y0, [] , size(Y0,3));

%% load and adjust data

T1 = readmatrix('NISTequivalentofRealData_MatrixTXray128chan_42_MinEn20.35_MaxEn161.15.csv');
T1 = 0.015*T1;

% add quartz extracted from data to dictionary
T1(43,:) = mean(recsParticles(:,130:140,130:140),[2,3]);

% number of materials (to be reconstructed)
k = 4;

%restrict data and dictionary to 64 spectral bins
Y = Y0(:,30:94);
T = T1(:,30:94);
T_gt = T1(:,30:94);

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

%% cJoint
Jopt.rho      = 1e-2;
Jopt.iterMax  = 500;
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

%% ADJUST-true

Dopt.rho        = 1e-2;
Dopt.iterMax    = 1000;
Dopt.innerIter  = 2;
Dopt.showIter   = 2;
[Ad_gt,Fd_gt,Rd_gt,hist_gt] = ADJUST(Y,W,k,T_gt,Dopt);

strProb.A = Ad_gt;
strProb.F = Fd_gt;

%% plots and tables

strOut = computeResults(strProb,strRU,strUR,strJ,strD);

%% save results

save([saveDir 'example_microCTData.mat'],'strOut');