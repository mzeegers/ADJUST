%%% Spectral Tomography
%%% Real hyperspectral data
%%% 4 materials 
%%% 128 channels 
%%% Full-view setting: 600 angles from 0 to 360
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

n         = 255;
phantom   = 'Real';

%% reconstruct gold particle

% data for gold particle
datapath = '/export/scratch3/zeegers/SpectralMicroCTDataset/CT_gold_particle/';

darkfield = double(imread(strcat(datapath, 'di000000.tif')));
flatfield = double(imread(strcat(datapath, 'io000000.tif')));
data = zeros(128,257,600);
datacor = zeros(128,257,600);
for i=1:size(data,3)
    data(:,:,i) = imread(strcat(datapath, 'scan_', sprintf('%06d',i-1), '.tif'));
    datacor(:,:,i) = log((flatfield-darkfield)./(data(:,:,i) - darkfield));
end


%% make reconstruction of every channel
n_ang = 600;
theta = linspace(0,2*pi,n_ang); 

DETSIZE = 0.80;
SOD = 12.750527;
SDD = 882.113593;
MAGN = SDD/SOD;
CONV = MAGN/DETSIZE;

volGeo = astra_create_vol_geom(255, 255);
projGeo = astra_create_proj_geom('fanflat', DETSIZE*CONV, 257, theta, SOD*CONV, (SDD-SOD)*CONV);
%W     = opTomo('cuda', projGeo, volGeo);

cor_shift = 257/2-142.373917
projGeo_cor = astra_geom_postalignment(projGeo, cor_shift)

recsGold = zeros(128,255,255);

for i=1:128   
    i
    ALG = 'SIRT_CUDA';
    ITER = 100;

    sinogram_id = astra_mex_data2d('create', '-sino', projGeo_cor, permute(datacor(i,:,:), [3 2 1]));
    rec_id = astra_mex_data2d('create', '-vol', volGeo);
    cfg = astra_struct(ALG);
    cfg.ReconstructionDataId = rec_id;
    cfg.ProjectionDataId = sinogram_id;

    alg_id = astra_mex_algorithm('create', cfg);

    astra_mex_algorithm('iterate', alg_id, ITER);

    recsGold(i,:,:) = astra_mex_data2d('get', rec_id);
    
end

%imshow(rescale(permute(recsGold(50,120:130,90:140), [2 3 1])))
%plot(recsGold(:,126,128))

%% reconstruct particles

% data for gold particle
datapath = '/export/scratch3/zeegers/SpectralMicroCTDataset/Spectral_CT_W_Au_Pb/';

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

cor_shift = 257/2-127.379570
projGeo_cor = astra_geom_postalignment(projGeo, cor_shift)

recsParticles = zeros(128,255,255);

for i=1:128   
    i
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


%imshow(rescale(permute(recsParticles(20,10:30,150:160), [2 3 1])))
%plot(recsParticles(:,24,155)) %gold
%plot(recsParticles(:,109,122)) %thungsten
%plot(recsParticles(:,95,160)) %lead, but not lead, I think it is gold
%plot(recsParticles(:,95,200)) %quartz
%plot(recsParticles(:,207,106)) %lead 1
%plot(recsParticles(:,33,160)) %lead 2
%plot(recsParticles(:,134,40)) %lead 3  = 1.2 lead3
%plot(recsParticles(:,96,184)) %lead 4

%%

k = 2;

n_ang  = 600;
theta  = linspace(0,2*pi,n_ang); 

DETSIZE = 0.80;
SOD = 57.999001;
SDD = 883.081116;
MAGN = SDD/SOD;
CONV = MAGN/DETSIZE;

volGeo = astra_create_vol_geom(255, 255);
projGeo = astra_create_proj_geom('fanflat', DETSIZE*CONV, 257, theta, SOD*CONV, (SDD-SOD)*CONV);

cor_shift = 257/2-127.379570
projGeo_cor = astra_geom_postalignment(projGeo, cor_shift)
W     = opTomo('cuda', projGeo_cor, volGeo);

%Y = permute(datacor, [3 2 1]);

%% adjustments to the data

recsParticles2 = recsParticles;
for i=1:255
    for j=1:255
        recsParticles2(:,i,j) = mean(recsParticles(:,130:140,130:140),[2,3]);
    end    
end
recsParticles2(:,200:210,100:110) = recsParticles(:,200:210,100:110); %lead 1
recsParticles2(:,30:40,155:165) = recsParticles(:,30:40,155:165); %lead 2
recsParticles2(:,130:140,45:55) = recsParticles(:,130:140,45:55); %lead 3  = 1.2 lead3
recsParticles2(:,90:100,180:190) = recsParticles(:,90:100,180:190); %lead 4
%%




%Y = reshape(Y, [] , size(Y,3));

n = 255;

T = zeros(3,128);

%T(1,:) = recsParticles(:,24,155); %gold
%T(2,:) = recsParticles(:,109,122); %thungsten
T(1,:) = recsParticles(:,95,200); %quartz
T(2,:) = recsParticles(:,207,106); %lead

Ynew = zeros(128,257*600);
for i=1:128
    Ynew(i,:) = W*reshape(recsParticles2(i,:,:), 255*255,1);
end

%%


%% Spectral Tomography algorithms

%%% Reconstruction and then Unmixing
[Aru,Fru] = RU(Y,W,k);

strRU.A   = Aru;
strRU.F   = Fru;

%%

%%% Unmixing then Reconstruction
[Aur,Fur] = UR(Y,W,k);

strUR.A   = Aur;
strUR.F   = Fur;

%%
%%% cJoint
Jopt.rho      = 1e-2;
Jopt.iterMax  = 1;
[Aj,Fj,histJ] = cJoint(Y,W,k,Jopt);

strJ.A = Aj;
strJ.F = Fj;

%% ADJUST
Dopt.rho        = 1e-2;
Dopt.iterMax    = 1000;
[Ad,Fd,Rd,hist] = ADJUST(Y,W,k,T,Dopt);

strD.A = Ad;
strD.F = Fd;

%%
strProb.A = Ad;
strProb.F = T;
strProb.k = k;
strProb.T = T;
strProb.n = [n n];

%% plots and tables

strOut = computeResults(strProb,strRU,strUR,strJ,strD);
save('results/example_realDatak4initial.mat', 'strOut');
