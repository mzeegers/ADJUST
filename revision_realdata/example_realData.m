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
%plot(recsParticles(:,134,40)) %lead 3  = 1.2*lead3
%plot(recsParticles(:,96,184)) %lead 4

%%

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

%% adjustments to the data

recsParticles2 = recsParticles;
% for i=1:255
%     for j=1:255
%         recsParticles2(:,i,j) = mean(recsParticles(:,130:140,130:140),[2,3]);
%     end    
% end

% recsParticles2(:,200:210,100:110) = recsParticles(:,200:210,100:110); %lead 1
% recsParticles2(:,30:40,155:165) = recsParticles(:,30:40,155:165); %lead 2
% recsParticles2(:,130:140,35:45) = recsParticles(:,130:140,35:45); %lead 3
% %= 1.2 lead3 , werkt deze wel?
% recsParticles2(:,90:100,180:190) = recsParticles(:,90:100,180:190); %lead 4




%%
%Particle locations:
% 18+1:23+1,138+1:144+1 % tungsten
% 49+1:51+1,143+1:147+1 % tungsten
% 70+1:74+1,127+1:132+1 % tungsten
% 105+1:110+1,117+1:123+1 % tungsten
% 91+1:95+1,157+1:161+1 % tungsten
% 105+1:111+1,164+1:169+1 % tungsten
% 139+1:144+1,115+1:120+1 % tungsten
% 151+1:157+1,215+1:220+1 % tungsten
% 190+1:195+1,167+1:171+1 % tungsten
% 202+1:206+1,164+1:168+1 % tungsten
% 195+1:201+1,210+1:216+1 % tungsten
% 219+1:224+1,179+1:185+1 % tungsten
% 224+1:228+1,136+1:141+1 % tungsten
% 179+1:182+1,82+1:85+1 % tungsten
% 186+1:189+1,33+1:37+1 % tungsten
% 164+1:168+1,107+1:112+1 % tungsten
% 29+1:33+1,127+1:132+1 % tungsten

% 19+1:25+1,151+1:156+1 % gold
% 73+1:78+1,141+1:145+1 % gold
% 113+1:118+1,145+1:149+1 % gold
% 134+1:138+1,108+1:114+1 % gold
% 145+1:149+1,147+1:152+1 % gold
% 209+1:214+1,136+1:141+1 % gold
% 222+1:224+1,132+1:135+1 % gold
% 189+1:192+1,112+1:117+1 % gold
% 191+1:196+1,90+1:95+1 % gold
% 180+1:185+1,31+1:37+1 % gold
% 152+1:156+1,29+1:34+1 % gold 
% 149+1:153+1,17+1:20+1 % gold
% 117+1:122+1,46+1:51+1 % gold
% 93+1:97+1,34+1:37+1 % gold
% 84+1:89+1,25+1:28+1 % gold
% 58+1:66+1,33+1:38+1 % gold
% 40+1:46+1,55+1:60+1 % gold
% 42+1:46+1,62+1:67+1 % gold
% 41+1:44+1,114+1:121+1 % gold

% 29+1:33+1,157+1:160+1 % lead
% 92+1:96+1,180+1:185+1 % lead
% 167+1:172+1,168+1:173+1 % lead
% 217+1:220+1,156+1:160+1 % lead
% 202+1:207+1,102+1:106+1 % lead
% 155+1:158+1,121+1:123+1 % lead
% 130+1:134+1,37+1:41+1 % lead
% 61+1:65+1,88+1:91+1 % lead
% 52+1:56+1,89+1:102+1 % lead
% 20+1:23+1,117+1:120+1 % lead


%Hide the lead
% for i=29+1:33+2
%     for j = 157+1:160+2
%         recsParticles2(:,i,j) = mean(recsParticles(:,130:140,130:140),[2,3]);
%     end
% end
% for i=92+1:96+2
%     for j = 180+1:185+2
%         recsParticles2(:,i,j) = mean(recsParticles(:,130:140,130:140),[2,3]);
%     end
% end
% for i=167+1:172+2
%     for j = 168+1:173+2
%         recsParticles2(:,i,j) = mean(recsParticles(:,130:140,130:140),[2,3]);
%     end
% end
% for i=217+1:220+2
%     for j = 156+1:160+2
%         recsParticles2(:,i,j) = mean(recsParticles(:,130:140,130:140),[2,3]);
%     end
% end
% for i= 202+1:207+2
%     for j = 102+1:106+2
%         recsParticles2(:,i,j) = mean(recsParticles(:,130:140,130:140),[2,3]);
%     end
% end
% for i=155+1:158+2
%     for j = 121+1:123+2
%         recsParticles2(:,i,j) = mean(recsParticles(:,130:140,130:140),[2,3]);
%     end
% end
% for i=130+1:134+2
%     for j = 37+1:41+2
%         recsParticles2(:,i,j) = mean(recsParticles(:,130:140,130:140),[2,3]);
%     end
% end
% for i=61+1:65+2
%     for j = 88+1:91+2
%         recsParticles2(:,i,j) = mean(recsParticles(:,130:140,130:140),[2,3]);
%     end
% end
% for i=52+1:56+2
%     for j = 89+1:102+2
%         recsParticles2(:,i,j) = mean(recsParticles(:,130:140,130:140),[2,3]);
%     end
% end
% for i=20+1:23+2
%     for j = 117+1:120+2
%         recsParticles2(:,i,j) = mean(recsParticles(:,130:140,130:140),[2,3]);
%     end
% end

%Hide the gold
% for i=19+1:25+2
%     for j = 151+1:156+2
%         recsParticles2(:,i,j) = mean(recsParticles(:,130:140,130:140),[2,3]);
%     end
% end
% for i=73+1:78+2
%     for j = 141+1:145+2
%         recsParticles2(:,i,j) = mean(recsParticles(:,130:140,130:140),[2,3]);
%     end
% end
% for i=113+1:118+2
%     for j = 145+1:149+2
%         recsParticles2(:,i,j) = mean(recsParticles(:,130:140,130:140),[2,3]);
%     end
% end
% for i=134+1:138+2
%     for j = 108+1:114+2
%         recsParticles2(:,i,j) = mean(recsParticles(:,130:140,130:140),[2,3]);
%     end
% end
% for i=145+1:149+2
%     for j = 147+1:152+2
%         recsParticles2(:,i,j) = mean(recsParticles(:,130:140,130:140),[2,3]);
%     end
% end
% for i=209+1:214+2
%     for j = 136+1:141+2
%         recsParticles2(:,i,j) = mean(recsParticles(:,130:140,130:140),[2,3]);
%     end
% end
% for i=222+1:224+2
%     for j = 132+1:135+2
%         recsParticles2(:,i,j) = mean(recsParticles(:,130:140,130:140),[2,3]);
%     end
% end
% for i=189+1:192+2
%     for j = 112+1:117+2
%         recsParticles2(:,i,j) = mean(recsParticles(:,130:140,130:140),[2,3]);
%     end
% end
% for i=191+1:196+2
%     for j = 90+1:95+2
%         recsParticles2(:,i,j) = mean(recsParticles(:,130:140,130:140),[2,3]);
%     end
% end
% for i=180+1:185+2
%     for j = 31+1:37+2
%         recsParticles2(:,i,j) = mean(recsParticles(:,130:140,130:140),[2,3]);
%     end
% end
% for i=152+1:156+2
%     for j = 29+1:34+2
%         recsParticles2(:,i,j) = mean(recsParticles(:,130:140,130:140),[2,3]);
%     end
% end
% for i=149+1:153+2
%     for j = 17+1:20+2
%         recsParticles2(:,i,j) = mean(recsParticles(:,130:140,130:140),[2,3]);
%     end
% end
% for i=117+1:122+2
%     for j = 46+1:51+2
%         recsParticles2(:,i,j) = mean(recsParticles(:,130:140,130:140),[2,3]);
%     end
% end
% for i=93+1:97+2
%     for j = 34+1:37+2
%         recsParticles2(:,i,j) = mean(recsParticles(:,130:140,130:140),[2,3]);
%     end
% end
% for i=84+1:89+2
%     for j = 25+1:28+2
%         recsParticles2(:,i,j) = mean(recsParticles(:,130:140,130:140),[2,3]);
%     end
% end
% for i=58+1:66+2
%     for j = 33+1:38+2
%         recsParticles2(:,i,j) = mean(recsParticles(:,130:140,130:140),[2,3]);
%     end
% end
% for i=40+1:46+2
%     for j = 55+1:60+2
%         recsParticles2(:,i,j) = mean(recsParticles(:,130:140,130:140),[2,3]);
%     end
% end
% for i=42+1:46+2
%     for j = 62+1:67+2
%         recsParticles2(:,i,j) = mean(recsParticles(:,130:140,130:140),[2,3]);
%     end
% end
% for i=41+1:44+2
%     for j = 114+1:121+2
%         recsParticles2(:,i,j) = mean(recsParticles(:,130:140,130:140),[2,3]);
%     end
% end

% Hide the tungsten
% for i=18+1:23+2
%     for j = 138+1:144+2
%         recsParticles2(:,i,j) = mean(recsParticles(:,130:140,130:140),[2,3]);
%     end
% end
% for i=49+1:51+2
%     for j = 143+1:147+2
%         recsParticles2(:,i,j) = mean(recsParticles(:,130:140,130:140),[2,3]);
%     end
% end
% for i=70+1:74+2
%     for j = 127+1:132+2
%         recsParticles2(:,i,j) = mean(recsParticles(:,130:140,130:140),[2,3]);
%     end
% end
% for i=105+1:110+2
%     for j = 117+1:123+2
%         recsParticles2(:,i,j) = mean(recsParticles(:,130:140,130:140),[2,3]);
%     end
% end
% for i=91+1:95+2
%     for j = 57+1:161+2
%         recsParticles2(:,i,j) = mean(recsParticles(:,130:140,130:140),[2,3]);
%     end
% end
% for i=105+1:111+2
%     for j = 164+1:169+2
%         recsParticles2(:,i,j) = mean(recsParticles(:,130:140,130:140),[2,3]);
%     end
% end
% for i=139+1:144+2
%     for j = 115+1:120+2
%         recsParticles2(:,i,j) = mean(recsParticles(:,130:140,130:140),[2,3]);
%     end
% end
% for i=151+1:157+2
%     for j = 215+1:220+2
%         recsParticles2(:,i,j) = mean(recsParticles(:,130:140,130:140),[2,3]);
%     end
% end
% for i=190+1:195+2
%     for j = 167+1:171+2
%         recsParticles2(:,i,j) = mean(recsParticles(:,130:140,130:140),[2,3]);
%     end
% end
% for i=202+1:206+2
%     for j = 164+1:168+2
%         recsParticles2(:,i,j) = mean(recsParticles(:,130:140,130:140),[2,3]);
%     end
% end
% for i=195+1:201+2
%     for j = 210+1:216+2
%         recsParticles2(:,i,j) = mean(recsParticles(:,130:140,130:140),[2,3]);
%     end
% end
% for i=219+1:224+2
%     for j = 179+1:185+2
%         recsParticles2(:,i,j) = mean(recsParticles(:,130:140,130:140),[2,3]);
%     end
% end
% for i=224+1:228+2
%     for j = 136+1:141+1
%         recsParticles2(:,i,j) = mean(recsParticles(:,130:140,130:140),[2,3]);
%     end
% end
% for i=179+1:182+2
%     for j = 82+1:85+2
%         recsParticles2(:,i,j) = mean(recsParticles(:,130:140,130:140),[2,3]);
%     end
% end
% for i=186+1:189+2
%     for j = 33+1:37+2
%         recsParticles2(:,i,j) = mean(recsParticles(:,130:140,130:140),[2,3]);
%     end
% end
% for i=164+1:168+2
%     for j = 107+1:112+2
%         recsParticles2(:,i,j) = mean(recsParticles(:,130:140,130:140),[2,3]);
%     end
% end
% for i=29+1:33+2
%     for j = 127+1:132+2
%         recsParticles2(:,i,j) = mean(recsParticles(:,130:140,130:140),[2,3]);
%     end
% end


%%
Y2 = zeros(154200,128);
for i=1:128
    Y2(:,i) = W*reshape(recsParticles2(i,:,:),255*255,1);
end


%%
Y = permute(datacor, [3 2 1]);
Y = reshape(Y, [] , size(Y,3));

%%
k = 4;
n = 255;
T = zeros(k,128);

T(1,:) = recsParticles(:,24,155); %gold
T(2,:) = recsParticles(:,109,122); %tungsten
T(3,:) = mean(recsParticles(:,130:140,130:140),[2,3]); %quartz
T(4,:) = recsParticles(:,207,106); %recsParticles(:,96,184); %lead

%writematrix(T,'SpectraExtracted.csv')

%Combination gold and lead is fine, combination gold and tungsten is fine
%Combination tungsten and lead

Y0 = Y;
T0 = T;

Y = Y0(:,30:94);
T = T0(:,30:94);


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
Dopt.iterMax    = 10000;
Dopt.relTol     = 1e-11
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
