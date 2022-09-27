%%% X-ray test script (with dataset for F and T)
%%% Dictionary approach
%%% Shepp Logan phantom of size 128 x 128
%%% 6 materials 
%%% 100 channels 
%%% Full-view setting: 60 angles from 0 to 180


%clc;
clearvars; close all;

% add required folders
mydir  = pwd;
idcs   = strfind(mydir,'/');
rwd    = mydir(1:idcs(end)-1);
addpath(genpath([rwd '/mbin/']));
addpath(genpath([rwd '/Data/']));
addpath(genpath([rwd '/GeneratedData/']));

% Get current time in year-month-day-hour-minute-sec format
tfmt = 'yyyy-MM-dd_HH';
startTime = strrep(strrep(datestr(datetime('now', 'InputFormat', tfmt, 'Format', tfmt)), ' ', '_'), ':', '-');

% results directory
phantomId = 5;
resDir = [rwd '/results/ex4_Thorax_20to80kV' startTime '/'];
if ~exist(resDir,'dir'), mkdir(resDir); end

% setting random stream
myStream = RandStream('mt19937ar','Seed',10);
RandStream.setGlobalStream(myStream)

%% generate U

% load matrices
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
    Am(:,i) = reshape(imresize(reshape(A(:,i),n,n),2), n*n*4, 1);
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
sigma = norm(Y-Y0,'fro')/norm(Y0,'fro');

% condition T (remove dependent rows)
% T = transpose(licols(T',1e-9));

% sizes
[m,c] = size(U0);          % number of images x number of channels
k     = 3; % size(A,2);     % number of materials


fprintf('cond(Ftrue) = %.2f \n',cond(F));
Su = svd(U0);
fprintf('cond(U)     = %.2f \n',Su(1)/Su(k));
fprintf('size(U): %d x %d | k = %d \n',size(U0),k);
fprintf('size(W): %d x %d \n',size(W));
fprintf('size(T): %d x %d \n',size(T));

%% material images and their signatures

for i=1:size(A,2)
    fig1 = figure(1);
    imagesc(reshape(A(:,i),n,n));
    axis square;colormap gray;set(gca,'xtick',[],'ytick',[]);
    saveas(fig1,[resDir 'At-' num2str(i)],'epsc');
    saveas(fig1,[resDir 'At-' num2str(i)],'fig');
    export_fig Test.eps
end

fig2 = figure(2);
semilogy(F','LineWidth',2);pbaspect([3 1 1]);set(gca,'FontSize',16);
saveas(fig2,[resDir 'Ft'],'epsc');
saveas(fig2,[resDir 'Ft'],'fig');

%Save matrices
save([resDir 'At.mat'],'A');
save([resDir 'Ft.mat'],'F');

%% Reconstruction and then Unmixing

fprintf('******* reconstruction-then-unmixing \n')
Wty = W'*Y;
WtW = @(x) W'*(W*x) + 1e-3*x;

for i=1:size(Y,2)
   [Up(:,i),~] = pcg(WtW,Wty(:,i),1e-6,20);
end

Up(Up <0) = 0;

opt = statset('MaxIter',10,'Display','final');
[A10,F10] = nnmf(Up,k,'Replicates',10,'Options',opt,'Algorithm','mult','Replicates',10);

opt = statset('MaxIter',1000,'Display','final');
[Aru,Fru] = nnmf(Up,k,'W0',A10,'H0',F10,'Options',opt,'Algorithm','als','Replicates',10);

%%% material images and their signatures
AruMax = max(abs(Aru(:)));
for i=1:size(Aru,2)
    fig1 = figure(1);
    imagesc(reshape(Aru(:,i),n,n),[-AruMax AruMax]);
    axis square;colormap gray;set(gca,'xtick',[],'ytick',[]);
    saveas(fig1,[resDir 'Aru-' num2str(i)],'epsc');
    saveas(fig1,[resDir 'Aru-' num2str(i)],'fig');
end

figure(2);
plot(Fru','LineWidth',2);pbaspect([3 1 1]);set(gca,'FontSize',16);
saveas(fig2,[resDir 'Fru'],'epsc');
saveas(fig2,[resDir 'Fru'],'fig');

%Save matrices
save([resDir 'Aru.mat'],'Aru');
save([resDir 'Fru.mat'],'Fru');

%% Unmixing then Reconstruction

fprintf('******* unmixing-then-reconstruction \n')
opt = statset('MaxIter',10,'Display','final');
[A20,F20] = nnmf(Y,k,'Replicates',10,'Options',opt,'Algorithm','mult','Replicates',10);

opt = statset('MaxIter',1000,'Display','final');
[WAur,Fur] = nnmf(Y,k,'W0',A20,'H0',F20,'Options',opt,'Algorithm','als','Replicates',10);


Wtp = W'*WAur;
WtW = @(x) W'*(W*x) + 1e-3*x;

for i=1:size(WAur,2)
   [Aur(:,i),~] = pcg(WtW,Wtp(:,i),1e-6,20);
end

%%% material images and their signatures
AurMax = max(abs(Aur(:)));
for i=1:size(Aur,2)
    fig1 = figure(1);
    imagesc(reshape(Aur(:,i),n,n),[-AurMax AurMax]);
    axis square;colormap gray;set(gca,'xtick',[],'ytick',[]);
    saveas(fig1,[resDir 'Aur-' num2str(i)],'epsc');
    saveas(fig1,[resDir 'Aur-' num2str(i)],'fig');
end

fig2 = figure(2);plot(Fur','LineWidth',2);pbaspect([3 1 1]);set(gca,'FontSize',16);
saveas(fig2,[resDir 'Fur'],'epsc');
saveas(fig2,[resDir 'Fur'],'fig');

%Save matrices
save([resDir 'Aur.mat'],'Aur');
save([resDir 'Fur.mat'],'Fur');

%% joint - JUST

fprintf('******* joint\n')

kD = k;
A0 = rand(size(W,2),kD);
F0 = rand(kD,size(Y,2));

Jopt.rho     = 1e-2;
Jopt.iterMax = 500;
Jopt.simFlag = 1;
Jopt.EQ_A    = 0;
Jopt.lambda  = 0;

[Aj,Fj,histJ] = JUST(Y,W,kD,A0,F0,Jopt);

%%% material images and their signatures

for i=1:size(Aj,2)
    fig1 = figure(1);
    imagesc(reshape(Aj(:,i),n,n),[0 1]);
    axis square;colormap gray;set(gca,'xtick',[],'ytick',[]);
    saveas(fig1,[resDir 'Aj-' num2str(i)],'epsc');
    saveas(fig1,[resDir 'Aj-' num2str(i)],'fig');
end

fig2 = figure(2);
plot(Fj','LineWidth',2);pbaspect([3 1 1]);set(gca,'FontSize',16);
saveas(fig2,[resDir 'Fj'],'epsc');
saveas(fig2,[resDir 'Fj'],'fig');

%Save matrices
save([resDir 'Aj.mat'],'Aj');
save([resDir 'Fj.mat'],'Fj');

%% dictionary - ADJUST

fprintf('******* ADJUST \n')

kD = k;
A0 = rand(size(W,2),kD);
for i=1:size(A0,1), A0(i,:) = A0(i,:)/sum(A0(i,:)); end
R0 = rand(kD,size(T,1));

Dopt.rho     = 1e-2;
Dopt.iterMax = 1000;
Dopt.simFlag = 1;
Dopt.EQ_A    = 0;
Dopt.lambda  = 0;

[Ad,Fd,Rd,hist] = ADJUSTrev(Y,W,kD,T,A0,R0,Dopt);

%% material images and their signatures
for i=1:size(Ad,2)
    fig1 = figure(1);
    imagesc(reshape(Ad(:,i),n,n),[0 1]);
    axis square;colormap gray;set(gca,'xtick',[],'ytick',[]);
    saveas(fig1,[resDir 'Ad-' num2str(i)],'epsc');
    saveas(fig1,[resDir 'Ad-' num2str(i)],'fig');
end
%%
fig2 = figure(2);
plot(Fd','LineWidth',2);pbaspect([3 1 1]);set(gca,'FontSize',16);
saveas(fig2,[resDir 'Fd'],'epsc');
saveas(fig2,[resDir 'Fd'],'fig');

%Save matrices
save([resDir 'Ad.mat'],'Ad');
save([resDir 'Fd.mat'],'Fd');

%%

A1 = [A(:,1) A(:,2) A(:,3)+A(:,4)+A(:,5)];
CMAP  = colorcube(size(A1,2)+1);

[Atf,Ftf]           = computeFinalResultsA(A1,F,[resDir '/Atf'],[resDir '/Ftf'],CMAP,n);
[Adf,Fdf,eD,eMD]    = computeFinalResults(A1,Ad,Fd,[resDir '/Adf'],[resDir '/Fdf'],CMAP,n);
[Ajf,Fjf,eJ,eMJ]    = computeFinalResults(A1,Aj,Fj,[resDir '/Ajf'],[resDir '/Fjf'],CMAP,n);
[Aurf,Furf,eUR,eMUR]= computeFinalResults(A1,Aur,Fur,[resDir '/Aurf'],[resDir '/Furf'],CMAP,n);
[Aruf,Fruf,eRU,eMRU]= computeFinalResults(A1,Aru,Fru,[resDir '/Aruf'],[resDir '/Fruf'],CMAP,n);

fprintf('RU    :%.4f\nUR    :%.4f\nJoint :%.4f\nADJUST:%.4f\n',eRU,eUR,eJ,eD);

%%
str.True.A   = A;
str.True.A1  = A1;
str.True.F   = F;
str.True.Af  = Atf;
str.Joint.A  = Aj;
str.Joint.F  = Fj;
str.Joint.Af = Ajf;
str.Joint.Ff = Fjf;
str.Joint.err= eJ;
str.Joint.erM= eMJ;
str.ADJUST.A  = Ad;
str.ADJUST.F  = Fd;
str.ADJUST.Af = Adf;
str.ADJUST.Ff = Fdf;
str.ADJUST.err= eD;
str.ADJSUT.erM= eMD;
str.UR.A  = Aur;
str.UR.F  = Fur;
str.UR.Af = Aurf;
str.UR.Ff = Furf;
str.UR.err= eUR;
str.UR.erM= eMUR;
str.RU.A  = Aru;
str.RU.F  = Fru;
str.RU.Af = Aruf;
str.RU.Ff = Fruf;
str.RU.err= eRU;
str.RU.erM= eMRU;
str.n     = n;
str.nmat  = size(A,2);

save([resDir 'Full_results.mat'],'str');





