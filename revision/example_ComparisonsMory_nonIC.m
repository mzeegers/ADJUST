%%% X-ray test script (with dataset for F and T)
%%% Dictionary approach
%%% Experimental data

clearvars;

%Long or short run
LongRun = true;
%Apply scaling by density (and some factor of 0.01) as in Mory paper
Scaling = true;

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
% if(LongRun)
%     resDir = [rwd '/results/ex15E_Comparison_LongRun_' startTime '/'];
% else
%     resDir = [rwd '/results/ex15E_Comparison_' startTime '/'];
% end
% if ~exist(resDir,'dir'), mkdir(resDir); end

% setting random stream
myStream = RandStream('mt19937ar','Seed',10);
RandStream.setGlobalStream(myStream)
 
%% generate U

SpectraFromMoryPaper = true;

% Custom attenuation spectra
Dict = readmatrix('ComparisonPaper_MatrixTXray100chan_93_MinEn20_MaxEn119_mat_CORRECTED.csv');
M_Dict = cat(1, Dict(53,:), Dict(64,:), Dict(85,:));

% Spectra from the paper
load('incidentSpectrum.mat', 'incidentSpectrum');
load('materialAttenuations.mat', 'materialAttenuations');
load('detectorResponse.mat', 'detectorResponse');

if(LongRun)
    ObjectSize = 128;
else
    ObjectSize = 64;
end
NbPixelsPerProj = round(ObjectSize * sqrt(2));

% Get them under matrix form, consistent with Cai's paper
materialAttenuations = materialAttenuations(20:119,:);

if(SpectraFromMoryPaper)
    M = materialAttenuations;
else
    M = M_Dict';
end    

S_unbinned = incidentSpectrum.';
S = S_unbinned(:);
S = S(20:119);
S = diag(S);


% Generate a simple phantom
unit = ObjectSize/8;
water=zeros(ObjectSize);
water(unit+1:7*unit, unit+1:7*unit) = ones(6*unit, 6*unit);
iodine=zeros(ObjectSize);
iodine(2*unit+1:3*unit, 2*unit+1:3*unit) = ones(unit, unit); %  * 0.01 / 4.933;
gadolinium=zeros(ObjectSize);
gadolinium(4*unit+1:5*unit, 5*unit+1:6*unit) = ones(unit, unit); %  * 0.01 / 7.9;

% Generate a larger version of the phantom
unit2 = (ObjectSize*2)/8;
water2=zeros(ObjectSize*2);
water2(unit2+1:7*unit2, unit2+1:7*unit2) = ones(6*unit2, 6*unit2);
iodine2=zeros(ObjectSize*2);
iodine2(2*unit2+1:3*unit2, 2*unit2+1:3*unit2) = ones(unit2, unit2); %  * 0.01 / 4.933;
gadolinium2=zeros(ObjectSize*2);
gadolinium2(4*unit2+1:5*unit2, 5*unit2+1:6*unit2) = ones(unit2, unit2); %  * 0.01 / 7.9;


% Compute the number of projections required
% to keep a pixels-to-voxels ratio of 4
NbProj = ceil(ObjectSize * 4 / sqrt(2));
ForwardTheta = [0:NbProj-1] * 180 / NbProj;

n_ang = NbProj;
theta = linspace(0,pi,n_ang); 

% Generate forward projection matrix
volGeo = astra_create_vol_geom(ObjectSize*2, ObjectSize*2);
projGeo = astra_create_proj_geom('parallel', 2, NbPixelsPerProj, theta);
Wma = opTomo('linear', projGeo, volGeo);
Wm = paralleltomo(ObjectSize*2,ForwardTheta, round(sqrt(2)*ObjectSize*2), sqrt(2)*ObjectSize*2);

%default: p = round(sqrt(2)*N);d = sqrt(2)*N;

% (non-inverse-crime) forward operator to be used for inversion
volGeo = astra_create_vol_geom(ObjectSize, ObjectSize);
projGeo = astra_create_proj_geom('parallel', 1, NbPixelsPerProj, theta);
Wa = opTomo('linear', projGeo, volGeo);
W = paralleltomo(ObjectSize,ForwardTheta);

% Compute forward projections
object = cat(2, iodine(:), gadolinium(:), water(:));
object2 = cat(2, iodine2(:), gadolinium2(:), water2(:));

object(:,3) = object(:,3) - logical(object(:,1) + object(:,2));
if(Scaling)
    object(:,1) = object(:,1)* 0.01 / 4.933;
    object(:,2) = object(:,2)* 0.01 / 7.9;
end

object2(:,3) = object2(:,3) - logical(object2(:,1) + object2(:,2));
if(Scaling)
    object2(:,1) = object2(:,1)* 0.01 / 4.933;
    object2(:,2) = object2(:,2)* 0.01 / 7.9;
end


%%

% Compute forward projections
projsm = Wm * object2/2;
projs = W * object;

projsma = Wma * object2/2;
projsa = Wa * object;

projsdown = zeros(size(projs,1), size(projs,2));
for i=1:size(projs,2)
    projsdown(:,i) = reshape(imresize(reshape(projsm(:,i), NbPixelsPerProj*2, n_ang), [NbPixelsPerProj n_ang]), n_ang*NbPixelsPerProj, 1);
end

ym = M*projsdown';
y = M*projs';
%%
% load matrices
n = ObjectSize;
T = M';

% add an extra spectrum (constant) to the dictionary
%T(end+1,:) = (1e-3*max(T(:)))*ones(1,size(T,2));

% measurements with noise
Y1 = y';
Y    = ym';
for i=1:size(Y,2)
    Y(:,i) = astra_add_noise_to_sino(Y(:,i),1e5);
end
Y_methods = exp(-Y)*S;
%%
%For Cai's method
attenuationFactors = exp(-y);
counts_without_attenuation = S * ones(size(attenuationFactors));
Y_Cai = (S*exp(-Y'))./counts_without_attenuation;
Y_Cai = Y_Cai';
y_without_attenuation = poissrnd(counts_without_attenuation) ./ counts_without_attenuation;
k_d = var(y_without_attenuation(:)) / mean(y_without_attenuation(:));

noAttenuation = S * ones(size(M,1),1);
S_Cai = S ./ repmat(noAttenuation, [1, size(S, 2)]);

% condition T (remove dependent rows)
% T = transpose(licols(T',1e-9));

% sizes
m = ObjectSize*ObjectSize; % number of images
c = 100; % number of channels

k = 3; %size(A,2);     % number of materials


%fprintf('cond(Ftrue) = %.2f \n',cond(F));
%Su = svd(U0);
%fprintf('cond(U)     = %.2f \n',Su(1)/Su(k));
%fprintf('size(U): %d x %d | k = %d \n',size(U0),k);
fprintf('size(W): %d x %d \n',size(W));
fprintf('size(T): %d x %d \n',size(T));


%% Reconstruction and then Unmixing

fprintf('******* reconstruction-then-unmixing \n')
Wty = W'*Y;
WtW = @(x) W'*(W*x) + 1e-3*x;

Up = [];
for i=1:size(Y,2)
   [Up(:,i),~] = pcg(WtW,Wty(:,i),1e-6,2);
end

Up(Up <0) = 0;

opt = statset('MaxIter',100,'Display','final');
[A10,F10] = nnmf(Up,k,'Replicates',20,'Options',opt,'Algorithm','mult');

opt = statset('MaxIter',1000,'Display','final');
[Aru,Fru] = nnmf(Up,k,'W0',A10,'H0',F10,'Options',opt,'Algorithm','als','Replicates',20);


%% Unmixing then Reconstruction

fprintf('******* unmixing-then-reconstruction\n')
opt = statset('MaxIter',10,'Display','final');
[A20,F20] = nnmf(Y,k,'Replicates',10,'Options',opt,'Algorithm','mult','Replicates',10);

opt = statset('MaxIter',1000,'Display','final');
[WAur,Fur] = nnmf(Y,k,'W0',A20,'H0',F20,'Options',opt,'Algorithm','als','Replicates',10);


Wtp = W'*WAur;
WtW = @(x) W'*(W*x) + 1e-3*x;

for i=1:size(WAur,2)
   [Aur(:,i),~] = pcg(WtW,Wtp(:,i),1e-6,5);
end


%% joint - JUST

fprintf('******* unmixing-then-reconstruction\n')

kD = k;
A0 = rand(size(W,2),kD,'single');
F0 = rand(kD,size(Y,2),'single');

Jopt.rho     = 1e-2;
if(LongRun)
    Jopt.iterMax = 2000;
else
    Jopt.iterMax = 200;
end
Jopt.simFlag = 1;
Jopt.EQ_A    = 0;
Jopt.lambda  = 0;

[Aj,Fj,histJ] = cJoint(Y,W*eye(size(W,2)),kD,Jopt);


%% dictionary - ADJUST

fprintf('******* ADJUST \n')

kD = k;
A0 = rand(size(W,2),kD);
for i=1:size(A0,1), A0(i,:) = A0(i,:)/sum(A0(i,:)); end
R0 = rand(kD,size(T,1));

Dopt.rho     = 1e-2;
if(LongRun)
    Dopt.iterMax = 2000;
else
    Dopt.iterMax = 200;
end
Dopt.innerIter = 4;
Dopt.simFlag = 1;
Dopt.EQ_A    = 0;
Dopt.lambda  = 0;

[Ad,Fd,Rd,hist] = ADJUST(Y,W,kD,diag([0.01/4.933;0.01/7.9;1])*T,Dopt);

%% Cai (works)

%Fessler trick
inv_passage = (S * M )./(S * ones(size(M))); % the beta_bar matrix in equation 71 of Fessler's patent
passage = inv(inv_passage.' * inv_passage) * inv_passage.'; % The Moore-Penrose pseudo inverse, since the inv_passage matrix is not square
out = M * passage;
T_Cai = passage;
M_Cai = out;
S_Cai = repmat(S_Cai, [1 1 NbPixelsPerProj]);

if(LongRun)
    NbIters_Cai = 5000; %long run
else
    NbIters_Cai = 200; %short run
end
delta_huber = [0.001, 0.001, 0.1];
runOnGPU = true;
useNesterov = false;
regulWeights = [100000, 100000, 100];
[Cai2013_iterates, Cai2013_costs ]= Cai2013(Y_Cai, ObjectSize, W, M_Cai, S_Cai, k_d, NbIters_Cai, T_Cai, regulWeights, delta_huber, runOnGPU, useNesterov);


%% Weidinger (works, results are fine)

S_Wei = repmat(S, [1 1 NbPixelsPerProj]); % Same spectrum, repeated for all pixels. Just to test for pixel-dependent spectrum feature
T_Wei = eye(3);

% Perform reconstruction
if(LongRun)
    NbIters_Wei = 5000; %long run
else
    NbIters_Wei = 20; %short run
end
runOnGPU = true;
regulWeights = [10000, 10000, 10];
[Weidinger2016_iterates, Weidinger2016_costs ]= Weidinger2016(Y_methods, ObjectSize, W, M, S_Wei, T_Wei, NbIters_Wei, regulWeights, runOnGPU);


%% Long (works, results are fine)
S_Long = repmat(S, [1 1 NbPixelsPerProj]);
T_Long = eye(3);

if(LongRun)
    NbIters_Long = 5000; %long run
else
    NbIters_Long = 10; %short run
end
delta_huber = [0.001, 0.001, 0.1];
runOnGPU = true;
NbSubsets = 20; % For ordered subsets
regulWeights = [100000, 100000, 100];
[Long2014_iterates, Long2014_costs]= Long2014(Y_methods, ObjectSize, NbPixelsPerProj, W, M, S_Long, T_Long, NbIters_Long, regulWeights, NbSubsets, delta_huber, runOnGPU);


%% Mechlem (works, results are fine)
S_Mech = repmat(S, [1 1 NbPixelsPerProj]);
T_Mech = eye(3);

if(LongRun)
    NbIters_Mech = 200; %long run
else
    NbIters_Mech = 50; %short run
end
delta_huber = [0.001, 0.001, 0.1];
runOnGPU = true;
NbSubsets = 4; % More than that speeds up computation, but makes it unstable, and ultimately, divergent
regulWeights = [100000, 100000, 100]; % Huber is scaled by delta, so where there is a low delta, there should be a high lambda to compensate
[Mechlem2017_iterates, Mechlem2017_costs ]= Mechlem2017(Y_methods, ObjectSize, NbPixelsPerProj, W, M, S_Mech, T_Mech, NbIters_Mech, regulWeights, NbSubsets, delta_huber, runOnGPU);


%% Barber (works, results are fine)

%Normalize
normalized = M;
passage = eye(3);
temp_passage_mat = eye(3); temp_passage_mat(1,1)= 1/norm(normalized(:,1));
normalized = normalized * temp_passage_mat;
passage = passage * temp_passage_mat;

temp_passage_mat = eye(3); temp_passage_mat(2,2)= 1/norm(normalized(:,2));
normalized = normalized * temp_passage_mat;
passage = passage * temp_passage_mat;

temp_passage_mat = eye(3); temp_passage_mat(3,3)= 1/norm(normalized(:,3));
normalized = normalized * temp_passage_mat;
passage = passage * temp_passage_mat;

M_Barber = normalized;
T_Barber = passage;

if(LongRun)
    NbIters_Barber = 10000; %long run
else
    NbIters_Barber = 200; %short run
end

runOnGPU = false;
lambda = 0.0001; % Controls convergence
theta = 0.5; % Between 0 and 1, see Chambolle Pock
TVthresholds = [100,100,10000];
[Barber2016_iterates, Barber2016_gaps, Barber2016_costs]= Barber2016(Y_methods, ObjectSize, W, M_Barber, S, T_Barber, NbIters_Barber, lambda, theta, TVthresholds, runOnGPU);

%%
str.True.SpectraFromMoryPaper = SpectraFromMoryPaper;
str.True.LongRun = LongRun;
str.True.Scaling = Scaling;
str.True.A   = object;
str.True.F   = M;
str.True.T   = T;
str.UR.A  = Aur;
str.UR.F  = Fur;
str.RU.A  = Aru;
str.RU.F  = Fru;
str.Joint.A  = Aj;
str.Joint.F  = Fj;
str.ADJUST.A  = Ad;
str.ADJUST.F  = Fd;
str.Cai.A = Cai2013_iterates(:,:,end);
str.Weidinger.A = Weidinger2016_iterates(:,:,end);
str.Long.A = Long2014_iterates(:,:,end);
str.Mechlem.A = Mechlem2017_iterates(:,:,end);
str.Barber.A = Barber2016_iterates(:,:,end);
str.n     = n;

save('example_ComparisonsMory_nonIC.mat','str');

%%
strProb.A = str.True.A;
strProb.F = str.True.F;
strProb.k = k;
strProb.T = str.True.T;
strProb.n = [n n];

strRU.A   = Aru;
strRU.F   = Fru;
strUR.A   = Aur;
strUR.F   = Fur;
strJ.A = Aj;
strJ.F = Fj;
strD.A = Ad;
strD.F = Fd;

strCai.A = str.Cai.A;
strCai.F = str.True.F;
strWeidinger.A = str.Weidinger.A;
strWeidinger.F = str.True.F;
strLong.A = str.Long.A;
strLong.F = str.True.F;
strMechlem.A = str.Mechlem.A;
strMechlem.F = str.True.F;
strBarber.A = str.Barber.A;
strBarber.F = str.True.F;

strCai.A(:,1) = strCai.A(:,1)/(0.01/4.933)
strCai.A(:,2) = strCai.A(:,2)/(0.01/7.9) 
strWeidinger.A(:,1) = strWeidinger.A(:,1)/(0.01/4.933)
strWeidinger.A(:,2) = strWeidinger.A(:,2)/(0.01/7.9) 
strLong.A(:,1) = strLong.A(:,1)/(0.01/4.933)
strLong.A(:,2) = strLong.A(:,2)/(0.01/7.9) 
strMechlem.A(:,1) = strMechlem.A(:,1)/(0.01/4.933)
strMechlem.A(:,2) = strMechlem.A(:,2)/(0.01/7.9) 
strBarber.A(:,1) = strBarber.A(:,1)/(0.01/4.933)
strBarber.A(:,2) = strBarber.A(:,2)/(0.01/7.9) 




%%
strOut = computeResultsMory(strProb,strRU,strUR,strJ,strD,strCai,strWeidinger,strLong,strMechlem,strBarber);
