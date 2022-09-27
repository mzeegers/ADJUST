% Load all results
load('/home/mory/data/MatlabOneStep/Cai/Cai2013_simulations_fessler_100000_100000_30.mat', 'Cai2013_iterates', 'Cai2013_costs');
load('/home/mory/data/MatlabOneStep/Long/Long2014_simulations_300000_300000_30.mat', 'Long2014_iterates', 'Long2014_costs');
load('/home/mory/data/MatlabOneStep/Weidinger/Weidinger2016_simulations_30000_30000_3.mat', 'Weidinger2016_iterates', 'Weidinger2016_costs');
load('/home/mory/data/MatlabOneStep/Barber/Barber2016_simulations_100_100_10000.mat', 'Barber2016_iterates', 'Barber2016_costs');
load('/home/mory/data/MatlabOneStep/Mechlem/Mechlem2017_simulations_30000_30000_10.mat', 'Mechlem2017_iterates', 'Mechlem2017_costs'); %Something like that

% Multiply by the material's density
Cai2013_iterates(:,1,:) = Cai2013_iterates(:,1,:) * 4.933;              Cai2013_iterates(:,2,:) = Cai2013_iterates(:,2,:) * 7.9; 
Long2014_iterates(:,1,:) = Long2014_iterates(:,1,:) * 4.933;            Long2014_iterates(:,2,:) = Long2014_iterates(:,2,:) * 7.9; 
Weidinger2016_iterates(:,1,:) = Weidinger2016_iterates(:,1,:) * 4.933;  Weidinger2016_iterates(:,2,:) = Weidinger2016_iterates(:,2,:) * 7.9; 
Barber2016_iterates(:,1,:) = Barber2016_iterates(:,1,:) * 4.933;        Barber2016_iterates(:,2,:) = Barber2016_iterates(:,2,:) * 7.9; 
Mechlem2017_iterates(:,1,:) = Mechlem2017_iterates(:,1,:) * 4.933;      Mechlem2017_iterates(:,2,:) = Mechlem2017_iterates(:,2,:) * 7.9; 

% Write the final iterates for each method
results=cell(5,2);
results{1,1}=Cai2013_iterates;          results{1,2}='Cai2013';
results{2,1}=Long2014_iterates;         results{2,2}='Long2014';
results{3,1}=Weidinger2016_iterates;    results{3,2}='Weidinger2016';
results{4,1}=Barber2016_iterates;       results{4,2}='Barber2016';
results{5,1}=Mechlem2017_iterates;      results{5,2}='Mechlem2017';

% Generate ground truth
ObjectSize = sqrt(size(results{1,1}, 1));
unit = ObjectSize/8;
water=zeros(ObjectSize);
water(unit+1:7*unit, unit+1:7*unit) = ones(6*unit, 6*unit);
iodine=zeros(ObjectSize);
iodine(2*unit+1:3*unit, 2*unit+1:3*unit) = ones(unit, unit) * 0.01;
gadolinium=zeros(ObjectSize);
gadolinium(4*unit+1:5*unit, 5*unit+1:6*unit) = ones(unit, unit) * 0.01;

% Initialize the Relative MSE table
RMSE_table = zeros(5,3);

for m=1:5
    result = results{m,1};
    result = result(:,:,end);
    
    % Compute MSE with the ground truth
    RMSE_table(m, 1) = sum((result(:,1) - iodine(:)).^2)     / sum(iodine(:).^2);
    RMSE_table(m, 2) = sum((result(:,2) - gadolinium(:)).^2) / sum(gadolinium(:).^2);
    RMSE_table(m, 3) = sum((result(:,3) - water(:)).^2)      / sum(water(:).^2);
    
    % Write one image per material
    result = reshape(result, [ObjectSize ObjectSize 3]);
    writeImageResult(result(:,:,1), ['/home/mory/Documents/ArticlesAndNotes/OneStepComparison/figures/iodine_',     results{m,2}, '.png'], 0, 0.015);
    writeImageResult(result(:,:,2), ['/home/mory/Documents/ArticlesAndNotes/OneStepComparison/figures/gadolinium_', results{m,2}, '.png'], 0, 0.015);
    writeImageResult(result(:,:,3), ['/home/mory/Documents/ArticlesAndNotes/OneStepComparison/figures/water_',      results{m,2}, '.png'], 0, 1.5);
end

% To determine the concentration in ROIs, crop by a few pixels
% in order to avoid partial volume effect on the borders
cropROI=2;

%% Get iodine concentration curve for each, and plot it
figure(1);
hold off

[y_Cai, x_Cai, mean_Cai, std_Cai] = getMaterialConcentration('iodine', Cai2013_iterates, 10, cropROI);
semilogx(x_Cai, y_Cai, '+','linewidth', 2)
hold on

[y_Long, x_Long, mean_Long, std_Long] = getMaterialConcentration('iodine', Long2014_iterates, 10, cropROI);
semilogx(x_Long, y_Long, '--', 'linewidth', 2)

[y_Weidinger, x_Weidinger, mean_Weidinger, std_Weidinger] = getMaterialConcentration('iodine', Weidinger2016_iterates, 10, cropROI);
semilogx(x_Weidinger, y_Weidinger, '-.', 'linewidth', 2)

[y_Barber, x_Barber, mean_Barber, std_Barber] = getMaterialConcentration('iodine', Barber2016_iterates, 10, cropROI);
semilogx(x_Barber, y_Barber, 'o', 'linewidth', 2)

[y_Mechlem, x_Mechlem, mean_Mechlem, std_Mechlem] = getMaterialConcentration('iodine', Mechlem2017_iterates, 1, cropROI);
semilogx(x_Mechlem, y_Mechlem, 'x', 'linewidth', 2)

target_x = 1:max([x_Cai(end), x_Long(end), x_Weidinger(end), x_Barber(end), x_Mechlem(end)]);
target_y = ones(size(target_x))*0.01;
semilogx(target_x, target_y,'linewidth', 2)

lgd = legend('Cai2013','Long2014','Weidinger2016','Barber2016','Mechlem2017', 'Ground truth', 'Location', 'southeast');
lgd.FontSize = 14;
set(gca,'linewidth',2)
xlabel('Iterations','fontsize',14)
ylabel('Iodine concentration in ROI','fontsize',14)
set(gcf, 'Position', [500, 500, 600, 450])
hold off

% Rescale the graph manually before running this line, otherwise the legend
% appears on top of some curves
saveas(gcf,'/home/mory/Documents/ArticlesAndNotes/OneStepComparison/figures/iodineConcentrations.png');

% Store the mean and std values in a table
MeanStdTable(:,1:2)=[mean_Cai, std_Cai; mean_Long, std_Long; mean_Weidinger, std_Weidinger; mean_Barber, std_Barber; mean_Mechlem, std_Mechlem ];

%% Get gadolinium concentration curve for each, and plot it
figure(2);
hold off

[y_Cai, x_Cai ,mean_Cai, std_Cai] = getMaterialConcentration('gadolinium', Cai2013_iterates, 10, cropROI);
semilogx(x_Cai, y_Cai, '+','linewidth', 2)
hold on

[y_Long,x_Long, mean_Long, std_Long] = getMaterialConcentration('gadolinium', Long2014_iterates, 10, cropROI);
semilogx(x_Long, y_Long, '--', 'linewidth', 2)

[y_Weidinger, x_Weidinger, mean_Weidinger, std_Weidinger] = getMaterialConcentration('gadolinium', Weidinger2016_iterates, 10, cropROI);
semilogx(x_Weidinger, y_Weidinger, '-.', 'linewidth', 2)

[y_Barber,x_Barber, mean_Barber, std_Barber] = getMaterialConcentration('gadolinium', Barber2016_iterates, 10, cropROI);
semilogx(x_Barber, y_Barber, 'o', 'linewidth', 2)

[y_Mechlem, x_Mechlem, mean_Mechlem, std_Mechlem] = getMaterialConcentration('gadolinium', Mechlem2017_iterates, 1, cropROI);
semilogx(x_Mechlem, y_Mechlem, 'x', 'linewidth', 2)

target_x = 1:max([x_Cai(end), x_Long(end), x_Weidinger(end), x_Barber(end), x_Mechlem(end)]);
target_y = ones(size(target_x))*0.01;
semilogx(target_x, target_y,'linewidth', 2)

lgd = legend('Cai2013','Long2014','Weidinger2016','Barber2016','Mechlem2017', 'Ground truth', 'Location', 'southeast');
lgd.FontSize = 14;
set(gca,'linewidth',2)
xlabel('Iterations','fontsize',14)
ylabel('Gadolinium concentration in ROI','fontsize',14)
set(gcf, 'Position', [500, 500, 600, 450])
hold off

% Rescale the graph manually before running this line, otherwise the legend
% appears on top of some curves
saveas(gcf,'/home/mory/Documents/ArticlesAndNotes/OneStepComparison/figures/gadoliniumConcentrations.png');

% Store the mean and std values in a table
MeanStdTable(:,3:4)=[mean_Cai, std_Cai; mean_Long, std_Long; mean_Weidinger, std_Weidinger; mean_Barber, std_Barber; mean_Mechlem, std_Mechlem ];
MeanStdTable(:,1:4) = MeanStdTable(:,1:4) * 1000; % Contrast in mg/ml, water in g/ml

%% Get water concentrations, but don't plot it (only for mean & std)

[y_Cai, x_Cai, mean_Cai, std_Cai] = getMaterialConcentration('water', Cai2013_iterates, 10, cropROI);
[y_Long, x_Long, mean_Long, std_Long] = getMaterialConcentration('water', Long2014_iterates, 10, cropROI);
[y_Weidinger, x_Weidinger, mean_Weidinger, std_Weidinger] = getMaterialConcentration('water', Weidinger2016_iterates, 10, cropROI);
[y_Barber, x_Barber, mean_Barber, std_Barber] = getMaterialConcentration('water', Barber2016_iterates, 10, cropROI);
[y_Mechlem, x_Mechlem, mean_Mechlem, std_Mechlem] = getMaterialConcentration('water', Mechlem2017_iterates, 1, cropROI);

% Store the mean and std values in a table
MeanStdTable(:,5:6)=[mean_Cai, std_Cai; mean_Long, std_Long; mean_Weidinger, std_Weidinger; mean_Barber, std_Barber; mean_Mechlem, std_Mechlem ];

%% Plot the costs functions
Cai2013_costs_toplot = Cai2013_costs(10:10:end);
Cai2013_costs_toplot = Cai2013_costs_toplot - min(Cai2013_costs_toplot);

Long2014_costs_toplot = Long2014_costs(10:10:end);
Long2014_costs_toplot = Long2014_costs_toplot - min(Long2014_costs_toplot);

Weidinger2016_costs_toplot = Weidinger2016_costs(10:10:end);
Weidinger2016_costs_toplot = Weidinger2016_costs_toplot - min(Weidinger2016_costs_toplot);

Barber2016_costs_toplot = Barber2016_costs;
Barber2016_costs_toplot = Barber2016_costs_toplot - min(Barber2016_costs_toplot);

Mechlem2017_costs_toplot = Mechlem2017_costs;
Mechlem2017_costs_toplot = Mechlem2017_costs_toplot - min(Mechlem2017_costs_toplot);

figure(3);
hold off
loglog(x_Cai,Cai2013_costs_toplot, '+','linewidth', 2)
hold on 
loglog(x_Long, Long2014_costs_toplot, '--','linewidth', 2)
loglog(x_Weidinger, Weidinger2016_costs_toplot, '-.','linewidth', 2)
load('/home/mory/data/MatlabOneStep/Barber/preps.mat', 'y', 'A', 'M', 'S', 'T');
% Barber2016_costs = zeros(size(Barber2016_iterates, 3), 1);
% for it=1:size(Barber2016_iterates, 3)
%     Barber2016_costs(it) = BarberComputeCosts(Barber2016_iterates(:,:,it), y, A, M, S);
% end
loglog(1:10000, Barber2016_costs_toplot, 'o','linewidth', 2)
loglog(x_Mechlem, Mechlem2017_costs_toplot, 'x','linewidth', 2)

lgd2 = legend('Cai2013','Long2014','Weidinger2016','Barber2016','Mechlem2017', 'Location', 'southwest');
lgd2.FontSize = 14;
set(gca,'linewidth',2)
xlabel('Iterations','fontsize',14)
ylabel('Costs - min(Costs)','fontsize',14)
set(gcf, 'Position', [500, 500, 600, 450])
hold off

% Rescale the graph manually before running this line, otherwise the legend
% appears on top of some curves
saveas(gcf,'/home/mory/Documents/ArticlesAndNotes/OneStepComparison/figures/CostFunctions.png');

%% Thresholds in the concentration of iodine and gadolinium
Cai_10_percent = GetFirstIterationMeetingTolerance(Cai2013_iterates, 10, 0.1, cropROI);
Long_10_percent= GetFirstIterationMeetingTolerance(Long2014_iterates, 10, 0.1, cropROI);
Weidinger_10_percent= GetFirstIterationMeetingTolerance(Weidinger2016_iterates, 10, 0.1, cropROI);
Barber_10_percent= GetFirstIterationMeetingTolerance(Barber2016_iterates, 10, 0.1, cropROI);
Mechlem_10_percent= GetFirstIterationMeetingTolerance(Mechlem2017_iterates, 1, 0.1, cropROI);

Cai_20_percent= GetFirstIterationMeetingTolerance(Cai2013_iterates, 10, 0.2, cropROI);
Long_20_percent= GetFirstIterationMeetingTolerance(Long2014_iterates, 10, 0.2, cropROI);
Weidinger_20_percent= GetFirstIterationMeetingTolerance(Weidinger2016_iterates, 10, 0.2, cropROI);
Barber_20_percent= GetFirstIterationMeetingTolerance(Barber2016_iterates, 10, 0.2, cropROI);
Mechlem_20_percent= GetFirstIterationMeetingTolerance(Mechlem2017_iterates, 1, 0.2, cropROI);

thresholdsTable = [Cai_20_percent,       Cai_10_percent; ...
                   Long_20_percent,      Long_10_percent; ... 
                   Weidinger_20_percent, Weidinger_10_percent; ...
                   Barber_20_percent,    Barber_10_percent; ...
                   Mechlem_20_percent,   Mechlem_10_percent];

timePerIter=[5.2; 265/200; 95/200; 204/200; 124/200];

% Where the threshold hasn't been reached, put value 10000 instead of 0
thresholdsTable(thresholdsTable==0)=10000;

% Get the times to reach the thresholds
timeTable=thresholdsTable .* timePerIter;

%% Compare Cai2013 with different mu-preconditioning

Cai2013_none = load('/home/mory/data/MatlabOneStep/Cai/Cai2013_simulations_none_100000_100000_30.mat');
Cai2013_none_iterates = Cai2013_none.Cai2013_iterates;
Cai2013_none_costs = Cai2013_none.Cai2013_costs;
Cai2013_normalize = load('/home/mory/data/MatlabOneStep/Cai/Cai2013_simulations_normalize_100000_100000_30.mat');
Cai2013_normalize_iterates = Cai2013_normalize.Cai2013_iterates;
Cai2013_normalize_costs = Cai2013_normalize.Cai2013_costs;
Cai2013_orthonormalize = load('/home/mory/data/MatlabOneStep/Cai/Cai2013_simulations_orthonormalize_100000_100000_30.mat');
Cai2013_orthonormalize_iterates = Cai2013_orthonormalize.Cai2013_iterates;
Cai2013_orthonormalize_costs = Cai2013_orthonormalize.Cai2013_costs;
Cai2013_fessler = load('/home/mory/data/MatlabOneStep/Cai/Cai2013_simulations_fessler_100000_100000_30.mat');
Cai2013_fessler_iterates = Cai2013_fessler.Cai2013_iterates;
Cai2013_fessler_costs = Cai2013_fessler.Cai2013_costs;

Mini = min(cat(1, Cai2013_none_costs(:), Cai2013_normalize_costs(:), Cai2013_orthonormalize_costs(:), Cai2013_fessler_costs(:)));

Cai2013_none_costs_toplot = Cai2013_none_costs(10:10:end);
Cai2013_none_costs_toplot = Cai2013_none_costs_toplot - Mini;

Cai2013_normalize_costs_toplot = Cai2013_normalize_costs(10:10:end);
Cai2013_normalize_costs_toplot = Cai2013_normalize_costs_toplot - Mini;

Cai2013_orthonormalize_costs_toplot = Cai2013_orthonormalize_costs(10:10:end);
Cai2013_orthonormalize_costs_toplot = Cai2013_orthonormalize_costs_toplot - Mini;

Cai2013_fessler_costs_toplot = Cai2013_fessler_costs(10:10:end);
Cai2013_fessler_costs_toplot = Cai2013_fessler_costs_toplot - Mini;

x_Cai = (1:size(Cai2013_none_iterates, 3)) * 10;

figure(4);
hold off
loglog(x_Cai,Cai2013_none_costs_toplot, '+','linewidth', 2)
hold on 
loglog(x_Cai, Cai2013_normalize_costs_toplot, '--','linewidth', 2)
loglog(x_Cai, Cai2013_orthonormalize_costs_toplot, '-.','linewidth', 2)
loglog(x_Cai, Cai2013_fessler_costs_toplot, 'x','linewidth', 2)

lgd3 = legend('none','normalize','orthonormalize','Fessler', 'Location', 'southwest');
lgd3.FontSize = 14;
set(gca,'linewidth',2)
xlabel('Iterations','fontsize',14)
ylabel('Costs - min(Costs)','fontsize',14)
% set(gcf, 'Position', [500, 500, 600, 450])
hold off

saveas(gcf,'/home/mory/Documents/ArticlesAndNotes/OneStepComparison/figures/mupreconditionings_costs.png');

%% Write a video with the simulation results

% First set all results to the same length, ie 10000 iterations, and 1000
% iterates saved (one every ten)
niterates = max([size(Cai2013_iterates, 3), size(Long2014_iterates, 3), size(Weidinger2016_iterates, 3), size(Barber2016_iterates, 3), size(Mechlem2017_iterates, 3)/10]);
nmat = 3;
videoCai2013_iterates = repmat(Cai2013_iterates(:,:,end), [1 1 niterates]);
videoCai2013_iterates(:,:,1:size(Cai2013_iterates, 3)) = Cai2013_iterates;
videoCai2013_iterates = permute(reshape(videoCai2013_iterates, [ObjectSize, ObjectSize, nmat, niterates]), [2 1 3 4]);

videoLong2014_iterates = repmat(Long2014_iterates(:,:,end), [1 1 niterates]);
videoLong2014_iterates(:,:,1:size(Long2014_iterates, 3)) = Long2014_iterates;
videoLong2014_iterates = permute(reshape(videoLong2014_iterates, [ObjectSize, ObjectSize, nmat, niterates]), [2 1 3 4]);

videoWeidinger2016_iterates = repmat(Weidinger2016_iterates(:,:,end), [1 1 niterates]);
videoWeidinger2016_iterates(:,:,1:size(Weidinger2016_iterates, 3)) = Weidinger2016_iterates;
videoWeidinger2016_iterates = permute(reshape(videoWeidinger2016_iterates, [ObjectSize, ObjectSize, nmat, niterates]), [2 1 3 4]);

videoBarber2016_iterates = repmat(Barber2016_iterates(:,:,end), [1 1 niterates]);
videoBarber2016_iterates(:,:,1:size(Barber2016_iterates, 3)) = Barber2016_iterates;
videoBarber2016_iterates = permute(reshape(videoBarber2016_iterates, [ObjectSize, ObjectSize, nmat, niterates]), [2 1 3 4]);

videoMechlem2017_iterates = repmat(Mechlem2017_iterates(:,:,end), [1 1 niterates]);
videoMechlem2017_iterates(:,:,1:size(Mechlem2017_iterates, 3)/10) = Mechlem2017_iterates(:,:,10:10:end);
videoMechlem2017_iterates = permute(reshape(videoMechlem2017_iterates, [ObjectSize, ObjectSize, nmat, niterates]), [2 1 3 4]);

% Concatenate
video_iterates = cat(1, videoCai2013_iterates, videoLong2014_iterates, videoWeidinger2016_iterates, videoBarber2016_iterates, videoMechlem2017_iterates);

% Scale to unity (apply window level)
video_iterates(:,:,1,:) = video_iterates(:,:,1,:) / 0.015;
video_iterates(:,:,2,:) = video_iterates(:,:,2,:) / 0.015;
video_iterates(:,:,3,:) = video_iterates(:,:,3,:) / 1.5;

%  reshape
video_iterates = reshape(video_iterates, [ObjectSize * 5, ObjectSize * nmat, niterates]);

% Write as video
% Video is too big to be written by Matlab at once. Split it into chunks of
% 100 iterates
% Then merge them using https://www.aconvert.com/video/merge/#
for k=0:(niterates/100)-1
    writeVideoResult(video_iterates(:,:,100*k+1:100*k+100), ['/home/mory/Documents/ArticlesAndNotes/OneStepComparison/supplementary/simulations_iterates_', num2str(k), '.avi'], 0, 1);
end

%% Write final iterates for the real data results

% Load all results
load('/home/mory/data/MatlabOneStep/Cai/Cai2013_rabbit_fessler_10000_10000_10.mat', 'Cai2013_iterates');
load('/home/mory/data/MatlabOneStep/Long/Long2014_rabbit_10000_10000_10.mat', 'Long2014_iterates');
load('/home/mory/data/MatlabOneStep/Weidinger/Weidinger2016_rabbit_3000_3000_1.mat', 'Weidinger2016_iterates');
load('/home/mory/data/MatlabOneStep/Barber/Barber2016_rabbit_200_200_20000.mat', 'Barber2016_iterates');
load('/home/mory/data/MatlabOneStep/Mechlem/Mechlem2017_rabbit_10000_10000_3.mat', 'Mechlem2017_iterates');

% Multiply by the material's density
Cai2013_iterates(:,1,:) = Cai2013_iterates(:,1,:) * 4.933;              Cai2013_iterates(:,2,:) = Cai2013_iterates(:,2,:) * 7.9; 
Long2014_iterates(:,1,:) = Long2014_iterates(:,1,:) * 4.933;            Long2014_iterates(:,2,:) = Long2014_iterates(:,2,:) * 7.9; 
Weidinger2016_iterates(:,1,:) = Weidinger2016_iterates(:,1,:) * 4.933;  Weidinger2016_iterates(:,2,:) = Weidinger2016_iterates(:,2,:) * 7.9; 
Barber2016_iterates(:,1,:) = Barber2016_iterates(:,1,:) * 4.933;        Barber2016_iterates(:,2,:) = Barber2016_iterates(:,2,:) * 7.9; 
Mechlem2017_iterates(:,1,:) = Mechlem2017_iterates(:,1,:) * 4.933;      Mechlem2017_iterates(:,2,:) = Mechlem2017_iterates(:,2,:) * 7.9; 

% Write the final iterates for each method
results=cell(5,2);
results{1,1}=Cai2013_iterates;          results{1,2}='Cai2013';
results{2,1}=Long2014_iterates;         results{2,2}='Long2014';
results{3,1}=Weidinger2016_iterates;    results{3,2}='Weidinger2016';
results{4,1}=Barber2016_iterates;       results{4,2}='Barber2016';
results{5,1}=Mechlem2017_iterates;      results{5,2}='Mechlem2017';

% Read the results from the scanner
water=read_mhd('/home/mory/data/Hamburg/data_2018_Feb_13.9714.1/trans_water.mhd');              water_slice=squeeze(water.data(:,2,:));
iodine=read_mhd('/home/mory/data/Hamburg/data_2018_Feb_13.9714.1/trans_iodine.mhd');            iodine_slice=squeeze(iodine.data(:,2,:));
gadolinium=read_mhd('/home/mory/data/Hamburg/data_2018_Feb_13.9714.1/trans_gadolinium.mhd');    gadolinium_slice=squeeze(gadolinium.data(:,2,:));
fovmask=read_mhd('/home/mory/data/Hamburg/data_2018_Feb_13.9714.1/trans_fovmask.mhd');          fov_slice=squeeze(fovmask.data(:,2,:));

% Write the scanner results as png with the same grayscale as the one step
% results
writeImageResult(fliplr(iodine_slice.'), '/home/mory/Documents/ArticlesAndNotes/OneStepComparison/figures/Rabbit_iodine_scanner.png', 0, 20);
writeImageResult(fliplr(gadolinium_slice.'), '/home/mory/Documents/ArticlesAndNotes/OneStepComparison/figures/Rabbit_gadolinium_scanner.png', 0, 7);
writeImageResult(fliplr(water_slice.'), '/home/mory/Documents/ArticlesAndNotes/OneStepComparison/figures/Rabbit_water_scanner.png', 0, 1500);

for m=1:5
    result = results{m,1};
    result = result(:,:,end);
    ObjectSize = sqrt(numel(result)/3);
    
    % Write one image per material
    result = reshape(result, [ObjectSize ObjectSize 3]);
    result = result .* fov_slice;
    writeImageResult(fliplr(result(:,:,1).'), ['/home/mory/Documents/ArticlesAndNotes/OneStepComparison/figures/Rabbit_iodine_',     results{m,2}, '.png'], 0, 0.02);
    writeImageResult(fliplr(result(:,:,2).'), ['/home/mory/Documents/ArticlesAndNotes/OneStepComparison/figures/Rabbit_gadolinium_', results{m,2}, '.png'], 0, 0.007);
    writeImageResult(fliplr(result(:,:,3).'), ['/home/mory/Documents/ArticlesAndNotes/OneStepComparison/figures/Rabbit_water_',      results{m,2}, '.png'], 0, 1.5);
end

%% Write a video with the real data results

% First set all results to the same length, ie 10000 iterations, and 1000
% iterates saved (one every ten)
niterates = max([size(Cai2013_iterates, 3), size(Long2014_iterates, 3), size(Weidinger2016_iterates, 3), size(Barber2016_iterates, 3), size(Mechlem2017_iterates, 3)/10]);
nmat = 3;
videoCai2013_iterates = repmat(Cai2013_iterates(:,:,end), [1 1 niterates]);
videoCai2013_iterates(:,:,1:size(Cai2013_iterates, 3)) = Cai2013_iterates;
videoCai2013_iterates = flip(permute(reshape(videoCai2013_iterates, [ObjectSize, ObjectSize, nmat, niterates]), [2 1 3 4]), 2);
videoCai2013_iterates = videoCai2013_iterates .* fov_slice;

videoLong2014_iterates = repmat(Long2014_iterates(:,:,end), [1 1 niterates]);
videoLong2014_iterates(:,:,1:size(Long2014_iterates, 3)) = Long2014_iterates;
videoLong2014_iterates = flip(permute(reshape(videoLong2014_iterates, [ObjectSize, ObjectSize, nmat, niterates]), [2 1 3 4]), 2);
videoLong2014_iterates = videoLong2014_iterates .* fov_slice;

videoWeidinger2016_iterates = repmat(Weidinger2016_iterates(:,:,end), [1 1 niterates]);
videoWeidinger2016_iterates(:,:,1:size(Weidinger2016_iterates, 3)) = Weidinger2016_iterates;
videoWeidinger2016_iterates = flip(permute(reshape(videoWeidinger2016_iterates, [ObjectSize, ObjectSize, nmat, niterates]), [2 1 3 4]), 2);
videoWeidinger2016_iterates = videoWeidinger2016_iterates .* fov_slice;

videoBarber2016_iterates = repmat(Barber2016_iterates(:,:,end), [1 1 niterates]);
videoBarber2016_iterates(:,:,1:size(Barber2016_iterates, 3)) = Barber2016_iterates;
videoBarber2016_iterates = flip(permute(reshape(videoBarber2016_iterates, [ObjectSize, ObjectSize, nmat, niterates]), [2 1 3 4]), 2);
videoBarber2016_iterates = videoBarber2016_iterates .* fov_slice;

videoMechlem2017_iterates = repmat(Mechlem2017_iterates(:,:,end), [1 1 niterates]);
videoMechlem2017_iterates(:,:,1:size(Mechlem2017_iterates, 3)/10) = Mechlem2017_iterates(:,:,10:10:end);
videoMechlem2017_iterates = flip(permute(reshape(videoMechlem2017_iterates, [ObjectSize, ObjectSize, nmat, niterates]), [2 1 3 4]), 2);
videoMechlem2017_iterates = videoMechlem2017_iterates .* fov_slice;

% Concatenate
video_iterates = cat(1, videoCai2013_iterates, videoLong2014_iterates, videoWeidinger2016_iterates, videoBarber2016_iterates, videoMechlem2017_iterates);

% Scale to unity (apply window level)
video_iterates(:,:,1,:) = video_iterates(:,:,1,:) / 0.02;
video_iterates(:,:,2,:) = video_iterates(:,:,2,:) / 0.007;
video_iterates(:,:,3,:) = video_iterates(:,:,3,:) / 1.5;

%  reshape
video_iterates = reshape(video_iterates, [ObjectSize * 5, ObjectSize * nmat, niterates]);

% Write as video
% Video is too big to be written by Matlab at once. Split it into chunks of
% 100 iterates
% Then merge them using https://www.aconvert.com/video/merge/#
for k=0:(niterates/100)-1
    writeVideoResult(video_iterates(:,:,100*k+1:100*k+100), ['/home/mory/Documents/ArticlesAndNotes/OneStepComparison/supplementary/rabbit_iterates_', num2str(k), '.avi'], 0, 1);
end
