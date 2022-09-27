lena = imread('/home/mory/Pictures/lena_std.tif');
lena = sum(double(lena), 3);
noisyLena = lena + (rand(size(lena)) - 0.5) * 100;
figure(1); Show(noisyLena);
[h, w] = size(lena);
regulStrength = 10000;
NbIters = 100;
delta = 1;
% denoisedLena = SQSDenoising(noisyLena, regulStrength, NbIters, delta);
% figure(2); Show(denoisedLena);

noisyLena = repmat(noisyLena, [1, 1, 3]);

denoisedLena = SQSDenoising(noisyLena, [10 100 1], 7, [0.1 0.1 0.1]);
figure(2); Show(denoisedLena(:,:,1));
figure(3); Show(denoisedLena(:,:,2));
figure(4); Show(denoisedLena(:,:,3));