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

strOut = load('example_Thorax_nonInversCrimeFinal.mat').strOut;
%%
strOut.strProb.A = A;
strOut.strProb.F = F;

save('example_Thorax_nonInversCrime.mat', 'strOut');