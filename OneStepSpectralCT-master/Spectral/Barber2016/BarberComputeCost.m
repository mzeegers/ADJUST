function [data_attachment_cost] = BarberComputeCost(x, measurements, A, M, S, diagS)
%BARBERCOMPUTECOST Computes cost function of Barber2016
%   Computes the value of the cost function of Barber2016. 

[NbBins, NbEnergies, NbPixelsPerProj] = size(S);
NbPixels = size(A,1);
NbProjs = NbPixels / NbPixelsPerProj;

% Equivalent of S * attenuationFactors, which also works when S is 3D
% (pixel-dependent incident spectrum, e.g. with a bowtie filter)
attenuationFactors = exp( - M * (A * x).');
expectedCounts = diagS * reshape(attenuationFactors, [NbEnergies * NbPixelsPerProj, NbProjs]);
expectedCounts = reshape(expectedCounts, [NbBins NbPixels]).';

% Data attachment
data_attachment_cost = gather(-sum(measurements(:) .* log(expectedCounts(:)) - expectedCounts(:)));

end