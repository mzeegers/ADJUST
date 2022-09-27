function [totalCost, data_attachment_cost, regularization_cost] = SQScomputeCost(x, measurements, A, M, S, diagS, T, lambda, delta, nrealmat, ObjectSize, method)
%SQSCOMPUTECOST Computes cost function of SQS methods
%   Computes the value of the cost function of SQS methods. 

    [NbBins, NbEnergies, NbPixelsPerProj] = size(S);
    NbPixels = size(A,1);
    NbProjs = NbPixels / NbPixelsPerProj;
   
    % Equivalent of S * attenuationFactors, which also works when S is 3D
    % (pixel-dependent incident spectrum, e.g. with a bowtie filter)
    attenuationFactors = exp( - M * (A * x).');
    expectedCounts = diagS * reshape(attenuationFactors, [NbEnergies * NbPixelsPerProj, NbProjs]);
    expectedCounts = reshape(expectedCounts, [NbBins NbPixels]).';
    
    % Data attachment
    data_attachment_cost = -sum(measurements(:) .* log(expectedCounts(:)) - expectedCounts(:));

    % Regularization
    tofilter = reshape(x * T.', [ObjectSize ObjectSize nrealmat]);
    diffs = SQSdiffWithNeighbors(tofilter);
    
    switch method
        case 'huber'
            regul_allmats = sum(Huber(diffs, delta, 0), 4);
        case 'hyperbola'
            regul_allmats = sum(hyperbola(diffs, delta, 0), 4);
        case 'green'
            regul_allmats = sum(GreenPrior(diffs, 0), 4);
        otherwise
            disp('Unknown approximation of L1 norm');
            return
    end

    regularization_cost = 0;
    for material = 1:nrealmat
        regul_mat = regul_allmats(:,:,material,:);
        regularization_cost = regularization_cost + lambda(material) * sum(regul_mat(:));
    end
    
    % Total
    totalCost = data_attachment_cost + regularization_cost;
    
    % Gather from GPU, if necessary (if already on CPU, gather does
    % nothing)
    totalCost = gather(totalCost);
    data_attachment_cost = gather(data_attachment_cost);
    regularization_cost = gather(regularization_cost);
end