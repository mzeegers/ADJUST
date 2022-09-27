function [forw, forwmu, forwmumu] = WeidingerForwardModel(x, A, M, S, diagS)
%GPU_WEIDINGERFORWARDMODEL Computes the forward model for Weidinger2016 and
%Mechlem2017
    
    % Various sizes
    NbPixels = size(A, 1);
    [NbBins, NbEnergies, NbPixelsPerProj] = size(S);
    NbProjs = NbPixels / NbPixelsPerProj;
    NbMats = size(M,2);

    % Compute forward projections
    projs = A * x;
    attenuationFactors = exp( - M * projs.');
    
    % Equivalent of S * attenuationFactors, which also works when S is 3D
    % (pixel-dependent incident spectrum, e.g. with a bowtie filter)
    reshaped_att = reshape(attenuationFactors, [NbEnergies * NbPixelsPerProj NbProjs]);
    forw = diagS * reshaped_att;
    forw = reshape(forw, [NbBins NbPixels]);
    
    % Compute the forwmu, i.e. the part of the gradient described in
    % equation (30). Again, equivalent of -S * attenuationFactorsMu, 
    % which works when S is 3D
    attenuationFactorsMu = attenuationFactors .* reshape(M, [NbEnergies 1 NbMats]);
    forwmu = -diagS * reshape(attenuationFactorsMu, [NbEnergies * NbPixelsPerProj, NbProjs * NbMats]);
    forwmu = reshape(forwmu, [NbBins NbPixels NbMats]);
    
    % Compute the forwmumu, i.e. the integral part of the hessian
    % (see equation (28)). Again, equivalent of S * attenuationFactorsMuMu, 
    % which works when S is 3D
    attenuationFactorsMuMu = reshape(M, [NbEnergies 1 1 NbMats]) .* attenuationFactorsMu;
    forwmumu = diagS * reshape(attenuationFactorsMuMu, [NbEnergies * NbPixelsPerProj, NbProjs * NbMats * NbMats]);
    forwmumu = reshape(forwmumu, [NbBins NbPixels NbMats NbMats]);
end