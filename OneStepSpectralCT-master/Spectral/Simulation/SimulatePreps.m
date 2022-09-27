function [y, k_d, A, M, S, T, old_M,x] = SimulatePreps(ObjectSize, noise, getBeamTransmissionRatios, syntheticMaterials)
% Spectral photon counting scanner simulator
% [y, k_d, A, M, T, old_M] = SimulatePreps(ObjectSize, noise, getBeamTransmissionRatios, syntheticMaterials)
% computes simulated preps through a simple phantom of water, iodine
% and gadolinium. 
%
% y is the set of preps, either as photon counts or as transmissionRatio
% (if getBeamTransmissionRatios is false or true, respectively), with
% Poisson noise if "noise" is true and without otherwise
%
% k_d is the ratio between transmissionRatios and their variance (useful
% only for Cai)
%
% A is the forward projection matrix, which computes line integrals
%
% M is the transformed attenuation matrix, which depends on the choice of
% "syntheticMaterials":
% - if syntheticMaterials == 'none', M is the original attenuations matrix
% - if syntheticMaterials == 'normalize', M is the normalized attenuations
% matrix
% - if syntheticMaterials == 'orthonormalize', M is the orthonormalized 
% attenuations matrix
% - if syntheticMaterials == 'fessler', M is obtained by a trick
% described by Fessler in a 1999 patent
%
% T is the transform matrix from the original attenuations matrix to M. It
% is also the transform to apply to the concentrations of synthetic
% materials to get back to the concentrations of real materials, since the
% transforms on material concentrations and on matrices are the inverse of one another.
%
% old_M is the original attenuations matrix


% Import Philips simulated spectrum, detector response and material
% attenuations
load('ScannerModel/incidentSpectrum.mat', 'incidentSpectrum');
load('ScannerModel/materialAttenuations.mat', 'materialAttenuations');
load('ScannerModel/detectorResponse.mat', 'detectorResponse');

% Get them under matrix form, consistent with Cai's paper
M = materialAttenuations;
S_unbinned = detectorResponse .* incidentSpectrum.';

% Define thresholds and get S binned
S = cat(1, sum(S_unbinned(30:50, :), 1), ...
           sum(S_unbinned(51:61, :), 1), ...
           sum(S_unbinned(62:71, :), 1), ...
           sum(S_unbinned(72:82, :), 1), ...
           sum(S_unbinned(83:end, :), 1));

% Generate a simple phantom
unit = ObjectSize/8;
water=zeros(ObjectSize);
water(unit+1:7*unit, unit+1:7*unit) = ones(6*unit, 6*unit);
iodine=zeros(ObjectSize);
iodine(2*unit+1:3*unit, 2*unit+1:3*unit) = ones(unit, unit) * 0.01 / 4.933;
gadolinium=zeros(ObjectSize);
gadolinium(4*unit+1:5*unit, 5*unit+1:6*unit) = ones(unit, unit) * 0.01 / 7.9;

% writeImageResult(iodine, 'iodineGT.png', 0, 0.01);
% writeImageResult(gadolinium, 'gadoliniumGT.png', 0, 0.01);
% writeImageResult(water, 'waterGT.png', 0, 1);

% Compute the number of projections required
% to keep a pixels-to-voxels ratio of 4
NbProj = ceil(ObjectSize * 4 / sqrt(2));
ForwardTheta = [0:NbProj-1] * 180 / NbProj;

% Generate forward projection matrix
A = paralleltomo(ObjectSize,ForwardTheta);

% Compute forward projections
x = cat(2, iodine(:), gadolinium(:), water(:));

% Apply forward model
tic
[y, k_d] = NoisyForwardModel(x, A, M, S, noise, getBeamTransmissionRatios);
if (getBeamTransmissionRatios)
    % Normalize the joint spectra matrix S
    noAttenuation = S * ones(size(M,1),1);
    S = S ./ repmat(noAttenuation, [1, size(S, 2)]);
end
toc
% ShowProjs(y);

% Normalize the attenuations matrix
old_M = M;
T = eye(3);
if (strcmp(syntheticMaterials,'normalize'))
    [M, T] = normalize_attenuations(M);
    disp('Using normalized attenuations')
end
if (strcmp(syntheticMaterials,'orthonormalize'))
    [M, T] = orthonormalize_attenuations(M);
    disp('Using orthonormalized attenuations')
end
if (strcmp(syntheticMaterials,'fessler'))
    [M, T] = fesslers_trick(M);
    disp('Using Fessler trick')
end
if (strcmp(syntheticMaterials,'barber'))
    [M, T] = barber(M);
    disp('Using Barber trick')
end

% Check that T is indeed the transform from old_M to M
temp = (old_M * T - M);
if(sum(abs(temp(:))) > 1e-10)
    sum(abs(temp(:)))
    disp('Some problem occurred during orthonormalization');
end

    % Construct a basis of artificial materials, proposed by Fessler 
    % in a 1999 patent
    function [out, passage] = fesslers_trick(in)
        inv_passage = (S * in )./(S * ones(size(in))); % the beta_bar matrix in equation 71 of Fessler's patent
        passage = inv(inv_passage.' * inv_passage) * inv_passage.'; % The Moore-Penrose pseudo inverse, since the inv_passage matrix is not square
        out = in * passage;
        
%         % Alternative solution
%         [u,s,v] = svd((S * in )./(S * ones(size(in))), 0);
%         passage = v * inv(s) * u.';
    end

    % Construct a basis of artificial materials, proposed by Barber
    % in her 2016 paper
    function [out, passage] = barber(in)
        SymReal = in.' * in;
        [V,D] = eig(SymReal);
        P = sqrt(D) * V.'; 
        P_inv = V * sqrt(inv(D));
        out = in * P_inv;
        passage = P_inv;
    end

    % Very crude way of normalizing the attenuation matrix,
    % while keeping track of the transformation (to invert it later)
    function [normalized, passage] = normalize_attenuations(normalized)
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
    end

    % Very crude way of orthonormalizing the attenuation matrix,
    % while keeping track of the transformation (to invert it later)
    function [orthonormal, passage] = orthonormalize_attenuations(orthonormal)
        passage = eye(3);
        temp_passage_mat = eye(3); temp_passage_mat(1,1)= 1/norm(orthonormal(:,1));
        orthonormal = orthonormal * temp_passage_mat;
        passage = passage * temp_passage_mat;

        dotp = (orthonormal(:,2).' * orthonormal(:,1));
        temp_passage_mat = eye(3); temp_passage_mat(1,2)=-dotp;
        orthonormal = orthonormal * temp_passage_mat;
        passage = passage * temp_passage_mat;
        temp_passage_mat = eye(3); temp_passage_mat(2,2)= 1/norm(orthonormal(:,2));
        orthonormal = orthonormal * temp_passage_mat;
        passage = passage * temp_passage_mat;

        dotp1 = (orthonormal(:,3).' * orthonormal(:,1));
        dotp2 = (orthonormal(:,3).' * orthonormal(:,2));
        temp_passage_mat = eye(3); temp_passage_mat(1,3)=-dotp1;temp_passage_mat(2,3)=-dotp2;
        orthonormal = orthonormal * temp_passage_mat;
        passage = passage * temp_passage_mat;
        temp_passage_mat = eye(3); temp_passage_mat(3,3)= 1/norm(orthonormal(:,3));
        orthonormal = orthonormal * temp_passage_mat;
        passage = passage * temp_passage_mat;
    end
end