function [iterates, costs] = Long2014(Y, ObjectSize, NbPixelsPerProj, A, M, S, T, maxIter, lambda, NbSubsets, delta, runOnGPU)

% Initialize
nmat = size(M, 2);
nrealmat = size(T, 1);
nenergies = size(M, 1);
nbins = size(S, 1);

% Make S a sparse block-diagonal matrix, instead of a 3D matrix
cellS = cell(size(S,3), 1);
for pix=1:size(S,3)
    cellS{pix} = sparse(squeeze(S(:,:,pix)));
end
diagS = blkdiag(cellS{:});

if (runOnGPU)
    % Reset the GPU
    gpuDevice(1);
    
    % Copy inputs to GPU
    M = gpuArray(M);
    S = gpuArray(S);
    T = gpuArray(T);
    lambda = gpuArray(lambda);
    delta = gpuArray(delta);
    x_n = zeros(ObjectSize^2, nmat, 'gpuArray');
    mumu = zeros(nenergies, nmat, nmat, 'gpuArray');
else
    x_n = zeros(ObjectSize^2, nmat);
    mumu = zeros(nenergies, nmat, nmat);
end

% Shuffle the projections order in the matrix A,
% for Ordered Subsets to be efficient
NbProj = size(A, 1) / NbPixelsPerProj;
Shuffled = randperm(NbProj)-1;
NewIndices = kron(Shuffled.' * NbPixelsPerProj, ones(NbPixelsPerProj, 1));
NewIndices = NewIndices + kron(ones(NbProj, 1), (1:NbPixelsPerProj).');
ShuffledA = A(NewIndices, :);

% Shuffle the measurements as well, with the same permutation
ShuffledY = Y(NewIndices, :);

% Compute the number of subsets
NbProjsPerSubset = ceil(NbProj / NbSubsets);
storeEveryNthIterate = 10;
iterates = zeros(ObjectSize^2, nrealmat, ceil(maxIter/storeEveryNthIterate));
costs = zeros(1, maxIter);

% Pre-compute the mu(e)mu_t(e) matrix
for e=1:nenergies
    mumu(e,:,:) = M(e,:).' * M(e,:);
end
mumu = reshape(mumu, [nenergies, nmat * nmat]);

% Build the submatrices once and for all
SubAs = cell(NbSubsets, 1);
SubA2_pis = cell(NbSubsets, 1);
SubYs = cell(NbSubsets, 1);
for subset=1:NbSubsets
    % Compute the first and last row's index for the current subset
    % and extract the corresponding sub matrix for forward projection
    firstRow = (subset - 1) * NbProjsPerSubset * NbPixelsPerProj + 1;
    lastRow = min(subset * NbProjsPerSubset * NbPixelsPerProj, size(ShuffledA, 1));
    SubA = ShuffledA(firstRow:lastRow, :);
    SubA2_pi = SubA .* sum(SubA, 2); % Warning: implicit extension

    if (runOnGPU)
        SubAs{subset} = gpuArray(SubA);
        SubA2_pis{subset} = gpuArray(SubA2_pi);
        SubYs{subset} = gpuArray(ShuffledY(firstRow:lastRow, :));
    else
        SubAs{subset} = SubA;
        SubA2_pis{subset} = SubA2_pi;
        SubYs{subset} = ShuffledY(firstRow:lastRow, :);        
    end
end

for iter=1:maxIter
    disp(strcat('In Long and Fessler 2014, starting iteration ', num2str(iter)));
    
    for subset=1:NbSubsets
        nraysInSubset = size(SubAs{subset}, 1);

        % Compute the gradient of the data-fidelity term
        totalLinearAttenuation = SubAs{subset} * x_n * M.';
        attenuation_factors = exp(-totalLinearAttenuation).';
        
        % Equivalent of attenuation_factors * S.', which also works when 
        % S is 3D
        reshaped_af = reshape(attenuation_factors, [nenergies * NbPixelsPerProj, nraysInSubset/NbPixelsPerProj]);
        y_bar_im_s_im = reshape(diagS * reshaped_af, [nbins nraysInSubset]).';
        
        ratio = 1 - SubYs{subset} ./ y_bar_im_s_im;
        c = optimal_curvature(totalLinearAttenuation);

        % Equivalent of -S * af_M, which also works when 
        % S is 3D
        af_M = reshape(attenuation_factors, [nenergies, nraysInSubset, 1]) .* reshape(M, [nenergies, 1, nmat]); % Warning: implicit extension
        grad_si_y_bar_im_s_imn = -diagS * reshape(af_M, [nenergies * NbPixelsPerProj, nraysInSubset/NbPixelsPerProj * nmat]);
        grad_si_y_bar_im_s_imn = reshape(grad_si_y_bar_im_s_imn, [nbins nraysInSubset nmat]);
        insum = squeeze(sum(ratio.' .* grad_si_y_bar_im_s_imn, 1));
        
        % Compute the curvature matrix
        reshaped_c = reshape(c, [1 NbPixelsPerProj nraysInSubset/NbPixelsPerProj nenergies]);
        reshaped_mumu = reshape(mumu, [1 1 1 nenergies nmat^2]);
        c_mumu = reshaped_c .* reshaped_mumu;
        c_mumu = permute(c_mumu, [4 2 1 3 5]);
        BigC = diagS * reshape(c_mumu, [nenergies * NbPixelsPerProj, nraysInSubset/NbPixelsPerProj * nmat^2]);
        BigC = reshape(BigC, [nbins nraysInSubset nmat^2]);
        BigC = squeeze(sum(BigC, 1));

        L_dot_js = SubAs{subset}.' * insum;

        % Compute the hessian of the data-fidelity term
        hessians_data = SubA2_pis{subset}' * BigC;

        % - Compute the gradient and hessian of the regularization term
        % - Perform one step of Newton's method to find the minimum of the
        % separable quadratic surrogate
        [R_dot_js, hessian_regul] = SQSregul(x_n, T, ObjectSize, lambda, delta, nrealmat, NbSubsets, 'hyperbola');
        gradients = L_dot_js + R_dot_js * T;
        if(runOnGPU)
            x_n = x_n - GPU_GetNewtonUpdate(hessians_data, hessian_regul, gradients, T);
        else
            x_n = x_n - GetNewtonUpdate(hessians_data, hessian_regul, gradients, T);
        end
    end
    
    % Store the iterate and cost
    if (mod(iter, storeEveryNthIterate)==0)
        iterates(:,:,round(iter/storeEveryNthIterate)) = gather(x_n * T.');
    end
    costs(iter) = SQScomputeCost(x_n, Y, A, M, S, diagS, T, lambda, delta, nrealmat, ObjectSize, 'hyperbola');
end

    function out = optimal_curvature(in)
        out = 2 * (1 - exp(-in) - in .* exp(-in))./(in.^2);
        out(abs(in)<eps) = 1;
        out(out<0)=0;
    end

end