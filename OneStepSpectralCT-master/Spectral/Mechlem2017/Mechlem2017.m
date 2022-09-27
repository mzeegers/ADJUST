function [iterates, costs] = Mechlem2017(measurements, ObjectSize, NbPixelsPerProj, A, M, S, T, maxIter, lambda, NbSubsets, delta, runOnGPU)

[nrealmat, nmat] = size(T);

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
    diagS = gpuArray(diagS);
    T = gpuArray(T);
    lambda = gpuArray(lambda);
    delta = gpuArray(delta);
    alpha_k = zeros(ObjectSize^2, nmat, 'gpuArray');
    
    % Nesterov
    t_nesterov = ones(1,maxIter * NbSubsets, 'gpuArray');
    sum_nesterov = zeros(1,maxIter * NbSubsets + 1, 'gpuArray');
else
    alpha_k = zeros(ObjectSize^2, nmat);
    
    % Nesterov
    t_nesterov = ones(1,maxIter * NbSubsets);
    sum_nesterov = zeros(1,maxIter * NbSubsets + 1);
end

% Initialize, and precompute what can be precomputed
v_k = alpha_k;
z_k = alpha_k;
iterates = zeros(ObjectSize^2, nrealmat, maxIter);
costs = zeros(1, maxIter);

% Shuffle the projections order in the matrix A,
% for Ordered Subsets to be efficient
% Only then, copy the shuffled projection matrix 
% to the GPU (advanced indexing is not available on GPU)
NbProj = size(A, 1) / NbPixelsPerProj;
Shuffled = randperm(NbProj)-1;
NewIndices = kron(Shuffled.' * NbPixelsPerProj, ones(NbPixelsPerProj, 1));
NewIndices = NewIndices + kron(ones(NbProj, 1), (1:NbPixelsPerProj).');
ShuffledA = A(NewIndices, :);

% Shuffle the measurements as well, with the same permutation
ShuffledMeasurements = measurements(NewIndices, :);

% Compute the number of projections per subset
NbProjsPerSubset = ceil(NbProj / NbSubsets);

% Initialize Nesterov coefficients
for k=2:(maxIter * NbSubsets + 1)
    t_nesterov(1,k) = 0.5*(1+sqrt(1+ 4 * t_nesterov(1,k-1)^2));
    sum_nesterov(1,k) = sum_nesterov(1,k-1) + t_nesterov(1,k);
end
ratios_nesterov = t_nesterov ./ sum_nesterov;
ratios_nesterov(1, 1) = 0; % To remove the NaN, although it seems it is never used

% Build the submatrices once and for all
SubAs = cell(NbSubsets, 1);
SubMeasurements = cell(NbSubsets, 1);
Sum_aijs = cell(NbSubsets, 1);
for subset=1:NbSubsets
    % Compute the first and last row's index for the current subset
    % and extract the corresponding sub matrix for forward projection
    firstRow = (subset - 1) * NbProjsPerSubset * NbPixelsPerProj + 1;
    lastRow = min(subset * NbProjsPerSubset * NbPixelsPerProj, size(ShuffledA, 1));
    SubA = ShuffledA(firstRow:lastRow, :);

    if (runOnGPU)
        SubAs{subset} = gpuArray(SubA);
        SubMeasurements{subset} = gpuArray(ShuffledMeasurements(firstRow:lastRow, :));
        Sum_aijs{subset} = gpuArray(full(sum(SubA, 2)));
    else
        SubAs{subset} = SubA;
        SubMeasurements{subset} = ShuffledMeasurements(firstRow:lastRow, :);
        Sum_aijs{subset} = full(sum(SubA, 2));
    end
end

% Start iterations
for iter=1:maxIter
    disp(['In Mechlem2017, starting iteration ', num2str(iter)]);
    for subset=1:NbSubsets
        
        % Run the forward model
        [forw, forwmu, forwmumu] = WeidingerForwardModel(alpha_k, SubAs{subset}, M, S, diagS);
        ratios = SubMeasurements{subset}.' ./ (forw + eps);

        % Compute the gradient of the data-attachment surrogate function at alpha_k
        ToBackProject = (1 - ratios) .* forwmu;
        ToBackProject = squeeze(sum(ToBackProject, 1));
        gradients = SubAs{subset}.' * ToBackProject;
        
        % Compute the hessian of the surrogate function at alpha_k
        ToBackProject = squeeze(sum(forwmumu, 1));
        ToBackProject = ToBackProject .* Sum_aijs{subset};
        ToBackProject = reshape(ToBackProject, [size(ToBackProject, 1) nmat * nmat]);
        hessians = SubAs{subset}.' * ToBackProject;

        % - Compute the gradient and hessian of the regularization surrogate
        % - Perform one step of Newton's method to find the minimum of the
        % separable quadratic surrogate
        [regul_grad, regul_hess] = SQSregul(alpha_k, T, ObjectSize, lambda, delta, nrealmat, NbSubsets, 'huber');
        gradients = gradients + regul_grad * T;
        if(runOnGPU)
            g_k = GPU_GetNewtonUpdate(hessians,regul_hess, gradients, T);
        else
            g_k = GetNewtonUpdate(hessians,regul_hess, gradients, T);
        end

        % Compute index in Nesterov pre-computed coefficients vectors
        i = (iter - 1) * NbSubsets + subset;

        % Perform the actual Nesterov update
        alpha_k = z_k - g_k;
        v_k = v_k - t_nesterov(1, i) * g_k;
        z_k = alpha_k + ratios_nesterov(1, i+1) * (v_k - alpha_k);
    end

    % Store the iterate
    iterates(:,:,iter) = gather(alpha_k * T.');
    costs(iter) = SQScomputeCost(alpha_k, measurements, A, M, S, diagS, T, lambda, delta, nrealmat, ObjectSize, 'huber');
end

end