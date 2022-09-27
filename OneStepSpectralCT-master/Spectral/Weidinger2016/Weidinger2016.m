function [iterates, costs] = Weidinger2016(measurements, ObjectSize, A, M, S, T, maxIter, lambda, runOnGPU)

% Initialize, and precompute what can be precomputed
nmat = size(M, 2);
nrealmat = size(T, 1);
sum_aij = full(sum(A, 2));
storeEveryNthIterate = 10;
iterates = zeros(ObjectSize^2, nrealmat, ceil(maxIter/storeEveryNthIterate));
costs = zeros(1, maxIter);

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
    measurements = gpuArray(measurements);
    A = gpuArray(A);
    M = gpuArray(M);
    S = gpuArray(S);
    diagS = gpuArray(diagS);
    T = gpuArray(T);
    lambda = gpuArray(lambda);
    f_n = zeros(ObjectSize^2, nmat, 'gpuArray');
else
    f_n = zeros(ObjectSize^2, nmat);
end

% Start iterations
for iter=1:maxIter
    disp(strcat('In Weidinger2016, starting iteration ', num2str(iter)));

    % Every three hundred iterations, defragment the GPU memory
    if((mod(iter, 300)==0) && runOnGPU)
        defragGPU(who());
    end
    
    [forw, forwmu, forwmumu] = WeidingerForwardModel(f_n, A, M, S, diagS);
    ratios = measurements.' ./ forw;
    
    % Compute the gradient of the surrogate function at f_n
    ToBackProject = (1 - ratios) .* forwmu;
    ToBackProject = squeeze(sum(ToBackProject, 1));
    gradients = A.' * ToBackProject;
     
    % Compute the hessian of the surrogate function at f_n
    ToBackProject = squeeze(sum(forwmumu, 1));
    ToBackProject = ToBackProject .* sum_aij;
    ToBackProject = reshape(ToBackProject, [size(ToBackProject, 1) nmat * nmat]);
    hessians = A.' * ToBackProject;

    % Compute the gradient and hessian of the regularization surrogate
    % Perform one step of Newton's method to find the minimum of the
    % separable quadratic surrogate
    % Store the iterate and cost
    [regul_grad, regul_hess] = SQSregul(f_n, T, ObjectSize, lambda, 0, nrealmat, 1, 'green');
    gradients = gradients + regul_grad * T;
    if (runOnGPU)
        f_n = f_n - GPU_GetNewtonUpdate(hessians,regul_hess, gradients, T);
    else
        f_n = f_n - GetNewtonUpdate(hessians,regul_hess, gradients, T);
    end
    
    % Store one iterate every few
    if (mod(iter, storeEveryNthIterate)==0)
        iterates(:,:,round(iter/storeEveryNthIterate)) = gather(f_n * T.');
    end
    
    % Store the cost
    costs(iter) = SQScomputeCost(f_n, measurements, A, M, S, diagS, T, lambda, 0, nrealmat, ObjectSize, 'green');
end

end