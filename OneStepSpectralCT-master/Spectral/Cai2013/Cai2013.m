function [iterates, costs] = Cai2013(y, ObjectSize, A, M, S, k_d, maxIter, T, regulWeights, delta_h, runOnGPU, useNesterov)

% Initialize conjugate gradient
nmat = size(M, 2);
nrealmat = size(T, 1);
[NbBins, NbEnergies, NbPixelsPerProj] = size(S);
NbRays = size(A, 1);

if (runOnGPU)
    % Reset the GPU
    gpuDevice(1);
    
    % Copy inputs to GPU
    y = gpuArray(y);
    A = gpuArray(A);
    M = gpuArray(M);
    S = gpuArray(S);
    T = gpuArray(T);
    regulWeights = gpuArray(regulWeights);
    delta_h = gpuArray(delta_h);

    x_k = zeros(ObjectSize^2, nmat, 'gpuArray');
else
    x_k = zeros(ObjectSize^2, nmat);
end

d_k = x_k;
g_k_minus_one = x_k;

% Initialize Nesterov intermediate variables
if (useNesterov)
    t_nesterov = ones(1,maxIter, 'like', x_k);
    sum_nesterov = zeros(1,maxIter + 1, 'like', x_k);
    v_k = x_k;
    z_k = x_k;
    
    for k=2:(maxIter + 1)
        t_nesterov(1,k) = 0.5*(1+sqrt(1+ 4 * t_nesterov(1,k-1)^2));
        sum_nesterov(1,k) = sum_nesterov(1,k-1) + t_nesterov(1,k);
    end
    ratios_nesterov = t_nesterov ./ sum_nesterov;
    ratios_nesterov(1, 1) = 0; % To remove the NaN, although it seems it is never used
end

% Initialize iterates and costs
% They contain maxIter or maxIter/storeEveryNthIterate elements, not including the initial guess
storeEveryNthIterate = 10;
iterates = zeros(ObjectSize^2, nrealmat, ceil(maxIter/storeEveryNthIterate));
costs = zeros(1, maxIter);
[y_bar, Q] = ForwardModel(x_k);
descentDirectionJustReset = false;

% sigma_x = sqrt(1./regulWeights); % Not used, just reminded here to show
% the link between the variable used in the paper and the one implemented

for iter=1:maxIter
    disp(strcat('In Cai2013, starting iteration ', num2str(iter)));

    % The forward model is computed at the end of the iteration,
    % so it can be used to compute the cost function and determine
    % whether the heuristic step size gives a descent step
    
    % Compute the gradient
    z = (y.^2 - y_bar.^2) ./ (k_d * y_bar.^2) - 1./y_bar; % NB: there is a typo in the paper, only the appendix is correct
    zS = reshape(z, [NbPixelsPerProj NbRays/NbPixelsPerProj NbBins]) .* permute(S, [3 4 1 2]);
    zS = reshape(sum(zS, 3), [NbRays NbEnergies]);
    g_k = A.' * (Q .* zS) * M;
    
    % Regularize, taking into account passage matrix to real volumes
    % and its transpose
    x_reshaped = reshape(x_k * T.', [ObjectSize, ObjectSize, nrealmat]);
    regul = div_space(Huber(grad_space(x_reshaped), delta_h, 1));
    regul_reshaped = reshape(regul, [ObjectSize^2, nrealmat]);
    for mat=1:nrealmat
        regul_reshaped(:,mat) = - regulWeights(mat) * regul_reshaped(:,mat);
    end
    g_k = g_k + regul_reshaped * T;

    % Compute the descent direction
    if (iter==1)
        beta = 0;
    else
        beta = dot(g_k(:) - g_k_minus_one(:), g_k_minus_one(:)) / dot(g_k_minus_one(:), g_k_minus_one(:));
    end
    if (beta<0)
        beta = 0;
    end
    d_k = -g_k + beta * d_k;
    g_k_minus_one = g_k;

    % Compute the descent step
    % Start by computing d_kt_H_d_k
    v = (2 * y.^2) ./ (k_d * y_bar.^3) - 1./y_bar.^2;
    Bd = A * d_k * M.';
    QSt = reshape(Q, [NbPixelsPerProj NbRays/NbPixelsPerProj NbEnergies]) .* permute(S, [3 4 2 1]);
    QSt = reshape(sum(QSt, 3), [NbRays NbBins]);
    brace = reshape(v .* QSt - z, [NbPixelsPerProj NbRays/NbPixelsPerProj NbBins]) .* permute(S, [3 4 1 2]);
    brace = reshape(sum(brace, 3), [NbRays NbEnergies]);
    Q_brace_Bd = Q .* brace .* Bd;
    if (sum(isnan(Q_brace_Bd(:))) > 0)
        disp('NaN detected');
    end
    Bt_Q_brace_Bd = A.' * Q_brace_Bd * M;
    dt_H_d = dot(d_k(:),Bt_Q_brace_Bd(:));
        
    % Compute the part due to regularization, again taking into account 
    % the passage matrix to real volumes and its transpose
    Td = d_k * T.';
    Td_reshaped = reshape(Td, [ObjectSize, ObjectSize, nrealmat]);
    Tx_reshaped = reshape(x_k * T.', [ObjectSize, ObjectSize, nrealmat]);
    regul = Td_reshaped .* div_space(Huber(grad_space(Tx_reshaped), delta_h, 2) .* grad_space(Td_reshaped));
    regul_reshaped = reshape(regul, [ObjectSize^2, nrealmat]);
    for mat=1:nrealmat
        regul_reshaped(:,mat) = - regulWeights(mat) * regul_reshaped(:,mat);
    end
    dt_H2_d = regul_reshaped * T;
    dt_H_d = dt_H_d + sum(dt_H2_d(:));

    alpha_k = - sum(g_k(:) .* d_k(:)) / dt_H_d;
    candidateFound = (iter==1);
    u_k = -alpha_k * d_k;

    % Compute the cost function for this candidate alpha. If it is larger than the
    % current cost function, divide alpha_k by 2 and try again
    % If we end up dividing alpha_k by 2 too many times, the problem is
    % probably the direction rather than the step length. In that case,
    % reinitialize the direction to -g_k
    NbTimesAlphaHalved = 0;
    while(not(candidateFound))
        
        % Compute the potential update
        u_k = -alpha_k * d_k;

        if (useNesterov)
            % Get the candidate
            x_kCandidate = z_k - u_k;
        else
            x_kCandidate = x_k - u_k;
        end
 
        % Compute the cost in order to determine whether the candidate
        % actually decreases the cost function or not
        [y_bar, Q] = ForwardModel(x_kCandidate);
        candidateCost = ComputeCost(y_bar, x_kCandidate);
        if (candidateCost < costs(iter-1))
            costs(iter) = gather(candidateCost);
            candidateFound = true;
            x_k = x_kCandidate;
            descentDirectionJustReset = false;
        else
            disp(strcat('In iteration ', num2str(iter), '. Dividing alpha by 2'));
            alpha_k = alpha_k / 2;
            NbTimesAlphaHalved = NbTimesAlphaHalved + 1;
        end
        if (NbTimesAlphaHalved==10)
            if(descentDirectionJustReset)
                return
            end
            alpha_k = 0; 
            candidateFound = true;
            d_k = -g_k;
            descentDirectionJustReset = true;
        end
    end
    
    if (useNesterov)
        % Compute the rest of the Nesterov's update
        v_k = v_k - t_nesterov(1, iter) * u_k;
        z_k = x_k + ratios_nesterov(1, iter+1) * (v_k - x_k);
    end
    
    % Store the iterate and cost
    if (mod(iter, storeEveryNthIterate)==0)
        iterates(:,:,round(iter/storeEveryNthIterate)) = gather(x_k * T.');
    end
    costs(iter) = gather(ComputeCost(y_bar, x_k));
end

    function [cost] = ComputeCost(y_bar, x)
        % Data-attachment
        cost = log(y_bar) + (y_bar - y).^2 ./ (k_d * y_bar);
        cost = sum(cost(:));
        
        % Regularization
        xs = reshape(x * T.', [ObjectSize, ObjectSize, nrealmat]);
        regul_allmats = Huber(grad_space(xs), delta_h, 0);
        for material = 1:nrealmat
            regul_mat = regul_allmats(:,:,material,:);
            cost = cost + regulWeights(material) * sum(regul_mat(:));
        end
        
    end

    function [y_bar, e_mB_x] = ForwardModel(x)
        % Compute forward projections
        projs = A * x;
        
        % Compute photon counts as seen by the detector
        attenuationFactors = exp(- projs * M.');  % attenuationFactors
        e_mB_x = attenuationFactors;
        y_bar = reshape(attenuationFactors, [NbPixelsPerProj NbRays/NbPixelsPerProj NbEnergies]) .* permute(S, [3 4 2 1]);
        y_bar = reshape(sum(y_bar, 3), [NbRays NbBins]);
    end

end
