function [iterates, gaps, costs] = Barber2016(measurements, ObjectSize, X, M, S, T, maxIter, lambda, theta, gammas, runOnGPU)

% Initialize, and precompute what can be precomputed
N_l = size(X, 1); % Number of pixels
N_w = size(measurements, 2); % Number of bins
N_i = size(S, 2); % Number of energies
nmat = size(M, 2); % Number of (possibly synthetic) materials, also denoted N_m
nrealmat = size(T, 1); % Number of real materials
NbPixelsPerProj = size(S,3);

% Make S a sparse block-diagonal matrix, instead of a 3D matrix
cellS = cell(NbPixelsPerProj, 1);
for pix=1:NbPixelsPerProj
    cellS{pix} = sparse(squeeze(S(:,:,pix)));
end
diagS = blkdiag(cellS{:});

if (runOnGPU)
    % Reset the GPU
    gpuDevice(1);
        
    % Copy inputs onto GPU
    measurements = gpuArray(measurements);
    X = gpuArray(X);
    M = gpuArray(M);
    S = gpuArray(S);
    T = gpuArray(T);
    
    f_n = zeros(ObjectSize^2, nmat, 'gpuArray');
    y_sino_n = zeros(size(measurements), 'gpuArray');
    y_grad_n = zeros(ObjectSize, ObjectSize, nmat, 2, 'gpuArray');
else
    f_n = zeros(ObjectSize^2, nmat);
    y_sino_n = zeros(size(measurements));
    y_grad_n = zeros(ObjectSize, ObjectSize, nmat, 2);
end
barf_n = f_n;
barf_nm1 = f_n;
y_sino_nm1 = y_sino_n;

storeEveryNthIterate = 10;
iterates = zeros(ObjectSize^2, nrealmat, ceil(maxIter/storeEveryNthIterate));
gaps = zeros(1,maxIter);
primals = zeros(1,maxIter);
duals = zeros(1,maxIter);
costs = zeros(1, maxIter);

% Start iterations
for iter=1:maxIter
    disp(strcat('In Barber2016, starting iteration ', num2str(iter)));
    
    f_0 = barf_n;

    % Compute Z f_0
    Zf_0 = M * (X * f_0)';
    AF = softexp(-Zf_0); %Attenuation Factors
    reshaped_AF = reshape(AF, [N_i*NbPixelsPerProj N_l/NbPixelsPerProj]);
    reshaped_EC = diagS * reshaped_AF; %Expected counts
    EC = reshape(reshaped_EC, [N_w N_l])';
    numeratorA = softexp_dot(-Zf_0);
%     ShowProjs(EC)
        
    % Compute Sigma
    Sigma_sino = K_times(ones(size(f_0)), numeratorA, EC, true);
    Sigma_grad = absgrad(reshape(ones(size(f_0)) * T.', [ObjectSize ObjectSize nmat])); % WRONG if T has negative values
    Sigma_sino = 1 ./ (lambda * Sigma_sino + eps);
    Sigma_grad = 1 ./ (lambda * Sigma_grad + eps);
    W = sqrt(Sigma_grad);
    
    Tho = Kt_times(ones(size(measurements)), numeratorA, EC, true) + reshape(absdiv(ones(size(Sigma_grad))), [ObjectSize^2 nmat]) * T.'; % WRONG if T is not diagonal
    Tho = lambda ./ (Tho + eps);
    
    % Compute z_0
    K_barf_nm1 = K_times(barf_nm1, numeratorA, EC, false);
    z_0 = (y_sino_nm1 - y_sino_n + Sigma_sino .* K_barf_nm1) ./ Sigma_sino;
    
    % Compute y_sino_np1
    r = measurements - EC;
    E = max(-r, zeros(size(measurements)));
    b = (EC - E) .* K_times(f_0, numeratorA, EC, false) - r;
    y_sino_np1 = (EC .* (y_sino_n + Sigma_sino .* K_times(barf_n, numeratorA, EC, false)) - Sigma_sino .* (b + E .* z_0) ) ./ (EC + Sigma_sino);
    
    % Compute y_grad_np1
    y_plus_grad = y_grad_n + Sigma_grad .* grad(reshape(barf_n * T.', [ObjectSize ObjectSize nmat]));
    g_plus = abs(y_plus_grad ./ (W + eps));
    hat_g_plus = (y_plus_grad ./ (W + eps)) ./ (g_plus + eps);
    y_grad_np1 = y_plus_grad - W .* hat_g_plus .* projOntoWeightedL1Ball(g_plus, W, gammas);
    
    % Compute f_np1
    f_np1 = barf_n - Tho .* (Kt_times(y_sino_np1, numeratorA, EC, false) - reshape(div(y_grad_np1), [ObjectSize^2 nmat]) * T.');
    
    % Compute barf_np1
    % Chambolle-Pock notation. Theta must be between 0 and 1. In Barber 2016, it is set to 1
    barf_np1 = f_np1 + theta * (f_np1 - f_n);
    
    % Compute the primal-dual gap to check for convergence
    Kf = K_times(barf_np1, numeratorA, EC, false);
    v = y_sino_np1(:) + b(:) + E(:) .* z_0(:);
    primals(iter) = gather(dot(Kf(:), 0.5 * EC(:) .* Kf(:) - b(:) - E(:) .* z_0(:)) + dot(z_0(:), b(:)));
    duals(iter) = gather(0.5 * sum(v ./ EC(:) .* v));
    gaps(iter) = primals(iter) - duals(iter);
    
    % Update iteration number on temporary variables
    f_n = f_np1;
    barf_nm1 = barf_n;
    barf_n = barf_np1;
    y_sino_nm1 = y_sino_n;
    y_sino_n = y_sino_np1;
    
    % Store the iterate
    if (mod(iter, storeEveryNthIterate)==0)
        iterates(:,:,round(iter/storeEveryNthIterate)) = gather(f_n * T.');
    end
    costs(iter) = BarberComputeCost(f_n, measurements, X, M, S, diagS);
    
    % Check if NaN
    if(sum(isnan(f_n(:)))>0)
        return
    end
end

% figure(2); 
% hold off
% plot(primals);
% hold on
% plot(duals);
% legend('primal', 'dual');

    function out = K_times(in, af, ec, abso)
        X_in = X * in;
        SubAMs = GetSubAMs(af, ec, abso);
        
        if (runOnGPU)
            % Prepare pagefun
            argu1 = permute(SubAMs, [1 3 2]);
            argu2 = reshape(X_in.', [nmat 1 N_l]);
            out = squeeze(pagefun(@mtimes, argu1, argu2)).';
        else
            out = zeros(N_l, N_w);
            % Compute pixelwise
            for l=1:N_l
                % Multiply by SubAM
                out(l, :) = X_in(l,:) * squeeze(SubAMs(:,l,:)).';
            end
        end
    end

    function out = Kt_times(in, af, ec, abso)
        SubAMs = GetSubAMs(af, ec, abso); 

        if (runOnGPU)
            % Prepare pagefun
            argu2 = permute(SubAMs, [1 3 2]);
            argu1 = reshape(in.', [1, N_w, N_l]);
            proj_out = squeeze(pagefun(@mtimes, argu1, argu2)).';
        else
            proj_out = zeros(N_l, nmat);
            for l=1:N_l
                % Multiply by SubAM
                proj_out(l,:) = in(l,:) * squeeze(SubAMs(:,l,:));
            end            
        end
        % Compute back projection
        out = X.' * proj_out;
    end

    function SubAMs = GetSubAMs(af, ec, abso)
        % Vectorized way of computing SubAMs
        SubAs = repmat(permute(S, [2 3 1]), [1 N_l/NbPixelsPerProj 1]) .* reshape(af, [N_i, N_l, 1]) ./ reshape(ec, [1, N_l, N_w]); % Warning: implicit extensions
        temp = reshape(permute(SubAs, [3 2 1]), [N_w * N_l, N_i]) * M;
        SubAMs = reshape(temp, [N_w, N_l, nmat]);
        
        % For step lengths, an absolute value is required
        if (abso)
            SubAMs = abs(SubAMs);
        end
    end

    function res = softexp(x)
        mask = (x < 0);
        res = exp(x) .* mask + (1 + x) .* (1 - mask);
    end

    function res = softexp_dot(x)
        mask = (x < 0);
        res = exp(x) .* mask + (1 - mask);
    end

    function g = grad(in)
        [h,w,~] = size(in);
        
        % Pad with zeros, then take the gradient
        pad_in = zeros(h+2,w+2,nrealmat, 'like', in);
        pad_in(2:end-1,2:end-1,:) = in;
        pad_gx = zeros(size(pad_in), 'like', pad_in);
        pad_gy = zeros(size(pad_in), 'like', pad_in);
        
        % Process one material at a time
        for m=1:size(in, 3)
            [pad_gx(:,:,m), pad_gy(:,:,m)] = gradient(pad_in(:,:,m));
        end
        g = cat(4, pad_gx(2:end-1, 2:end-1, :), pad_gy(2:end-1, 2:end-1, :));
    end

    function res = div(in)
        [h,w,~] = size(in);
        
        % Pad with zeros, then take the divergence
        pad_in = zeros(h+2,w+2,nrealmat, 2, 'like', in);
        pad_in(2:end-1,2:end-1,:,:) = in;
        pad_res = zeros(size(squeeze(pad_in(:,:,:,1))), 'like', in);
        
        % Process one material at a time
        for m=1:size(in, 3)
            [pad_res(:,:,m)] = divergence(squeeze(pad_in(:,:,m, 1)), squeeze(pad_in(:,:,m, 2)));
        end
        res = pad_res(2:end-1, 2:end-1, :);
    end

    function g = absgrad(in)
        [h,w,~] = size(in);
        
        % Pad with zeros, then take the gradient
        pad_in = zeros(h+2,w+2,nrealmat, 'like', in);
        pad_in(2:end-1,2:end-1,:) = in;
        pad_gx = zeros(size(pad_in), 'like', pad_in);
        pad_gy = zeros(size(pad_in), 'like', pad_in);
        
        % Process one material at a time
        for m=1:size(in, 3)
            % Same as [pad_gx, pad_gy]=gradient(pad_in), but with only
            % positive coefficients
            pad_gx(:,:,m) = conv2(pad_in(:,:,m), [0.5, 0, 0.5], 'same');
            pad_gy(:,:,m) = conv2(pad_in(:,:,m), [0.5; 0; 0.5], 'same');
        end
        g = cat(4, pad_gx(2:end-1, 2:end-1, :), pad_gy(2:end-1, 2:end-1, :));
    end

    function res = absdiv(in)
        [h,w,~] = size(in);
        
        % Pad with zeros, then take the divergence
        pad_in = zeros(h+2,w+2,nrealmat, 2, 'like', in);
        pad_in(2:end-1,2:end-1,:,:) = in;
        pad_res = zeros(size(squeeze(pad_in(:,:,:,1))), 'like', in);
        
        % Process one material at a time
        for m=1:size(in, 3)
            [pad_res(:,:,m)] = conv2(squeeze(pad_in(:,:,m, 1)), [0.5, 0, 0.5], 'same') ...
                             + conv2(squeeze(pad_in(:,:,m, 2)), [0.5; 0; 0.5], 'same');
        end
        res = pad_res(2:end-1, 2:end-1, :);
    end


    function out = projOntoWeightedL1Ball(g, weights, gammas)
        % Initialize the output
        out = zeros(size(g), 'like', g);
        
        for material = 1:size(g, 3)
            % Select the current material
            g_m = g(:,:,material, :);
            g_m = g_m(:);
            w_m = weights(:,:,material, :);
            w_m = w_m(:);
            
            % If we're already in the ball, do nothing
            l1norm = sum(abs(g_m ./ (w_m + eps)));
            if( l1norm <= gammas(material))
                out(:,:,material, :) = g(:,:,material, :);
                
            % Otherwise, project onto the sphere
            else
                disp('TV constraint activated');
                
                % Find the correct alpha by dichotomy
                min_alpha = 0;
                max_alpha = max(g_m .* w_m);
                current_alpha = (max_alpha + min_alpha)/2;

                % Run 20 iterations of dichotomic search
                for dich_iter=1:20
                    current_g_m = max(g_m - current_alpha ./ (w_m + eps), 0);
                    l1norm = sum(current_g_m);

                    if (l1norm > gammas(material)) % Need larger alpha
                        min_alpha = current_alpha;
                    else
                        max_alpha = current_alpha; % Need smaller alpha
                    end
                    current_alpha = (max_alpha + min_alpha)/2;
                end
                current_g_m = max(g_m - current_alpha ./ (w_m + eps), 0);
                out(:,:,material, :) = reshape(current_g_m, [ObjectSize ObjectSize 1 2]);
            end
        end
    end


end