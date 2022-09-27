function updates = GPU_GetNewtonUpdate(data_hessians,regul_hessians, gradients, T)
%GPU_GETNEWTONUPDATE Computes the update step of Newton's method on GPU
%   Computes the update step of Newton's method on GPU, using pagefun to 
%   invert many small hessians

    [nrealmat, nmat] = size(T);
    nvox = size(gradients, 1);

    % Equivalent to calling diag on each row of regul_hess
    regul_hess_for_pagefun = regul_hessians(:, 1);
    for j=2:nrealmat
        regul_hess_for_pagefun = cat(2, regul_hess_for_pagefun, zeros(nvox, nrealmat, 'gpuArray'), regul_hessians(:,j));
    end
    
    % Prepare for pagefun multiplication by T and T.'
    regul_hess_for_pagefun = reshape(regul_hess_for_pagefun.', [nrealmat nrealmat nvox]);
    T_for_pagefun = repmat(T, [1 1 nvox]);
    Tt_for_pagefun = repmat(T.', [1 1 nvox]);
    
    % Perform multiplications
    regul_hess_for_pagefun = pagefun(@mtimes, regul_hess_for_pagefun, T_for_pagefun);
    regul_hess_for_pagefun = pagefun(@mtimes, Tt_for_pagefun, regul_hess_for_pagefun);
    
    % Prepare full hessian for pagefun
    data_hess_for_pagefun = reshape(data_hessians.', [nmat nmat nvox]);
    H_for_pagefun = data_hess_for_pagefun  + regul_hess_for_pagefun;
    
    % Prepare gradients for pagefun
    gradients_for_pagefun = reshape(gradients.', [nmat 1 nvox]);
    
    % Compute update by pagefun
    updates = squeeze(pagefun(@mldivide, H_for_pagefun, gradients_for_pagefun)).';
end

