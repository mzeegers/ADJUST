function updates = GetNewtonUpdate(data_hessians,regul_hessians, gradients, T)
%GETNEWTONUPDATE Computes the update step of Newton's method
%   Computes the update step of Newton's method, voxel-per-voxel, 
%   given the hessians and the gradients

    nmat = size(T, 2);
    updates = zeros(size(gradients));

    % Compute f_n+1, voxel by voxel
    for vox=1:size(data_hessians, 1)
        H = reshape(data_hessians(vox, :), [nmat nmat]) + T.' * diag(regul_hessians(vox, :)) * T;
        updates(vox, :) = gradients(vox, :) / H; % Equivalent to "* inv(H)" but faster
    end

end

