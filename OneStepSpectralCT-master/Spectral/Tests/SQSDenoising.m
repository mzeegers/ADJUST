function f = SQSDenoising(f_0, lambda, maxIter, delta)

    [ObjectSize, ~, nchannels] = size(f_0);
    f_0 = reshape(f_0, [ObjectSize^2, nchannels]);
    f = f_0;
       
    % Start iterations
    for iter=1:maxIter
        disp(strcat('In SQSDenoising, starting iteration ', num2str(iter)));

        [regul_grad, regul_hess, regul_cost] = SQSregul(f, eye(nchannels), ObjectSize, lambda, delta, nchannels, 1, 'huber');
        disp(['Regularization term is ', num2str(regul_cost)]);
        
        % Compute the gradient of the surrogate function at alpha_k
        data_attachment_grad = 2 * (f - reshape(f_0, [ObjectSize^2, nchannels]));
        g = data_attachment_grad + regul_grad;

        % Compute the hessian of the surrogate function at alpha_k
        data_attachment_hess = 2 * ones(size(f));
        h = data_attachment_hess + regul_hess;

        % Compute update
        update = -g./h;

        % Compute f_k+1
        f = f + update;

        % Check if NaN
        if(sum(isnan(f(:)))>0)
            return
        end
    end
    
    f = reshape(f, [ObjectSize,ObjectSize, nchannels]);
end