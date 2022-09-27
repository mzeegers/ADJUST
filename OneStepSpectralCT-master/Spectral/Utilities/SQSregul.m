function [grad_cost, hess_cost, value_cost] = SQSregul(input, T, ObjectSize, lambda, delta, nrealmat, NbSubsets, method)
%SQSREGUL Regularization function used in all SQS methods

    % Compute the gradient and the hessian of the regularization surrogate
    tofilter = reshape(input * T.', [ObjectSize ObjectSize nrealmat]);
    diffs = SQSdiffWithNeighbors(tofilter);

    switch method
        case 'huber'
            value_cost = sum(Huber(diffs, delta, 0), 4);
            grad_cost = sum(Huber(diffs, delta, 1), 4);
            hess_cost = sum(Huber(diffs, delta, 2), 4);
        case 'hyperbola'
            value_cost = sum(hyperbola(diffs, delta, 0), 4);
            grad_cost = sum(hyperbola(diffs, delta, 1), 4);
            hess_cost = sum(hyperbola(diffs, delta, 2), 4);
        case 'green'
            value_cost = sum(GreenPrior(diffs, 0), 4);
            grad_cost = sum(GreenPrior(diffs, 1), 4);
            hess_cost = sum(GreenPrior(diffs, 2), 4);            
        otherwise
            disp('Unknown approximation of L1 norm');
            return
    end

    % Multiply by the regularization weights
    value_cost = value_cost .* reshape(lambda, [1 1 nrealmat]);
    grad_cost = 2 * grad_cost .* reshape(lambda, [1 1 nrealmat]) / NbSubsets;
    hess_cost = 4 * hess_cost .* reshape(lambda, [1 1 nrealmat]) / NbSubsets;
    
    % Reshape for the output
    value_cost = sum(value_cost(:));
    grad_cost = reshape(grad_cost, [ObjectSize^2 nrealmat]);
    hess_cost = reshape(hess_cost, [ObjectSize^2 nrealmat]);
end