function [ out ] = GreenPrior( t, order )
%GREENPRIOR(t, order) is the order-th derivative of the prior defined in 
% "Bayesian Reconstructions From Emission Tomography Data Using a Modified
% EM Algorithm", by Peter J. Green, IEEE TMI 1990
%
%   t is the argument
%   order=0 means Huber function, order=1 is the first derivative, ...

c1 = 27/128;
c2 = 16/(3 * sqrt(3));
switch(order)
    case 0
        out = c1 * log(cosh(c2 * t));
    case 1
        out = c1 * c2 * tanh(c2 * t);
    case 2
        out = c1 * c2^2 ./ (cosh(c2 * t).^2);
    otherwise
        disp('Unhandled order')
        out = zeros(size(t), 'like', t);
        return
end

