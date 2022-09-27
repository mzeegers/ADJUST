function [ out ] = Huber( t, delta, order )
%HUBER Huber(t, delta, order) is the order-th derivative of the Huber
%function
%   t is the argument
%   delta is the parameter (cutoff)
%   order=0 means Huber function, order=1 is the first derivative, ...

% Prepare delta for implicit extension
[~, ~, channels, ~] = size(t);
if not((channels == numel(delta)))
    disp('Inconsistent delta and input');
    return
end

d = reshape(delta, [1 1 channels 1]); % Prepare for implicit extension

switch(order)
    case 0
        out = 2 * d .* abs(t) - d.^2;
        out(abs(t)<d) = t(abs(t)<d).^2;
    case 1
        out = 2 * d .* sign(t);
        out(abs(t)<d) = 2 * t(abs(t)<d);
    case 2
        out = zeros(size(t), 'like', t);
        out(abs(t)<d) = 2;
    otherwise
        disp('Unhandled order')
        out = zeros(size(t), 'like', t);
        return
end

