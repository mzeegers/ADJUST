function y = hyperbola(x, delta, order)
%HYPERBOLA Hyperbola(t, delta, order) is the order-th derivative of the 
% hyperbola function proposed in Long & Fessler 2014
%   t is the argument
%   delta is the parameter (cutoff)
%   order=0 means hyperbola function, order=1 is the first derivative, ...

% Prepare delta for implicit extension
[~, ~, channels, ~] = size(x);
if not((channels == numel(delta)))
    disp('Inconsistent delta and input');
    return
end

d = reshape(delta, [1 1 channels 1]); % Prepare for implicit extension

B = 1 + 3 * (x ./ d).^2;
switch(order)
    case 0
        y = (d.^2)/3 .* (sqrt(B) - 1);
    case 1
        y = B.^(-1/2) .* x;
    case 2
        y = B.^(-3/2);
    otherwise
        disp('Unhandled order')
        y = zeros(size(x), 'like', x);
        return
end