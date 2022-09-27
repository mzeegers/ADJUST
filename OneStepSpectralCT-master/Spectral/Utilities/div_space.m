function [ out ] = div_space( grad_input )
%DIV 2D discrete divergence with correct boundary conditions for TV
%minimization

    [h, w, nchannels, ~] = size(grad_input);       
    pinx = zeros(h+1, w, nchannels, 'like', grad_input);
    pinx(2:end-1, :, :) = grad_input(1:end-1,:,:,1);
    piny = zeros(h, w+1, nchannels, 'like', grad_input);
    piny(:, 2:end-1, :) = grad_input(:,1:end-1,:,2);
    out = diff(pinx, 1, 1) + diff(piny, 1, 2);
end