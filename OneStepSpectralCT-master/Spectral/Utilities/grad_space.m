function [ out ] = grad_space( input )
%GRAD_SPACE discrete uncentered gradient with boundary conditions as
%described in "Chambolle 2004"

    [h, w, nchannels]=size(input);
    out = zeros(h, w, nchannels, 2, 'like', input);

    pad_in = zeros(h+1,w+1,nchannels, 'like', input);
    pad_in(1:end-1,1:end-1,:) = input;

    pgx = diff(pad_in, 1, 1); pgx(end,:,:,:)=0;
    pgy = diff(pad_in, 1, 2); pgy(:,end,:,:)=0;
    out(:, :, :, 1) = pgx(:,1:end-1,:);
    out(:, :, :, 2) = pgy(1:end-1,:,:);
end