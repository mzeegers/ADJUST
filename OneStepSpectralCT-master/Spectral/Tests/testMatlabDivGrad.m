ObjectSize = 512;

% Generate u and compute its gradient
u = rand(ObjectSize);
pu = padarray(u, [1 1], 0, 'both'); 
[pgx_u, pgy_u] = gradient(pu);
gx_u =  pgx_u(2:end-1, 2:end-1);
gy_u =  pgy_u(2:end-1, 2:end-1);

% Generate the gradient of v and compute its divergence
gx_v = rand(ObjectSize); gy_v = rand(ObjectSize);
pgx_v = padarray(gx_v, [1 1], 0, 'both'); 
pgy_v = padarray(gy_v, [1 1], 0, 'both'); 
pdiv_v = divergence(pgx_v, pgy_v);
div_v = pdiv_v(2:end-1, 2:end-1);

% Compute the dot products
dotprod1 = dot(gx_u(:), gx_v(:)) + dot(gy_u(:), gy_v(:));
dotprod2 = dot(u(:),-div_v(:));

% Compute non-adjointness index
non_adjointness = 1 - dotprod1 / dotprod2



% Compare gradient by matlab and by convolution
lena = imread('/home/mory/Pictures/lena_std.tif');
lena = sum(double(lena), 3);
lena = lena(150:200,150:200);
lena = padarray(lena, [1 1], 0, 'both'); 

% Compute the gradient by Matlab
[gx_lena, gy_lena] = gradient(lena);
gx_lena = gx_lena(2:end-1, 2:end-1);
gy_lena = gy_lena(2:end-1, 2:end-1);

% Compute the gradient by convolution with a kernel
gx_lena_by_conv = conv2(lena, [0.5, 0, -0.5], 'same');
gy_lena_by_conv = conv2(lena, [0.5; 0; -0.5], 'same');
gx_lena_by_conv = gx_lena_by_conv(2:end-1, 2:end-1);
gy_lena_by_conv = gy_lena_by_conv(2:end-1, 2:end-1);

figure(1); Show(gx_lena);
figure(2); Show(gx_lena_by_conv);
figure(3); Show(gx_lena - gx_lena_by_conv);
figure(4); Show(lena)


% Compare divergence by Matlab and by convolution
pgx_lena = padarray(gx_lena, [1 1], 0, 'both');
pgy_lena = padarray(gy_lena, [1 1], 0, 'both');

% Compute the divergence by Matlab
pdiv_lena = divergence(pgx_lena, pgy_lena);
div_lena = pdiv_lena(2:end-1,2:end-1);

% Compute the divergence by convolution
pdiv_lena_by_conv = conv2(pgx_lena, [0.5, 0, -0.5], 'same') + conv2(pgy_lena, [0.5; 0; -0.5], 'same');
div_lena_by_conv = pdiv_lena_by_conv(2:end-1,2:end-1);

figure(1); Show(div_lena);
figure(2); Show(div_lena_by_conv);










