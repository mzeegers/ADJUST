function [A,F] = UR(Y,W,k)
%UR Unmixing-then-Reconstruction for spectral tomography
% 
% We solve the problem 
%                       W * A * F = Y 
% by two-step approach
%   (Unmixing)          P * F = Y
%   (Reconstruction)    W * A = P 
%
% Input:
%   Y - Spectral tomographic projections (size: m x c)
%   W - Tomographic projection operator (size: m x n)
%   k - Number of materials (note: k << min(m,c))
%
% Here, m is number of projections per spectral bin, c is the number of
% spectral bins, n is the size of image
%
% Output:
%   A - spatial distribution of materials (size: n x k)
%   F - spectral distribution of materials (size: k x c)
%
% Authors:
%   Ajinkya Kadu,
%       Centrum Wiskunde & Informatica, Amsterdam (aak@cwi.nl)
%   Mathé Zeegers, 
%       Centrum Wiskunde & Informatica, Amsterdam (M.T.Zeegers@cwi.nl)

% get sizes
[m,c] = size(Y);
[m,n] = size(W);

%% main loop

fprintf('==========================================================\n');
fprintf('     Unmixing - then - Reconstruction  \n');
fprintf('==========================================================\n');

%%% unmixing - Non-negative least squares (Mult followed by ALS)
opt     = statset('MaxIter',10,'Display','final');
[A0,F0] = nnmf(Y,k,'Options',opt,'Algorithm','mult','Replicates',10);

opt   = statset('MaxIter',1000,'Display','final');
[P,F] = nnmf(Y,k,'W0',A0,'H0',F0,'Options',opt,'Algorithm','als','Replicates',10);


%%% reconstruction with CG and Tikhonov regularization
regT = 1e-3;        % Tikhonov reg parameter 
Wtp  = W'*P;        % backprojection
WtW  = @(x) W'*(W*x) + regT*x;  % regularized forward operator

A    = zeros(n,k);
tol  = 1e-6;
iter = 20;
for i=1:k           % run for each material
   [A(:,i),~] = pcg(WtW,Wtp(:,i),tol,iter);
end

fprintf('==========================================================\n');

end
