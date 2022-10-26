function [sO] = computeResults(strProb, strRU, strUR, strJ, strD, sliceNo)
% computeResults Plot various comparison plots and tables for methods 
%   1) RU - Reconstruction-then-Unmixing
%   2) UR - Unmixing-the-Reconstruction
%   3) cJoint - classical joint reconstruction algorithm
%   4) ADJUST - proposed method 
%
% Input:
%   strProb - problem structure containing A (spatial maps), F (spectral 
%             maps) and n (size of spatial images) 
%   strRU - numerical results of RU method that contains A and F
%   strUR - numerical results of UR method that contains A and F
%   strJ  - numerical results of classical Joint method that contains A 
%           and F
%   strD  - numerical results of ADJUST method that contains A and F
%   sliceNo- slice number for visualizing results on 3d phantom
% 
% Output:
%   sO - output structure containing sO.UR, sO.RU, sO.cJoint, sO.ADJUST
%   that have the following details
%       A    : reconstructed spatial maps aligned with true maps 
%              (size: n x k)
%       F    : reconstructed spectral maps aligned with true maps
%              (size: k x c)
%       psnr : peak signal-to-noise ratio for each reconstructed material 
%              spatial map w.r.t. true map (size: kx1)
%       ssim : structural similarity index for each reconstructed material 
%              spatial map w.r.t. true map (size: kx1)
%       imse : mean-squared error for each reconstructed material 
%              spatial map w.r.t true map (size: kx1)
%       psnrM: mean PSNR over all materials
%       ssimM: mean SSIM over all materials
%       imseM: mean MSE over all materials
% 
%  
% Authors:
%   Ajinkya Kadu,
%       Centrum Wiskunde & Informatica, Amsterdam (aak@cwi.nl)
%   Mathé Zeegers, 
%       Centrum Wiskunde & Informatica, Amsterdam (M.T.Zeegers@cwi.nl)

%%

n   = strProb.n;
A   = strProb.A;
F   = strProb.F;

Aur = strUR.A;
Fur = strUR.F;
Aru = strRU.A;
Fru = strRU.F;
Aj  = strJ.A;
Fj  = strJ.F;
Ad  = strD.A;
Fd  = strD.F;

sO.strProb = strProb; 
sO.GT.A    = A;
sO.GT.F    = F;
%% UR

[Aurf, Furf] = alignMatrices(A, Aur, Fur);
sO.UR   = computePerformance(Aurf, A, n);
sO.UR.A = Aurf;
sO.UR.F = Furf;

%% RU

[Aruf, Fruf] = alignMatrices(A, Aru, Fru);
sO.RU   = computePerformance(Aruf, A, n);
sO.RU.A = Aruf;
sO.RU.F = Fruf;

%% cJoint

[Ajf, Fjf]   = alignMatrices(A, Aj, Fj);
sO.cJoint   = computePerformance(Ajf, A, n);
sO.cJoint.A = Ajf;
sO.cJoint.F = Fjf;

%% ADJUST

[Adf, Fdf]   = alignMatrices(A, Ad, Fd);
sO.ADJUST   = computePerformance(Adf, A, n);
sO.ADJUST.A = Adf;
sO.ADJUST.F = Fdf;

%% print

fprintf('\n--------------------------------------------------------------------------\n');

fprintf('%12s | %12s | %12s | %12s | %12s \n','','RU','UR',...
    'cJoint','ADJUST');

fprintf('--------------------------------------------------------------------------\n');

fprintf('%12s | %12.4f | %12.4f | %12.4f | %12.4f \n','PSNR', sO.RU.psnrM,...
    sO.UR.psnrM, sO.cJoint.psnrM, sO.ADJUST.psnrM);

fprintf('%12s | %12.4f | %12.4f | %12.4f | %12.4f \n','SSIM', sO.RU.ssimM,...
    sO.UR.ssimM, sO.cJoint.ssimM, sO.ADJUST.ssimM);

fprintf('%12s | %12.4f | %12.4f | %12.4f | %12.4f \n','MSE', sO.RU.imseM,...
    sO.UR.imseM, sO.cJoint.imseM, sO.ADJUST.imseM);

fprintf('--------------------------------------------------------------------------\n');


%% plot figure

k = size(A,2);

if length(n) == 2
    fig = figure(100); colormap hot;
    for i = 1:k
        AtI  = reshape(A(:,i), n);
        AruI = reshape(Aurf(:, i), n);
        AurI = reshape(Aurf(:, i), n);
        AjI  = reshape(Ajf(:, i), n);
        AdI  = reshape(Adf(:, i), n);

        subplot(k, 5, 5*(i-1)+1);
        imshow(AtI, [0 1]); title(sprintf('True-%d', i));
        subplot(k, 5, 5*(i-1)+2);
        imshow(AruI, [0 1]);title(sprintf('RU-%d', i));
        subplot(k, 5, 5*(i-1)+3);
        imshow(AurI, [0 1]); title(sprintf('UR-%d', i));
        subplot(k, 5, 5*(i-1)+4);
        imshow(AjI, [0 1]); title(sprintf('cJoint-%d', i));
        subplot(k, 5, 5*(i-1)+5);
        imshow(AdI, [0 1]); title(sprintf('ADJUST-%d', i));
        pause(0.1);
    end
    h = axes(fig,'visible','off');
    colorbar(h,'Position',[0.93 0.168 0.022 0.7]);

elseif length(n) == 3
    figure(100);
    for i = 1:k
        AtI  = reshape(A(:, i), n);
        AruI = reshape(Aurf(:, i), n);
        AurI = reshape(Aurf(:, i), n);
        AjI  = reshape(Ajf(:, i), n);
        AdI  = reshape(Adf(:, i), n);
        
        if nargin < 7, sliceNo = floor(n(3)/2); end

        subplot(k, 5, 5*(i-1)+1);
        imshow(AtI(:, :, sliceNo), []); title(sprintf('True-%d', i));
        subplot(k, 5, 5*(i-1)+2);
        imshow(AruI(:, :, sliceNo), []); title(sprintf('RU-%d', i));
        subplot(k, 5, 5*(i-1)+3);
        imshow(AurI(:, :, sliceNo), []); title(sprintf('UR-%d', i));
        subplot(k, 5, 5*(i-1)+4);
        imshow(AjI(:, :, sliceNo), []); title(sprintf('cJoint-%d', i));
        subplot(k, 5, 5*(i-1)+5);
        imshow(AdI(:, :, sliceNo), []); title(sprintf('ADJUST-%d', i));
        pause(0.1);
    end

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Supporting functions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [strOut] = computePerformance(A, At, n)
% computes the performance measures algorithm by comparing reconstructed 
% spatial maps with ground truth spatial maps
% We look at PSNR, SSIM and IMMSE
% Input - 
%   A - reconstructed spatial maps
%   At- ground truth spatial maps
%   n - size of images
%
% Output:
%   strOut - structure containing
%       psnr - peak signal-to-noise-ratio of individual spatial map
%       ssim - structural similarity index of individual spatial map
%       imse - mean-squared-error of individual spatial map
%       psnrM, ssimM, imseM - mean of these measures
%
k = size(At, 2);

psnrI = zeros(k, 1);
ssimI = zeros(k, 1);
imseI = zeros(k, 1);

% compute for each spatial map 
for i = 1:k
    Ix = reshape(A(:, i), n);
    Iy = reshape(At(:, i), n);
    Ix = double(Ix);
    Iy = double(Iy);
    psnrI(i, 1) = psnr(Ix, Iy);
    ssimI(i, 1) = ssim(Ix, Iy);
    imseI(i, 1) = immse(Ix, Iy);
end

strOut.psnr = psnrI;
strOut.ssim = ssimI;
strOut.imse = imseI;

% compute the averages
strOut.psnrM = mean(psnrI);
strOut.ssimM = mean(ssimI);
strOut.imseM = mean(imseI);

end


function [A2f, F2f, cumError, ErrorMatrix] = alignMatrices(A, A1, F1)
% align the reconstructed spatial and spectral maps based on the ground
% truth or given spatial and spectral maps
% 
% Input:
%   A - ground truth or given spatial maps 
%   A1 - spatial maps that needs to be aligned 
%   F1 - spectral maps that needs to be aligned
%
% Output:
%   A2f - aligned spatial maps
%   F2f - aligned spectral maps
%   cumError - cumulative errors
%   ErrorMatrix - mean-squared error for each pair (spatial maps only).

Amax = max(A1(:));
A2 = A1 / Amax;
F2 = F1 * Amax;

[cumError, ErrorMatrix, minrowI, mincolI] = computeError(A2, A);

[~, idD]= sort(mincolI);
A2f = A2(:, minrowI(idD));
F2f = F2(minrowI(idD), :);

end


function [cumError, ErrorMatrix, minrowI, mincolI] = computeError(A, Atrue)
% Function for calculating the error measure per image
%
% Input:
%  A -  reconstructed spatial maps
%  Atrue - ground truth spatial maps
%
% Output:
%   cumError - cumulative error for each material
%   ErrorMatrix - mean-squared-error for each material spatial maps pair
%   minrowI - 
%   mincolI - 
% 

nmat = size(Atrue, 2);
nPix = size(Atrue, 1);

% First compute the errors for all image pairs
ErrorMatrix = zeros(nmat, nmat);
for i = 1:nmat
    for j = 1:nmat
        ErrorMatrix(i, j) = norm(A(:, i) - Atrue(:, j));
    end
end

% ErrorMatrix

% Now do the matching and iteratively removing rows and colums
cumError = 0;
minrowI = zeros(nmat, 1);
mincolI = zeros(nmat, 1);

% psnrv   = zeros(nmat,1);
% ssimv   = zeros(nmat,1);
% immsev  = zeros(nmat,1);

for k = 1:nmat

    %Find minimum value and indices in error array
    [minvalue, minindex] = min(ErrorMatrix(:));
    [minrow, mincol] = ind2sub(size(ErrorMatrix), minindex);

    % trueI = Atrue(:,mincol);
    % recI  = A(:,minrow);

    errorPair = ErrorMatrix(minrow, mincol);
    cumError  = cumError + errorPair;

    % fprintf("Matched image %d to true image %d with error %f \n", minrow, mincol, errorPair)

    minrowI(k, 1) = minrow;
    mincolI(k, 1) = mincol;

    %Remove the corresponding row and column from the array

    ErrorMatrix(minrow, :) = Inf;
    ErrorMatrix(:, mincol) = Inf;

end

% fprintf("The total error is %f \n", cumError)

end

