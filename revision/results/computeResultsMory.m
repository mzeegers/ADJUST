function [sO] = computeResults(strProb,strRU,strUR,strJ,strD,strCai,strWeidinger,strLong,strMechlem,strBarber,sliceNo)
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
ACai  = strCai.A;
FCai  = strCai.F;
AWeidinger  = strWeidinger.A;
FWeidinger  = strWeidinger.F;
ALong  = strLong.A;
FLong  = strLong.F;
AMechlem  = strMechlem.A;
FMechlem  = strMechlem.F;
ABarber  = strBarber.A;
FBarber  = strBarber.F;

sO.GT.A = A;
sO.GT.F = F;


%% UR

[Aurf,Furf] = alignMatrices(A,Aur,Fur,n);
sO.UR       = computePerformance(Aurf,A,n);
sO.UR.A     = Aurf;
sO.UR.F     = Furf;

%% RU

[Aruf,Fruf] = alignMatrices(A,Aru,Fru,n);
sO.RU       = computePerformance(Aruf,A,n);
sO.RU.A     = Aruf;
sO.RU.F     = Fruf;

%% cJoint

[Ajf,Fjf]   = alignMatrices(A,Aj,Fj,n);
sO.cJoint   = computePerformance(Ajf,A,n);
sO.cJoint.A = Ajf;
sO.cJoint.F = Fjf;

%% ADJUST

[Adf,Fdf]   = alignMatrices(A,Ad,Fd,n);
sO.ADJUST   = computePerformance(Adf,A,n);
sO.ADJUST.A = Adf;
sO.ADJUST.F = Fdf;

%% Cai

[ACaif,FCaif]   = alignMatrices(A,ACai,FCai,n);
sO.Cai   = computePerformance(ACaif,A,n);
sO.Cai.A = ACaif;
sO.Cai.F = FCaif;

%% Weidinger

[AWeidingerf,FWeidingerf]   = alignMatrices(A,AWeidinger,FWeidinger,n);
sO.Weidinger   = computePerformance(AWeidingerf,A,n);
sO.Weidinger.A = AWeidingerf;
sO.Weidinger.F = FWeidingerf;

%% Long

[ALongf,FLongf]   = alignMatrices(A,ALong,FLong,n);
sO.Long   = computePerformance(ALongf,A,n);
sO.Long.A = ALongf;
sO.Long.F = FLongf;

%% Mechlem

[AMechlemf,FMechlemf]   = alignMatrices(A,AMechlem,FMechlem,n);
sO.Mechlem   = computePerformance(AMechlemf,A,n);
sO.Mechlem.A = AMechlemf;
sO.Mechlem.F = FMechlemf;

%% Barber

[ABarberf,FBarberf]   = alignMatrices(A,ABarber,FBarber,n);
sO.Barber   = computePerformance(ABarberf,A,n);
sO.Barber.A = ABarberf;
sO.Barber.F = FBarberf;


%% print

fprintf('\n--------------------------------------------------------------------------\n');

fprintf('%12s | %12s | %12s | %12s | %12s | %12s | %12s | %12s | %12s | %12s \n','','RU','UR',...
    'cJoint','ADJUST', 'Cai', 'Weidinger', 'Long', 'Mechlem', 'Barber');

fprintf('--------------------------------------------------------------------------\n');

fprintf('%12s | %12.4f | %12.4f | %12.4f | %12.4f | %12.4f | %12.4f | %12.4f | %12.4f | %12.4f \n','PSNR',sO.RU.psnrM,...
    sO.UR.psnrM,sO.cJoint.psnrM,sO.ADJUST.psnrM, sO.Cai.psnrM, sO.Weidinger.psnrM, sO.Long.psnrM, sO.Mechlem.psnrM, sO.Barber.psnrM);

fprintf('%12s | %12.4f | %12.4f | %12.4f | %12.4f | %12.4f | %12.4f | %12.4f | %12.4f | %12.4f \n','SSIM',sO.RU.ssimM,...
    sO.UR.ssimM,sO.cJoint.ssimM,sO.ADJUST.ssimM, sO.Cai.ssimM, sO.Weidinger.ssimM, sO.Long.ssimM, sO.Mechlem.ssimM, sO.Barber.ssimM);

fprintf('%12s | %12.4f | %12.4f | %12.4f | %12.4f | %12.4f | %12.4f | %12.4f | %12.4f | %12.4f \n','MSE',sO.RU.imseM,...
    sO.UR.imseM,sO.cJoint.imseM,sO.ADJUST.imseM, sO.Cai.imseM, sO.Weidinger.imseM, sO.Long.imseM, sO.Mechlem.imseM, sO.Barber.imseM);

fprintf('--------------------------------------------------------------------------\n');


%% plot figure

k = size(A,2);

if length(n) == 2
    fig = figure(100);colormap hot;
    for i=1:k
        AtI  = reshape(A(:,i),n);
        AruI = reshape(Aurf(:,i),n);
        AurI = reshape(Aurf(:,i),n);
        AjI  = reshape(Ajf(:,i),n);
        AdI  = reshape(Adf(:,i),n);

        subplot(k,5,5*(i-1)+1);imagesc(AtI,[0 1]);title(sprintf('True-%d',i));
        axis image;set(gca,'XTick',[],'YTick',[]);
        subplot(k,5,5*(i-1)+2);imagesc(AruI,[0 1]);title(sprintf('RU-%d',i));
        axis image;set(gca,'XTick',[],'YTick',[]);
        subplot(k,5,5*(i-1)+3);imagesc(AurI,[0 1]);title(sprintf('UR-%d',i));
        axis image;set(gca,'XTick',[],'YTick',[]);
        subplot(k,5,5*(i-1)+4);imagesc(AjI,[0 1]);title(sprintf('cJoint-%d',i));
        axis image;set(gca,'XTick',[],'YTick',[]);
        subplot(k,5,5*(i-1)+5);imagesc(AdI,[0 1]);title(sprintf('ADJUST-%d',i));
        axis image;set(gca,'XTick',[],'YTick',[]);
        pause(1);
    end
    h = axes(fig,'visible','off');
    colorbar(h,'Position',[0.93 0.168 0.022 0.7]);

elseif length(n) == 3
    figure(100);
    for i=1:k
        AtI  = reshape(A(:,i),n);
        AruI = reshape(Aurf(:,i),n);
        AurI = reshape(Aurf(:,i),n);
        AjI  = reshape(Ajf(:,i),n);
        AdI  = reshape(Adf(:,i),n);
        
        if nargin < 7 
            sliceNo = floor(n(3)/2);
        end

        subplot(k,5,5*(i-1)+1);imshow(AtI(:,:,sliceNo));title(sprintf('True-%d',i));
        subplot(k,5,5*(i-1)+2);imshow(AruI(:,:,sliceNo));title(sprintf('RU-%d',i));
        subplot(k,5,5*(i-1)+3);imshow(AurI(:,:,sliceNo));title(sprintf('UR-%d',i));
        subplot(k,5,5*(i-1)+4);imshow(AjI(:,:,sliceNo));title(sprintf('cJoint-%d',i));
        subplot(k,5,5*(i-1)+5);imshow(AdI(:,:,sliceNo));title(sprintf('ADJUST-%d',i));
        pause(1);
    end

end

end

function [strOut] = computePerformance(A,At,n)

k = size(At,2);

psnrI = zeros(k,1);
ssimI = zeros(k,1);
imseI = zeros(k,1);


for i=1:k
    Ix = reshape(A(:,i),n);
    Iy = reshape(At(:,i),n);
    Ix = double(Ix);
    Iy = double(Iy);
    psnrI(i,1) = psnr(Ix,Iy);
    ssimI(i,1) = ssim(Ix,Iy);
    imseI(i,1) = immse(Ix,Iy);
end

strOut.psnr = psnrI;
strOut.ssim = ssimI;
strOut.imse = imseI;

strOut.psnrM = mean(psnrI);
strOut.ssimM = mean(ssimI);
strOut.imseM = mean(imseI);

end


function [A2f,F2f,cumError,ErrorMatrix] = alignMatrices(A,A1,F1,n)

    Amax = max(A1(:));
    A2 = A1/Amax;
    F2 = F1*Amax;

    [cumError,ErrorMatrix,minrowI,mincolI] = computeError(A2,A);

    [~,idD]= sort(mincolI);
    A2f = A2(:, minrowI(idD));
    F2f = F2(minrowI(idD),:);

end


%%% Function for calculating the error measure per image
function [cumError,ErrorMatrix,minrowI,mincolI] = computeError(A,Atrue)

    B = A;
    Btrue = Atrue;
    
    nmat = size(Atrue,2);
    nPix = size(Atrue,1);
    
    %First compute the errors for all image pairs
    ErrorMatrix = zeros(nmat, nmat);
    for i = 1:nmat
        for j = 1:nmat
            ErrorMatrix(i,j) = norm(B(:,i) - Btrue(:,j));
        end
    end
    
    % ErrorMatrix
    
    %Now do the matching and iteratively removing rows and colums
    cumError = 0;
    minrowI = zeros(nmat,1);
    mincolI = zeros(nmat,1);
    psnrv   = zeros(nmat,1);
    ssimv   = zeros(nmat,1);
    immsev  = zeros(nmat,1);
    
    for k = 1:nmat

        %Find minimum value and indices in error array
        [minvalue, minindex] = min(ErrorMatrix(:));
        [minrow, mincol] = ind2sub(size(ErrorMatrix), minindex);
        
        trueI = Btrue(:,mincol);
        recI  = B(:,minrow);
        
        errorPair = ErrorMatrix(minrow,mincol);
        cumError = cumError + errorPair;
        
        % fprintf("Matched image %d to true image %d with error %f \n", minrow, mincol, errorPair)
        
        minrowI(k,1) = minrow;
        mincolI(k,1) = mincol;
        
        %Remove the corresponding row and column from the array
        
        ErrorMatrix(minrow,:) = Inf;
        ErrorMatrix(:,mincol) = Inf;
    
    end
    
    % fprintf("The total error is %f \n", cumError) 
        
end

function [Af] = buildSingleImageforA(A,idx)

[m,n] = size(A);

if nargin < 2
    idx = linspace(1,n,n);
end

vdx = linspace(1,n,n);

Af = A(:,idx)*vdx';


end

