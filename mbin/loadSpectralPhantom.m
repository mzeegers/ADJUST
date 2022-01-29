function [A,F,T,q] = loadSpectralPhantom(phantomId,n)
% loadSpectralPhantom gives the spatial and material maps for various
% phantoms. It also gives the spectral dictionary that may be required for
% ADJUST algorithm. 

if nargin < 2
    if(phantomId == 'Thorax')
        n = 512;
    else
        n = 128;
    end
end

switch phantomId
    case 'SheppLogan'
        A = spectralSheppLogan(n);
        T = readmatrix('MatrixTXrayHard100chan_93mat_20MinEn_119MaxEn.csv');
        F = T(0+(1:size(A,2)),:);
        q = readmatrix('VectorIXrayMatMoEn100Min5Max35.csv');
    case 'Disks'
        k = 8;
        A = spectralPhantomCircles(n,k);
        T = readmatrix('MatrixTXrayHard100chan_93mat_20MinEn_119MaxEn.csv');
        F = T(10+(1:size(A,2)),:);
        q = readmatrix('VectorIXrayMatMoEn100Min5Max35.csv');
    case 'Thorax'
        A0 = spectralThoraxPhantom();   
        A  = [A0(:,1)+A0(:,2) A0(:,3)+A0(:,4) A0(:,5)+A0(:,7) A0(:,6) A0(:,8)];
        T  = readmatrix('Thorax_MatrixTXray100chan_60mat_MinEn20_MaxEn80.csv');
        F  = readmatrix('Thorax_MatrixFXray100chan_6matBoneBloodIodineSofttissueBloodLungsIodine_MinEn20_MaxEn80.csv');
        F  = F(1:5,:); %remove the pure iodine from the dictionary 
        q  = readmatrix('VectorIXrayMatWEn100Min20Max80.csv');
    case 'MixedDisk'
        k = 5;
        A = spectralPhantomCirclesMixed(n, k);
        T = readmatrix('MatrixTXrayHard100chan_93mat_20MinEn_119MaxEn.csv');
        F = T(0+(1:size(A,2)),:);
        q = readmatrix('VectorIXrayMatMoEn100Min5Max35.csv');
    case 'SheppLogan3D'
        A = spectralSheppLogan3D(n);
        T = readmatrix('MatrixTXrayHard100chan_93mat_20MinEn_119MaxEn.csv');
        F = T(0+(1:size(A,2)),:);
        q = readmatrix('VectorIXrayMatMoEn100Min5Max35.csv');   
    otherwise
        error('phantom can either be SheppLogan, Disks, Thorax or Mixed \n');
end

end
