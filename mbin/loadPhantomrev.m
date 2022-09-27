function [A,F] = loadPhantomrev(phantomId,n,T,F)

if phantomId < 1 || phantomId > 7
    error('phantomId must be one of 1-7\n');
end


switch phantomId
    case 1
        A = spectralSheppLogan(n);
        F = T(0+(1:size(A,2)),:);
    case 2
        k = 8;
        Disk = false;
        A = spectralPhantomCircles(2*n,k,Disk);
        F = T(10+(1:size(A,2)),:);
    case 3
        k = 5;
        Disk = false;   %return argument to True
        A = spectralPhantomCircles(2*n,k,Disk);
        F = T(18+(1:size(A,2)),:);   %return 0 to 18
    case 4
        A = spectralSheppLogan3D(n);
        F = T(0+(1:size(A,2)),:);
    case 5
        A0 = spectralThoraxPhantom();    %Need to be 512
        A  = [A0(:,1)+A0(:,2) A0(:,3)+A0(:,4) A0(:,5)+A0(:,7) A0(:,6) A0(:,8)];
        F  = F(1:5,:); % T(0+(1:size(A,2)),:);   %Needs to be changed/discussed
    case 6
        k = 5;
        A = spectralPhantomCirclesMixed(2*n, k, false);
        F = T(0+(1:size(A,2)),:);
    case 7
        k = 11;
        A = spectralPhantomCircles(2*n, k, false);
        F = T(0+(1:size(A,2)),:);
end

end
