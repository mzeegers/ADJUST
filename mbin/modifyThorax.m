function [Ao,Fo] = modifyThorax(A,F)
%modifyThorax function only required to convert the Thorax phantom for
% computational purposes. Thorax phantoms consists of 5 materials, of which
% 2 are hard materials and 3 soft materials. For our spectral
% reconstruction techniques, we would like to combine 3 soft materials into
% one. 
%   
% Input:
%    A - The thorax phantom with 5 materials
%    F - The spectra fo the 5 materials
% 
% Output:
%    Ao - The thorax phantom with 3 (combined) materials
%    Fo - The spectra for the 3 (combined) materials
% 
%
% Authors:
%   Ajinkya Kadu,
%       Centrum Wiskunde & Informatica, Amsterdam (aak@cwi.nl)
%   Mathé Zeegers, 
%       Centrum Wiskunde & Informatica, Amsterdam (M.T.Zeegers@cwi.nl)
                        

[n,k] = size(A);
[k,c] = size(F);

Ao = zeros(n,k-2);
Fo = zeros(k-2,c);

% hard materials
for i=1:k-3
    Ao(:,i) = A(:,i);
    Fo(i,:) = F(i,:);
end

% combine soft materials
Asoft = zeros(n,1);
Fsoft = zeros(1,c);
for i=k-2:k
    Asoft = Asoft + A(:,i);
    Fsoft = Fsoft + F(i,:);
end

Ao(:,k-2) = Asoft;
Fo(k-2,:) = Fsoft;


end
