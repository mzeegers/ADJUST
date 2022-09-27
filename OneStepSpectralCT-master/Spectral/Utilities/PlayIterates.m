function [] = PlayIterates(iterates,mat_inf, mat_sup)
[ObjectSize, nmat, niterates] = size(iterates);
ObjectSize = sqrt(ObjectSize);

% Normalize each material to uint8
mat_rescaled = uint8(zeros(size(iterates)));
for mat=1:nmat
    mat_rescaled(:,mat,:) = uint8((iterates(:,mat,:) - mat_inf(mat)) * 255 / (mat_sup(mat) - mat_inf(mat)));
end

% Play sequences side by side
implay(reshape(mat_rescaled, [ObjectSize, ObjectSize * nmat, niterates]));

% % Write it
% writeVideoResult(double(reshape(mat_rescaled, [ObjectSize, ObjectSize * nmat, niterates])), 'WeidingerNoRegul256_2000iter.avi', 0,255);
end