function [] = PlayIteratesTransposed(iterates,mat_inf, mat_sup)
[ObjectSize, nmat, niterates] = size(iterates);
ObjectSize = sqrt(ObjectSize);

% Normalize each material to uint8
tmp = iterates;
tmp(:,1,:) = tmp(:,1,:) * 4.933;
tmp(:,2,:) = tmp(:,2,:) * 7.9;

mat_rescaled = zeros(size(iterates));
for mat=1:nmat
    mat_rescaled(:,mat,:) = (tmp(:,mat,:) - mat_inf(mat)) / (mat_sup(mat) - mat_inf(mat));
end

% Play sequences side by side
mat_rescaled=permute(reshape(mat_rescaled, [ObjectSize, ObjectSize, nmat, niterates]), [2 1 3 4]);

implay(reshape(mat_rescaled, [ObjectSize, ObjectSize*nmat, niterates]));

% % Write it
% writeVideoResult(double(reshape(mat_rescaled, [ObjectSize, ObjectSize * nmat, niterates])), 'WeidingerNoRegul256_2000iter.avi', 0,255);
end