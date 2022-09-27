function [out, passage] = mu_preconditioning(in, syntheticMaterials, spectrum)

    if (strcmp(syntheticMaterials,'normalize'))
        disp('Using normalized attenuations')

        % Very crude way of normalizing the attenuation matrix,
        % while keeping track of the transformation (to invert it later)
        passage = eye(3);
        out = in;
        temp_passage_mat = eye(3); temp_passage_mat(1,1)= 1/norm(out(:,1));
        out = out * temp_passage_mat;
        passage = passage * temp_passage_mat;

        temp_passage_mat = eye(3); temp_passage_mat(2,2)= 1/norm(out(:,2));
        out = out * temp_passage_mat;
        passage = passage * temp_passage_mat;

        temp_passage_mat = eye(3); temp_passage_mat(3,3)= 1/norm(out(:,3));
        out = out * temp_passage_mat;
        passage = passage * temp_passage_mat;
    end
    if (strcmp(syntheticMaterials,'orthonormalize'))
        disp('Using orthonormalized attenuations')
        
        % Very crude way of orthonormalizing the attenuation matrix,
        % while keeping track of the transformation (to invert it later)
        passage = eye(3);
        out = in;
        temp_passage_mat = eye(3); temp_passage_mat(1,1)= 1/norm(out(:,1));
        out = out * temp_passage_mat;
        passage = passage * temp_passage_mat;

        dotp = (out(:,2).' * out(:,1));
        temp_passage_mat = eye(3); temp_passage_mat(1,2)=-dotp;
        out = out * temp_passage_mat;
        passage = passage * temp_passage_mat;
        temp_passage_mat = eye(3); temp_passage_mat(2,2)= 1/norm(out(:,2));
        out = out * temp_passage_mat;
        passage = passage * temp_passage_mat;

        dotp1 = (out(:,3).' * out(:,1));
        dotp2 = (out(:,3).' * out(:,2));
        temp_passage_mat = eye(3); temp_passage_mat(1,3)=-dotp1;temp_passage_mat(2,3)=-dotp2;
        out = out * temp_passage_mat;
        passage = passage * temp_passage_mat;
        temp_passage_mat = eye(3); temp_passage_mat(3,3)= 1/norm(out(:,3));
        out = out * temp_passage_mat;
        passage = passage * temp_passage_mat;
    end
    if (strcmp(syntheticMaterials,'fessler'))
        disp('Using Fessler trick')
        
        if (numel(size(spectrum))==3) % If spectrum is pixel-dependent
            S = spectrum(:,:,floor(size(spectrum, 3)/2)); % Take the spectrum value in the middle
        else
            S = spectrum;
        end

        % Construct a basis of artificial materials, proposed by Fessler in a 1999 patent
        inv_passage = (S * in )./(S * ones(size(in))); % the beta_bar matrix in equation 71 of Fessler's patent
        passage = inv(inv_passage.' * inv_passage) * inv_passage.'; % The Moore-Penrose pseudo inverse, since the inv_passage matrix is not square
        out = in * passage;
    end
    if (strcmp(syntheticMaterials,'barber'))
        disp('Using Barber trick')

        % Construct a basis of artificial materials, proposed by Barber in her 2016 paper
        SymReal = in.' * in;
        [V,D] = eig(SymReal);
        P_inv = V * sqrt(inv(D));
        out = in * P_inv;
        passage = P_inv;
    end

end