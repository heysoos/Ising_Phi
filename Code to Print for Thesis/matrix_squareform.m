function [ sqr_mat ] = matrix_squareform( mat )
%Matrix Squareform: Apply squareform on the first 2 dimensions of a 3D
%matrix or the first dimension of a 2D matrix

if length(size(mat)) > 2
    for i = 1:size(mat,3)
        sqr_mat(:,i) = squareform(squeeze(mat(:,:,i)));
    end
else
    for i = 1:size(mat,2)
        sqr_mat(:,:,i) = squareform(mat(:,i));
    end
end

end