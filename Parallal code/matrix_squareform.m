function [ sqr_mat ] = matrix_squareform( mat )
%Matrix Squareform: Apply squareform on the first 2 dimensions of a 3D
%matrix or the first dimension of a 2D matrix

if ndims(mat) == 3 % for (2+1)D matrices
    for i = 1:size(mat,3)
        
        x = squeeze(mat(:,:,i));
        x(1:length(x)+1:end) = 0;
        
        sqr_mat(:,i) = squareform(x);
    end
elseif ndims(mat) ==  2 % for (1+1)D vectors
    for i = 1:size(mat,2)
        sqr_mat(:,:,i) = squareform(mat(:,i));
    end
else
   msg = 'Matrix must be 2 or 3 dimensional!';
   error(msg)
end

end                                                    