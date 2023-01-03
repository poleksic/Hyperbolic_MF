function K = ExtractAboveZscoreCutoff(B, zcut)
    % Keep only significantly large (with score > zcut) entries of B
    % B is assumed to be a similarity matrix (diagonal = 1)

    [m,n] = size(B);
    if m ~= n
        error('ExtractAboveZscoreCutoff: not a square matrix\n');
    end
    B(1:m+1:m*m) = realmax;
    T = B(:);
    T( T>= realmax) = [];
    mn = mean(T);
    B(1:m+1:m*m) = mn;
    Z = zscore(B(:));
    A = reshape(Z,[m,m]);
    K = B;
    K(A < zcut) = 0;
    K(1:m+1:m*m) = 1;
end