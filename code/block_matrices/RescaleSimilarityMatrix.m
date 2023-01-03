function K = RescaleSimilarityMatrix(B, largest)
    [m,~] = size(B);
    K = B;
    K(1:m+1:m*m) = 0;
    K = RescaleMatrix(K,largest);
    K(1:m+1:m*m) = 1;
end