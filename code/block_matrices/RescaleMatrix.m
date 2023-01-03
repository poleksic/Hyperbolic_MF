function K = RescaleMatrix(B, largest)
    % Rescale B so that scores are between 0 and 1
    K = largest * B ./ max(B(:));
    K(isnan(K)) = 0;
end