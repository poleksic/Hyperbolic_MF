function [ P ] = GetP(X)
    P = exp(-X);
    P = 1 ./ (1 + P);
end