function [ P ] = h_getSquaredLoentzianProb(U,V,theta)
    X = -1 * h_squared_Lorentzian_distance(U,V) + theta;

    P = exp(X);
    P = P ./ (1 + P);
end
