function [ P ] = h_Get_Weighted_Nrlmf_Prob(U,WU,V,WV,theta)
    W = WU * WV';
    X = h_product(U,V) + theta*W; 
    P = exp(-X);
    P = 1 ./ (1 + P);
end
