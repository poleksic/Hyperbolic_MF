function [ P ] = h_Get_Nrlmf_Prob(U,V,theta)
    X = h_product(U,V) + theta; 
    P = exp(-X);
    P = 1 ./ (1 + P);
end
