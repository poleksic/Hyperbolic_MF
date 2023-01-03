function [ P ] = h_ngn_getProb(U_t,U_d,V_t,V_d,theta)
    X = 0.5* (h_product(U_d,V_t) + h_product(V_d,U_t)) + theta; 
    P = exp(-X);
    P = 1 ./ (1 + P);
end
