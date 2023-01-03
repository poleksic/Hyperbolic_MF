function [ P ] = h_Get_Dnilmf_Prob(U, V, M, N, alpha, beta, gamma, theta)
    X = alpha * h_product(U,V) + beta * M * h_product(U,V) + gamma * h_product(U,V) * N' + theta;
    P = exp(-X);
    P = 1 ./ (1 + P);
end