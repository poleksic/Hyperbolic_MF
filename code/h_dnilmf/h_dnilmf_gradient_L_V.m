function GRAD = h_dnilmf_gradient_L_V(R, U, V, M, N, c, alphaU, alpha, beta, gamma, theta)
    
    % Input: 
    % W - weight matrix, R - relationship matrix, P - matrix of probabilities
    % U, V - matrices of latent vectors, alpha - trainable parameter
    
    % Output: 
    % Matrix of gradient vectors (rows) of the loss function
    % Note: This is the true gradient in Minkowski space (as the last column is already multiplied by -1)
    small_positive = 1e-12;
    
    [n, r] = size(V);    
    In = eye(n);
    d = r - 1;  % dimension of hyperbolic space
    
    arccosh = zeros(n,1);
    der_arccosh = zeros(n,1);
    
    for j=1:n
        if abs(V(j,r) - 1) < small_positive
            arccosh(j) = 0;
            der_arccosh(j) = 0;
        else
            arccosh(j) = acosh(V(j,r));
            der_arccosh(j) = sqrt(V(j,r)^2 -1);
        end
    end
    
    P = h_Get_Dnilmf_Prob(U,V,M,N,alpha,beta,gamma,theta);
    
    Q = (1 + c*R - R) .* P;
    
    Rt = R';
    Nt = N';
    Qt = Q';
    
    GRAD = (alpha * In + gamma * Nt) * Qt * U - c * (alpha * In + gamma * N) * Rt * U - beta * (c * Rt - Qt) * M * U;
    
    for j=1:n
        if abs(V(j,r) - 1) < small_positive
            GRAD(j,r) = GRAD(j,r) - (d-1)/3;
        else
            GRAD(j,r) = GRAD(j,r) - 2*alphaU*arccosh(j)/der_arccosh(j) - (d-1)* (V(j,r)*arccosh(j)-der_arccosh(j))/(arccosh(j)*der_arccosh(j)^2);
        end
    end
    
    % linear
%    GRAD = GRAD - 4*beta*(H*V);

    % Euclid regularization (discard)
%      GRAD = GRAD + 4*beta*(V .* (H*ones(n,r)) - H*V);

    % distance regularization (use this)
    % adding to the loss function the sum of the terms s_ij*d^2(Ui,Uj)
    % where d is the distance (not Lorentzian but ordinary) between the 
    % points in the hyperbolic space (d(x,y) = arccosh(-<x,y>))
%     HP = h_product(V,V) - theta;
%     B = H.*acosh(-1*HP)./ sqrt(HP.*HP -ones(n,n));
     % B(i,i) = 0/0 but is equal to t_ii in the limit
%     B(1:n+1:end) = diag(H);    
%     GRAD = GRAD - 4*beta*(B*V);
end
