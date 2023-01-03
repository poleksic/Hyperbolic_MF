function GRAD = h_dnilmf_gradient_L_U(R, U, V, M, N, c, alphaU, alpha, beta, gamma, theta)
    
    % Input: 
    % W - weight matrix, R - relationship matrix, P - matrix of probabilities
    % U, V - matrices of latent vectors, alpha - trainable parameter
    
    % Output: 
    % Matrix of gradient vectors (rows) of the loss function
    % Note: This is the true gradient in Minkowski space (as the last column is already multiplied by -1)
    
    small_positive = 1e-12;
    
    [m, ~] = size(R);
    Im = eye(m);
    [m, r] = size(U);
    d = r - 1;  % dimension of hyperbolic space
    
    arccosh = zeros(m,1);
    der_arccosh = zeros(m,1);
    
    for i=1:m
        if abs(U(i,r) - 1) < small_positive
            arccosh(i) = 0;
            der_arccosh(i) = 0;
        else
            arccosh(i) = acosh(U(i,r));
            der_arccosh(i) = sqrt(U(i,r)^2 -1);
        end
    end
    Mt = M';
    Nt = N';
    
    P = h_Get_Dnilmf_Prob(U,V,M,N,alpha,beta,gamma,theta);
    Q = (1 + c*R - R) .* P;
    
    GRAD = (alpha * Im + beta * Mt) * Q * V - c * (alpha * Im + beta * Mt) * R * V - gamma * (c * R - Q) * Nt * V;
    
    for i=1:m
        if abs(U(i,r) - 1) <= small_positive
            GRAD(i,r) = GRAD(i,r) - (d-1)/3;
        else
            GRAD(i,r) = GRAD(i,r) - 2*alphaU*arccosh(i)/der_arccosh(i) - (d-1)*(U(i,r)*arccosh(i)-der_arccosh(i))/(arccosh(i)*der_arccosh(i)^2);
        end
    end
end
