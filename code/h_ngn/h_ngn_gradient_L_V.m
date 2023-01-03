function GRAD = h_ngn_gradient_L_V(A, P, U, V, DHH, c, alpha, lN)
    
    % Input: 
    % W - weight matrix, R - relationship matrix, P - matrix of probabilities
    % U, V - matrices of latent vectors, alpha - trainable parameter
    
    % Output: 
    % Matrix of gradient vectors (rows) of the loss function
    % Note: This is the true gradient in Minkowski space (as the last column is already multiplied by -1)
    small_positive = 1e-12;
    
    [n, r] = size(V);
    [m,r] = size(U);
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
    
    GRAD = P'*U + (c-1)*(A'.*P')*U - c*A'*U + lN.*(DHH*V);
        
    for j=1:n
        if abs(V(j,r) - 1) < small_positive
            GRAD(j,r) = GRAD(j,r) - (d-1)/3;
        else
            GRAD(j,r) = GRAD(j,r) - 2*alpha*arccosh(j)/der_arccosh(j) - (d-1)* (V(j,r)*arccosh(j)-der_arccosh(j))/(arccosh(j)*der_arccosh(j)^2);
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
%     HP = h_product(V,V) - small;
%     B = nN.*acosh(-1*HP)./ sqrt(HP.*HP - ones(n,n));
%      % B(i,i) = 0/0 but is equal to t_ii in the limit
%     %B(1:n+1:end) = diag(H);    
%     GRAD = GRAD - 4*lN*(B*V);
end
