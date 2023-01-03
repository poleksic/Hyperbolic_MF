function GRAD = h_ngn_gradient_L_U(A, P, U, V, DHH, c, alpha, lM)
    
    % Input: 
    % W - weight matrix, R - relationship matrix, P - matrix of probabilities
    % U, V - matrices of latent vectors, alpha - trainable parameter
    
    % Output: 
    % Matrix of gradient vectors (rows) of the loss function
    % Note: This is the true gradient in Minkowski space (as the last column is already multiplied by -1)
    
    small_positive = 1e-12;
    
    [m,n] = size(A);
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
    
    GRAD = P*V + (c-1)*(A.*P)*V - c*A*V + lM.*(DHH*U);
        
    for i=1:m
        if abs(U(i,r) - 1) <= small_positive
            GRAD(i,r) = GRAD(i,r) - (d-1)/3;
        else
            GRAD(i,r) = GRAD(i,r) - 2*alpha*arccosh(i)/der_arccosh(i) - (d-1)*(U(i,r)*arccosh(i)-der_arccosh(i))/(arccosh(i)*der_arccosh(i)^2);
        end
    end


    % linear
%    GRAD = GRAD - 4*beta*(H*U);

    % Lorentzian regularization (this is identical to Euclidean)
%     GRAD = GRAD + 4*beta*(U .* (H*ones(m,r)) - H*U);

    % distance regularization (use this)
    % adding to the loss function the sum of the terms s_ij*d^2(Ui,Uj)
    % where d is the distance (not Lorentzian but ordinary) between the 
    % points in the hyperbolic space (d(x,y) = arccosh(-<x,y>))
%     HP = h_product(U,U) - small;
%     B = nM.*acosh(-1*HP)./ sqrt(HP.*HP - ones(m,m));
% %      B(i,i) = 0/0 but is equal to s_ii in the limit
% %     B(1:m+1:end) = diag(H);    
%     GRAD = GRAD - 4*lM*(B*U);

end
