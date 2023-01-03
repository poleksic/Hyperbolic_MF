function M = h_product(U,V)
   % hyperbolic product of two matrices 
   % Note: not symmetric
   [~,k] = size(U);
   [~,k] = size(V);
   
   M = U*V'-2*(U(:,k)*V(:,k)');
end