function L = h_loss_function_third(U,V,S,T,betaU,betaV)
   [m,r] = size(U);
   [n,r] = size(V);
   d = r-1;
   
   LU = 0;
   for i=1:m
       for j=1:m
           LU = LU + S(i,j)*h_inner(U(i,:)-U(j,:),U(i,:)-U(j,:));
       end
   end
   LU = betaU * LU;
   LV = 0;
   for i=1:n
       for j=1:n
            LV = LV + T(i,j)*h_inner(V(i,:)-V(j,:),V(i,:)-V(j,:));
       end
   end
   LV = betaV * LV;
   L = LU + LV;
end