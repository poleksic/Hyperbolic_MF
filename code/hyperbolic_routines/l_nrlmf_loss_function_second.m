function L = h_loss_function_second(U,V,alphaU,alphaV)
   [m,r] = size(U);
   [n,r] = size(V);
   d = r-1;
   
   L = 0;
   for i=1:m
       sum = dot(U(i,:),U(i,:)) - U(i,r)^2;
       L = L + alphaU * sum * acosh(U(i,r))^2 / (U(i,r)^2 - 1) + (d-1) * log(sqrt(U(i,r)^2 - 1)/acosh(U(i,r)));
   end
   for j=1:n
       sum = dot(V(j,:),V(j,:)) - V(j,r)^2;
       L = L + alphaV * sum * acosh(V(j,r))^2 / (V(j,r)^2 - 1) + (d-1) * log(sqrt(V(j,r)^2 - 1)/acosh(V(j,r)));
   end
end