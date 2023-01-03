%Fermi-Dirac distribution
function L = h_loss_function_first(R,W,U,V,theta,t)
   [m,n] = size(R);
   
   L = 0;
   for i=1:m
       for j=1:n
          L = L + W(i,j) * (log(1+exp((h_inner(U(i,:),V(j,:)+theta)/t))) - R(i,j) * (h_inner(U(i,:),V(j,:)+theta)/t));
       end
   end
end