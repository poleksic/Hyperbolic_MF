function a = isNumberMatrix(X)
   sm = sum(sum(~isreal(X))) + sum(sum(isnan(X))) + sum(sum(isinf(X))); 
   a = 1 - logical(sm);
end