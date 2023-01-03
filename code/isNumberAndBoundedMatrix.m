function a = isNumberAndBoundedMatrix(X,bound)
   a = isNumberMatrix(X);
   if max(abs(X(:))) > bound
       a = 0;
   end
end