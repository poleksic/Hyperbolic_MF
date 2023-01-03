function [S] = Jaccard_Mtx(X,Y)

    [m,n] = size(X);
    [p,k] = size(Y);
    if n~=k 
        error('Error in Jaccard_Mtx.')
    end
    
    S = zeros(m,p);
    for i=1:m
        for j=1:p
            S(i,j) = Jaccard(X(i,:),Y(j,:));
        end
    end

end

