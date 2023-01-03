function Y_new = AmmendProfile(Y, S, J)
    [m,~] = size(Y);
    Y_new = Y;
    [indexZeros,~] = find(sum(Y,2)==0);
    numIndexZeros = sum(logical(indexZeros));

    SR = S;
    
    SR = SR - diag(diag(SR));
    SR(:,indexZeros) = 0;

    [~,ind] = sort(SR,2,'descend');
    TOP_HORIZ = zeros(m,m);
    for i = 1:m
        for j = 1:J
            col = ind(i,j);
            TOP_HORIZ(i,col) = SR(i,col);
        end
    end
               
    for x=1:numIndexZeros
        i = indexZeros(x);
        Y_new(i,:) = TOP_HORIZ(i,:) * Y ./ sum(TOP_HORIZ(i,:));
    end
    
end


