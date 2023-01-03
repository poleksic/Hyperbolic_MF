function [K] = NeighborSimMtx(S,J)
    [m,~] = size(S);
    K = zeros(m,m);

    S(1:m+1:m*m) = 0;
    % sort rows in descending order
    [~, ind] = sort(S,2,'descend');
    
    for i = 1:m
        for j = 1:J
            % next largest element in row i is at position col
            col = ind(i,j);
            K(i,col) = S(i,col);
        end
    end
end
