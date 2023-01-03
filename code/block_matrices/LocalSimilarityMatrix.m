function [L] = LocalSimilarityMatrix(P, J)
% L(i,j) = P(i,j) / T  if j is one of J (J=4 by default) nearest neighbors
% Here T = Sum P(i,l) over all J nearest neighbors l.
% (All other elements of L are zeros.)

    [m,~] = size(P);
    [~, ind] = sort(P,2,'descend');
    L = zeros(m,m);
    
    for i = 1:m
        T = 0;
        for j = 1:J
            % next largest element in row i is at position col
            col = ind(i,j);
            T = T + P(i,col);
        end
        for j = 1:J
            % next largest element in row i is at position col
            col = ind(i,j);
            if T == 0
                L(i,col) = 0;
            else
                L(i,col) = P(i,col) / T;
            end
        end
    end
end