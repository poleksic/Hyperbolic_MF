function [N] = NormalizeRows(P)
    
%     This is how DNILMF (demo) does it, but it's probably wriong
%     P = P ./ sum(P,2);
%     P(isnan(P)) = 0;
    
% NEW WAY (USED??):
% If i ~= j then N(i,j) = P(i,j) / (2*T) where T = Sum P(i,j) over all j!=i
% If i = j then N(i,j) = 1/2

    [m,~] = size(P);
    N = zeros(m,m);

    for i = 1:m
        E = 2 * (sum(P(i,:)) - P(i,i));
        for j = 1:m
            if i ~= j
                if E ~= 0
                    N(i,j) = P(i,j) / E;
                end
            else
                N(i,j) = 1/2;
            end
        end
    end
 
end