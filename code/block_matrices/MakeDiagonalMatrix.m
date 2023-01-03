function [G, G_bar] = MakeDiagonalMatrix(H)
    % As per the NGN routine
    [m,~] = size(H);
    D1 = sum(H,2);
    D2 = sum(H,1);
    G = zeros(m,m);
    G_bar = zeros(m);
    for i=1:m
        G(i,i) = D1(i);
        G_bar(i,i) = D2(i);
    end
end

