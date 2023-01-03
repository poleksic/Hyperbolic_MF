function [S] = MakePositiveDefinite(S)
% Check if S is positive definite
% If not, S = S + small * I, where small is a small positive value

    [m,~] = size(S);
    small = 0.01;
    [~,p] = chol(S);
    i=0;
    while (p > 0)
        [~,p] = chol(S);
        S = S + small * eye(m);
    end
end

