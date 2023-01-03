function [S] = Symmetrize(S)
% make the input matrix symmetric
    S = (S + S')/2; 
end

