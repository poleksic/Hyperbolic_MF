% Form matrix consisting of ExcludedRows of U; 
% among all matrix columns, find one with the smallest range of
% values of U and return the index "ind" of that column
function [mk,ind] = h_find_smallest_variance_coord(U, ExcludedRows)
    U = U(ExcludedRows,:);
    K = max(U) - min(U);
    [mk, ind] = min(K);
end
