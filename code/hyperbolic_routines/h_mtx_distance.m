function DIST = h_mtx_distance(X,Y)
    DIST = acosh(-1 * h_product(X,Y));
end