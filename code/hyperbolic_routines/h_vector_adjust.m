function new_x = h_vector_adjust(x)
    % adjust slight deviation of x to place it on hyperboloid
    new_x = (sqrt(-1/h_inner(x,x)))*x;
end