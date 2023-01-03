function dist = h_distance(x,y)
    % returns distance of two points in hyperboloid model
    dist = acosh(-1 * h_inner(x,y));
end