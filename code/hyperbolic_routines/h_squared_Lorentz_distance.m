function dist = h_squared_Lorentz_distance(x,y)
    % returns distance of two points in hyperboloid model
    dist = h_inner(x,x) + h_inner(y,y) - 2 * h_inner(x,y);
end