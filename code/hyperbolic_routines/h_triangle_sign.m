function s = h_triangle_sign(x,y,m)

%     positive if triangle inequality holds
    s = h_squared_Lorentz_distance(x,m) + h_squared_Lorentz_distance(y,m) - h_squared_Lorentz_distance(x,y);
    s = s / (h_inner(x,m) * h_inner(y,m));

end