function cosh = h_cosh_angle(a,b)

    cosh = -h_inner(a,b) / (h_norm(a)*h_norm(b));
end