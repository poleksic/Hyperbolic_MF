function proj = h_projection(mu,x)
    % orthogonal projection of vector(point) x onto tangent plane at mu
    proj = x + h_inner(mu,x)*mu;
end