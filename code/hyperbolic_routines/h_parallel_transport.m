function y = h_parallel_transport(mu,x)
    % y in tangent space at mu is obtained by parallel transport of x along
    % the geodesic from mu_0 to mu

    l = length(mu);
    mu_0 = zeros(1,l);
    mu_0(l) = 1;
    
    numer = h_inner(mu + h_inner(mu_0,mu) * mu_0, x);
    denom = 1-h_inner(mu_0,mu);
    
    y = x + (numer/denom)*(mu_0 + mu);
end

