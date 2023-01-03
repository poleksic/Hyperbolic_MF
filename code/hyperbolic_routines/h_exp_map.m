function exp_mu = h_exp_map(mu,x)
    % returns the value of exponential map:  exp_mu = Exp_mu(x)
    % this value has the property hdist(mu,exp_mu)=norm(x)
    norm_x = h_norm(x);
    exp_mu = cosh(norm_x)*mu + (sinh(norm_x)/norm_x)*x;
    exp_mu = h_vector_check_adjust(exp_mu);
end