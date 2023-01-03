function [xnorm,xhyper] = h_pseudo_hyp_Gauss(mu,sigma)
    % samples x from wrapped
    % normal distribution
    % Input 
    % mu - point on hyperboloid
    % sigma - stdev
    % Output:
    % v - Gaussian
    % x - Pseudo Hyperbolic Gaussian
   
    k = length(mu);
    u = normrnd(0,sigma,[1,k]);
    u(k) = 0;
    xnorm = h_parallel_transport(mu,u);
    xhyper = h_exp_map(mu,xnorm);
end

