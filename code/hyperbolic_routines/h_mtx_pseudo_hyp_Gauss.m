function [XNORM,XHYPER] = h_mtx_pseudo_hyp_Gauss(num_rows,mu,sigma)
    XNORM = zeros(num_rows,length(mu));
    XHYPER = zeros(num_rows,length(mu));
    
    for i=1:num_rows
        [XNORM(i,:), XHYPER(i,:)] = h_pseudo_hyp_Gauss(mu,sigma);
    end
end