function [U V] = h_dnilmf_WeightImputeLogFactorization(R,M,N,c,theta,pseudo_std,alphaU,alphaV,alpha,beta,gamma,eta0,iter,rnk,cutoff)
    rng shuffle;
        
    [m, n] = size(R);
        
    mu_0 = zeros(1,rnk+1);
    mu_0(rnk+1) = 1;

    [~, U] = h_mtx_pseudo_hyp_Gauss(m,mu_0,pseudo_std);
    [~, V] = h_mtx_pseudo_hyp_Gauss(n,mu_0,pseudo_std);

    [U,V] = h_dnilmf_gradient_descent(R,U,V,M,N,c,theta,alphaU,alphaV,alpha,beta,gamma,eta0,iter,cutoff);
    [U,V] = h_dnilmf_gradient_descent(R,U,V,M,N,c,theta,alphaU,alphaV,alpha,beta,gamma,eta0/10,iter,cutoff);
    [U,V] = h_dnilmf_gradient_descent(R,U,V,M,N,c,theta,alphaU,alphaV,alpha,beta,gamma,eta0/100,iter,cutoff);

end
