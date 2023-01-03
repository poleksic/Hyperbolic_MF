function [U, V] = h_ngn_WeightImputeLogFactorization(A,DHH,c,theta,pseudo_std,alphaU,alphaV,lM,lN,eta0,iter,d,cutoff)
    rng shuffle;
                
    [large_m, large_n] = size(A);
    mu_0 = zeros(1,d+1);
    mu_0(d+1) = 1;

    [~, U] = h_mtx_pseudo_hyp_Gauss(large_m,mu_0,pseudo_std);
    [~, V] = h_mtx_pseudo_hyp_Gauss(large_n,mu_0,pseudo_std);
    
    [U,V] = h_ngn_gradient_descent(A,U,V,DHH,c,theta,alphaU,alphaV,lM,lN,eta0,iter,cutoff);
    [U,V] = h_ngn_gradient_descent(A,U,V,DHH,c,theta,alphaU,alphaV,lM,lN,eta0/10,iter,cutoff);
    [U,V] = h_ngn_gradient_descent(A,U,V,DHH,c,theta,alphaU,alphaV,lM,lN,eta0/100,iter,cutoff);
end
