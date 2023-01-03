function [U, V] = h_nrlmf_WeightImputeLogFactorization(R,DMM,DNN,c,theta,pseudo_std,alphaU,alphaV,lM,lN,eta0,iter,d,cutoff)
    rng shuffle;
        
    [m, n] = size(R);
        
    mu_0 = zeros(1,d+1);
    mu_0(d+1) = 1;

    [~, U] = h_mtx_pseudo_hyp_Gauss(m,mu_0,pseudo_std);
    [~, V] = h_mtx_pseudo_hyp_Gauss(n,mu_0,pseudo_std);

    [U,V] = h_nrlmf_gradient_descent(R,U,V,DMM,DNN,c,theta,alphaU,alphaV,lM,lN,eta0,iter,cutoff);
    [U,V] = h_nrlmf_gradient_descent(R,U,V,DMM,DNN,c,theta,alphaU,alphaV,lM,lN,eta0/10,iter,cutoff);
    [U,V] = h_nrlmf_gradient_descent(R,U,V,DMM,DNN,c,theta,alphaU,alphaV,lM,lN,eta0/100,iter,cutoff);
end
