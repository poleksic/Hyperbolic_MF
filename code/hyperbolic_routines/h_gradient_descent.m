function [U,V] = h_gradient_descent(R,U,V,DMM,DNN,c,theta,alphaU,alphaV,lM,lN,eta0,iter,cutoff)
    itr = 0;
    eta = eta0;
    last_loss = 0
    while itr < iter
        itr = itr + 1;
          first =  h_loss_function_first(R,W,U,V);
          second = h_loss_function_second(U,V,alphaU,alphaV);
          third = h_loss_function_third(U,V,S,T,betaU,betaV);
          loss = first + second + third;
        %   fprintf('first = %.5f   second = %.5f third = %.5f  score = %.5f\n',first, second,third,score);
        
        P = h_getProb(U,V,theta);

        % compute gradient wrt 
        GRAD_U = h_gradient_L_U(R, P, U, V, DMM, c, alphaU, lM);
        GRAD_V = h_gradient_L_V(R, P, U, V, DNN, c, alphaV, lN);        

        % project gradients to tangent vectors       
        GRAD_U = h_mtx_projection(U,GRAD_U);
        GRAD_V = h_mtx_projection(V,GRAD_V);

        % use Exponential Mapping to project to Hyperbolic space
        Utemp = h_mtx_exp_map(U,-1*eta*GRAD_U);
        Vtemp = h_mtx_exp_map(V,-1*eta*GRAD_V);
        if ~isNumberAndBoundedMatrix(Utemp,cutoff) || ~isNumberAndBoundedMatrix(Vtemp,cutoff)
            eta = eta/10;
            itr = itr - 1;
            continue;
        else
            U = Utemp;
            V = Vtemp;
            eta = eta0;
        end
    end
end
