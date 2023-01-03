function [U,V] = h_dnilmf_gradient_descent(R,U,V,M,N,c,theta,alphaU,alphaV,alpha,beta,gamma,eta0,iter,cutoff)
    
    itr = 0;
    eta = eta0;
    [previous_first_loss, previous_second_loss] =  h_dnilmf_Loss(R,U,V,M,N,c,alpha,beta,gamma,theta,alphaU,alphaV);
    previous_loss = previous_first_loss + previous_second_loss;
    while itr < iter
        itr = itr + 1;
             
        % compute gradient wrt 
        GRAD_U = h_dnilmf_gradient_L_U(R, U, V, M, N, c, alphaU, alpha, beta, gamma, theta);
        GRAD_V = h_dnilmf_gradient_L_V(R, U, V, M, N, c, alphaV, alpha, beta, gamma, theta);        

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
        [first_loss, second_loss] =  h_dnilmf_Loss(R,U,V,M,N,c,alpha,beta,gamma,theta,alphaU,alphaV);
        loss = first_loss + second_loss;
%         fprintf('first_loss = %.5f   second_loss = %.5f LOSS = %.5f\n',first_loss, second_loss,loss);
        delta_loss = (previous_loss - loss) / abs(previous_loss);
        if abs(delta_loss) < 0.00001
            break;
        end
        previous_loss = loss;
    end
end
