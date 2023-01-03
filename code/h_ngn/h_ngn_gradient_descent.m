function [U,V] = h_ngn_gradient_descent(A,U,V,DHH,c,theta,alphaU,alphaV,lM,lN,eta0,iter,cutoff)
    
    itr = 0;
    eta = eta0;
    [previous_first_loss, previous_second_loss, previous_third_loss] = h_ngn_Loss(A,U,V,DHH,c,theta,alphaU,alphaV,lM,lN);
    previous_loss = previous_first_loss + previous_second_loss + previous_third_loss;
    while itr < iter
        itr = itr + 1;

        P = h_Get_Nrlmf_Prob(U,V,theta);
        % compute gradient wrt 
        GRAD_U = h_ngn_gradient_L_U(A, P, U, V, DHH, c, alphaU, lM);
        GRAD_V = h_ngn_gradient_L_V(A, P, U, V, DHH, c, alphaV, lN);      
        
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
       
        [first_loss, second_loss, third_loss] = h_ngn_Loss(A,U,V,DHH,c,theta,alphaU,alphaV,lM,lN);
        loss = first_loss + second_loss + third_loss;
%         fprintf('itr = %d first_loss = %.5f second_loss = %.5f third_loss = %.5f LOSS = %.5f\n',itr,first_loss, second_loss, third_loss, loss);
        delta_loss = (previous_loss - loss) / abs(previous_loss);
        if abs(delta_loss) < 0.00001
            break;
        end
        previous_loss = loss;      
    end
end
