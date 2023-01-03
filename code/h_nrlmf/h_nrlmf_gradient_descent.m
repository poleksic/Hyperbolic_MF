function [U,V] = h_nrlmf_gradient_descent(R,U,V,DMM,DNN,c,theta,alphaU,alphaV,lM,lN,eta0,iter,cutoff)
    
    itr = 0;
    eta = eta0;
    [previous_first_loss, previous_second_loss, previous_third_loss] = h_nrlmf_Loss(R,U,V,DMM,DNN,c,theta,alphaU,alphaV,lM,lN);
    previous_loss = previous_first_loss + previous_second_loss + previous_third_loss;
    while itr < iter
        itr = itr + 1;
        P = h_Get_Nrlmf_Prob(U,V,theta);
        
        % compute gradient wrt 
        GRAD_U = h_nrlmf_gradient_L_U(R, P, U, V, DMM, c, alphaU, lM);
        GRAD_V = h_nrlmf_gradient_L_V(R, P, U, V, DNN, c, alphaV, lN);        
        
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
       
        [first_loss, second_loss, third_loss] = h_nrlmf_Loss(R,U,V,DMM,DNN,c,theta,alphaU,alphaV,lM,lN);
        loss = first_loss + second_loss + third_loss;
%         fprintf('first_loss = %.5f second_loss = %.5f third_loss = %.5f LOSS = %.5f\n',first_loss, second_loss,third_loss, loss);
        delta_loss = (previous_loss - loss) / abs(previous_loss);
        if abs(delta_loss) < 0.00001
            break;
            ends
        previous_loss = loss;      
    end
end
