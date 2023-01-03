function [NEW_GRAD_U, NEW_GRAD_V, new_itr] = h_adjust_gradient(GRAD_U, GRAD_V, cutoff, itr)
    new_itr = itr;
    NEW_GRAD_U = GRAD_U;
    NEW_GRAD_V = GRAD_V;
    
    maxu = max(abs(GRAD_U(:)));
    maxv = max(abs(GRAD_V(:)));

    if maxu > cutoff || maxv > cutoff
%         warning('gradient too big');
        new_itr = itr - 1;
        
        [m,n] = size(GRAD_U);

        for i= 1:m
            maxu = max(abs(GRAD_U(i,:)));
            if maxu > cutoff
                NEW_GRAD_U(i,:) = (1.0/maxu)*cutoff*GRAD_U(i,:);
            end
        end

        [m,n] = size(GRAD_V);

        for i=1:m
            maxu = max(abs(GRAD_V(i,:)));
            if maxu > cutoff
                NEW_GRAD_V(i,:) = (1.0/maxu)*cutoff*GRAD_V(i,:);
            end
        end
        
    end
    

    
end