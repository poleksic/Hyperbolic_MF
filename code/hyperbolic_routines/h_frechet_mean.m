function y = h_frechet_mean(X,w,iter)
    % Input: 
    % t points in n-dimensional hyperbolic space (rows of X)
    % and weights (given in the row) w
    
    % Output: Frechet mean y (as a row) of points X 
    
    [t,~] = size(X);
    m = length(w);
    
    if t~=m
        error('dimensions error in Frechet mean');
    end
    
    y = X(1,:);
    for k=1:iter
        u = next_u(X,w,y);
        y = (1.0/sqrt(-1*h_inner(u,u)))*u;
    end
    
end

function u = next_u(X,w,y)
    n = length(y);
    t = length(w);
    u = zeros(1,n);
    for l=1:t
        if abs(h_inner(X(l,:),y) + 1) < 0.00000001
            second = 1;
        else
            second = acosh(-1*h_inner(X(l,:),y))/sqrt((h_inner(X(l,:),y))^2 - 1);
        end
        u = u + 2*w(l)*second*X(l,:);
    end
end