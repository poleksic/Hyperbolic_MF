function [first_loss,second_loss] = h_dnilmf_Loss(R,U,V,M,N,c,alpha,beta,gamma,theta,alphaU,alphaV)
    [m,n] = size(R);
    [m,r] = size(U);
    DOT = h_product(U,V);
    MDOT = M * DOT;
    DOTN = DOT * N;
    CORE = alpha * DOT + beta * MDOT + gamma * DOTN + theta;
    LL = (1 + c*R -R) .* log(1 + exp(CORE)) - c*R .* CORE;
    first_loss = sum(sum(LL));
    
    d = r-1;

    second_loss = 0;
    for i=1:m
       total = dot(U(i,:),U(i,:)) - U(i,r)^2;
       second_loss = second_loss + alphaU * total * acosh(U(i,r))^2 / (U(i,r)^2 - 1) + (d-1) * log(sqrt(U(i,r)^2 - 1)/acosh(U(i,r)));
    end
    for j=1:n
       total = dot(V(j,:),V(j,:)) - V(j,r)^2;
       second_loss = second_loss + alphaV * total * acosh(V(j,r))^2 / (V(j,r)^2 - 1) + (d-1) * log(sqrt(V(j,r)^2 - 1)/acosh(V(j,r)));
    end
end

