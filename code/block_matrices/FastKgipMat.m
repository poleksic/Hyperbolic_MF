function [Kmat] = FastKgipMat(Y, sigma0)

%     [m,n] = size(AA);
%     sigma = sigma0 * (sum(sum(AA)) / m);
%     An = sum(AA.*AA,2); % this is a vector column
%     C = -2 * (AA * AA');
%     C = C + An;
%     C = C + An';
%     Kmat = exp(-C/sigma);
        
% P1(i,j) = exp(-a||Y(i,:) - Y(j,:)||^2) where a = (a' * nd)/sum||Y(s,:)||^2 (over all s) where nd is number of drugs
% and a' is either 1 for all i,j 9same as the algorithm above) or, better yet (as per NGN):
% a'(i,j)= 1 / (sum( Y_i_k ^ Y_j_k ) + 1)  (the sum is over all k)
% IMPORTANT: the above has to be recalculated in cross-validation
    [m,n] = size(Y);

    Kmat = zeros(m,m);
    
    sum_norms = 0;
    for i=1:m
        sum_norms = sum_norms + dot(Y(i,:),Y(i,:));
    end
    for i=1:m
        for j=1:m
%            a_prime = 1 / ( dot(Y(i,:),Y(j,:)) + 1 );
            a_prime = 1;
            a = (a_prime * m) / sum_norms;
            Kmat(i,j) = exp( -a * dot(Y(i,:)-Y(j,:), Y(i,:)-Y(j,:)) );
        end
    end
    
end

