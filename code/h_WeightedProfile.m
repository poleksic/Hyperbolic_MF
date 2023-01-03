function nU = h_WeightedProfile(U, S, row_no, J)
    nU = U;
    [m,r] = size(U);
    for i = 1:length(row_no)
            A = S(row_no(i),:);
            A(row_no) = 0;
            [values, ind] = sort(A, 'descend');
             values = values(1:J);
             sum_values = sum(values);
             ind = ind(1:J);
             if  (sum_values == 0)
                 warning('There are no nearest neighbors for first latent matrix')
             end
             nU(row_no(i),:) = h_frechet_mean(U(ind,:),values,20);
    end
end
