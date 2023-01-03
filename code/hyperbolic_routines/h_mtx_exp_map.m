function EXP_MU = h_mtx_exp_map(MU,X)
    % apply exponent map to rows
    [m,n] = size(MU);
    EXP_MU = zeros(m,n);
    for i=1:m
        EXP_MU(i,:) = h_exp_map(MU(i,:),X(i,:));
    end
end