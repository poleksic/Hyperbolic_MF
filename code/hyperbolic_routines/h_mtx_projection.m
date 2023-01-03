function PROJ = h_mtx_projection(MU,X)
    % apply projection to rows
    [m,n] = size(MU);
    PROJ = zeros(m,n);
    for i=1:m
        PROJ(i,:) = h_projection(MU(i,:),X(i,:));
    end
end