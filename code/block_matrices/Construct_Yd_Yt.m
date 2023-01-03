function [Y_dt] = Construct_Yd_Yt(Y,K_d,K_t)
    [m,~] = size(K_d);
    [n,~] = size(K_t);
    % This routine builds a global interaction matrix
    Zi_Zj = logical(K_d) * logical(Y) * logical(K_t');
    Y_dt = (K_d * Y * K_t') ./ Zi_Zj;
    Y_dt(isnan(Y_dt)) = 0;
end
