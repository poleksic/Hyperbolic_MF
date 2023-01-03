function [Y_dt] = Construct_Y_dt(Y,K_dt)
    Zi_Zj = logical(K_dt) * logical(Y);
    Y_dt = (K_dt * Y) ./ Zi_Zj;
    Y_dt(isnan(Y_dt)) = 0;
end
