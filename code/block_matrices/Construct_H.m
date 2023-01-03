function [H] = Construct_H(K_d, K_t, Y);
    K_d_hat = logical(K_d);
    K_t_hat = logical(K_t);
    Y_hat = Y ./ sum(Y);
    Y_hat(isnan(Y_hat)) = 0;
    H = [K_t_hat, Y_hat';
        Y_hat, K_d_hat];
end