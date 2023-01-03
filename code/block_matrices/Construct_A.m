function [A] = Construct_A(K_d, K_t, Y, Y_d, Y_t, gamma);
    A = [K_t, gamma * Y' + (1 - gamma) * Y_t;
        gamma * Y + (1 - gamma) * Y_d, K_d'];
end