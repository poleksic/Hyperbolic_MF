% x_h = x ./ sqrt(-1 * h_inner(x,x))

function x_h = h_adjust_on_hyperboloid(x)

    if (h_inner(x,x) > 0)
        fprintf('inner product = %.2f\n',h_inner(x,x));
        error('h_adjust_on_hyperboloid: inner product positive')
    else
        x_h = x ./ sqrt(-1 * h_inner(x,x));
    end

end