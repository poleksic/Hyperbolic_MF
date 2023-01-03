function hinner = h_inner(x,y)
    % returns Minkowski bilinear form <x,y> = x1*y1+...+xn-1*yn-1 - xn*yn
    len = length(x);
    hinner = 0;
    for i=1:len-1
        hinner = hinner + x(i)*y(i);
    end
    hinner = hinner - x(len)*y(len);
end
