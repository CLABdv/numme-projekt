function [lower, upper] = phi_src(func, lower, upper, tau)
    f1 = func(lower);
    f2 = func(upper);
    y = (sqrt(5) - 1)/2 * (upper-lower) + lower;
    x = upper - (sqrt(5) - 1)/2 * (upper - lower);
    while upper - lower > tau
        if f1 < f2
        lower = x;
        x = upper - (sqrt(5) - 1)/2 * (upper - lower);
        f1 = func(x);
    else 
        upper = y;
        y = (sqrt(5) - 1)/2 * (upper - lower) + lower;
        f2 = func(y);
    end
    end
end