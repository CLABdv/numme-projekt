function X = gaussNewton(X0, F, J, tau)
% F är en funktion och J är dess Jacobian, tau är felgräns

    while true
        dx = - J(X0) \ F(X0);
        X0 = X0 + dx;
        if all (abs(dx) - tau < 0)
            break
        end
    end

    X = X0;
end