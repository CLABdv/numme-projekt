% f är funktionen som används för att finna lägsta relativa ljudnivån för
% TV'n
function A = f(xs, ys, a, S0, n)
    omega = 30;
    S = @(x, y) a * S0(x-xs, y-ys);
    [~, Sol] = hhsolver(omega, S, n);
    w = find(Sol.x<=0.25 & Sol.y>=0.5);
    A = max(abs(Sol.u(w)))/max(abs(Sol.u(:)));
end